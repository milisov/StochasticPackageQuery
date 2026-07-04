import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

ablation_dir = "/home/fm2288/StochasticPackageQuery/test/AblationStudyTerminationCondition"
rs_mad_dir = "/home/fm2288/StochasticPackageQuery/test/RS_MAD_upto_1024"
output_pdf_path = os.path.join(ablation_dir, "ablation_termination_condition.pdf")

stages = [
    ("stage1lcvar", os.path.join(ablation_dir, "stage1lcvar")),
    ("stage1random", os.path.join(ablation_dir, "stage1random")),
    ("RS_MAD", rs_mad_dir),
]


def analyze_stage(stage_dir):
    rows = []

    for thresh_dir in sorted(os.listdir(stage_dir)):
        thresh_path = os.path.join(stage_dir, thresh_dir)
        if not os.path.isdir(thresh_path):
            continue

        threshold = thresh_dir.split("_")[1]

        tablesByM = {}
        for table in os.listdir(thresh_path):
            if not table.endswith(".csv"):
                continue
            parts = table.split(".")[0].split("_")
            N = parts[1]
            M = parts[2]
            df = pd.read_csv(os.path.join(thresh_path, table))
            df['Table'] = table.split(".")[0]
            if M not in tablesByM:
                tablesByM[M] = []
            tablesByM[M].append(df)

        for M, dfs in tablesByM.items():
            df = pd.concat(dfs, ignore_index=True)
            feasible_mask = (df['deter_feas'] == 1) & (df['prob_feas'] == 1)

            numberOfFeasibleForTable = df.groupby(['Hardness', 'Table']).apply(
                lambda x: feasible_mask[x.index].sum()
            )
            numberOfFeasibleByHardness = numberOfFeasibleForTable.reset_index(name='numberOfFeasible')
            numFeasibleMean = round(numberOfFeasibleByHardness['numberOfFeasible'].mean(), 4)
            numFeasibleStd = round(numberOfFeasibleByHardness['numberOfFeasible'].std(), 4)

            objRatioMean = round(df[feasible_mask]['ObjRatio'].mean(), 4)
            objRatioStd = round(df[feasible_mask]['ObjRatio'].std(), 4)

            runtimeMean = round(df['Runtime'].mean(), 4)
            runtimeStd = round(df['Runtime'].std(), 4)

            ntuplesMean = round(df['NTuples'].mean(), 4)

            rows.append({
                'Threshold': threshold,
                'M': M,
                'NumberFeasibleMean': numFeasibleMean,
                'NumberFeasibleStd': numFeasibleStd,
                'ObjRatioMean': objRatioMean,
                'ObjRatioStd': objRatioStd,
                'RuntimeMean': runtimeMean,
                'RuntimeStd': runtimeStd,
                'NTuplesMean': ntuplesMean,
            })

    result_df = pd.DataFrame(rows)

    result_df['NumberFeasible'] = result_df['NumberFeasibleMean'].astype(str) + ' ± ' + result_df['NumberFeasibleStd'].astype(str)
    result_df['ObjRatio'] = result_df['ObjRatioMean'].astype(str) + ' ± ' + result_df['ObjRatioStd'].astype(str)
    result_df['Runtime'] = result_df['RuntimeMean'].astype(str) + ' ± ' + result_df['RuntimeStd'].astype(str)

    m_values = sorted(result_df['M'].unique(), key=lambda x: int(x))

    feasible_pivot = result_df.pivot_table(index='Threshold', columns='M', values='NumberFeasible', aggfunc='first').reset_index()
    feasible_pivot.columns.name = None
    feasible_pivot.rename(columns={m: f'NumberFeasible_{m}' for m in m_values}, inplace=True)

    objratio_pivot = result_df.pivot_table(index='Threshold', columns='M', values='ObjRatio', aggfunc='first').reset_index()
    objratio_pivot.columns.name = None
    objratio_pivot.rename(columns={m: f'ObjRatio_{m}' for m in m_values}, inplace=True)

    runtime_pivot = result_df.pivot_table(index='Threshold', columns='M', values='Runtime', aggfunc='first').reset_index()
    runtime_pivot.columns.name = None
    runtime_pivot.rename(columns={m: f'Runtime_{m}' for m in m_values}, inplace=True)

    ntuples_pivot = result_df.pivot_table(index='Threshold', columns='M', values='NTuplesMean', aggfunc='first').reset_index()
    ntuples_pivot.columns.name = None
    ntuples_pivot.rename(columns={m: f'NTuplesMean_{m}' for m in m_values}, inplace=True)

    display_df = (feasible_pivot
                  .merge(objratio_pivot, on='Threshold')
                  .merge(runtime_pivot, on='Threshold')
                  .merge(ntuples_pivot, on='Threshold'))

    sum_feasible = result_df.groupby('Threshold')['NumberFeasibleMean'].sum().reset_index(name='SumNumberFeasibleMean')
    sum_obj = result_df.groupby('Threshold')['ObjRatioMean'].sum().reset_index(name='SumObjRatioMean')
    display_df = display_df.merge(sum_feasible, on='Threshold').merge(sum_obj, on='Threshold')
    display_df['Threshold'] = display_df['Threshold'].astype(float)
    display_df.sort_values(by=['SumObjRatioMean', 'SumNumberFeasibleMean'], ascending=False, inplace=True)

    return display_df


with PdfPages(output_pdf_path) as pdf:
    for stage_name, stage_dir in stages:
        display_df = analyze_stage(stage_dir)
        print(f"\n=== {stage_name} ===")
        print(display_df.to_string(index=False))

        fig, ax = plt.subplots(figsize=(24, max(len(display_df) * 0.5 + 1.5, 3)))
        ax.axis('off')
        ax.set_title(stage_name, fontsize=12, fontweight='bold', pad=15)
        tbl = ax.table(cellText=display_df.values, colLabels=display_df.columns, loc='center', cellLoc='center')
        tbl.auto_set_font_size(False)
        tbl.set_fontsize(7)
        tbl.auto_set_column_width(col=list(range(len(display_df.columns))))
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()

print(f"\nPDF saved to {output_pdf_path}")
