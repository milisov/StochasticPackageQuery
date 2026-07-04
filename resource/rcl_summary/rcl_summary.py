import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

base_dir = "/home/fm2288/StochasticPackageQuery/stochastic-sketchrefine/results/ImputedDataToleranceApproximationBound/"
output_pdf_path = os.path.join(os.path.dirname(__file__), "rcl_summary.pdf")

result_df = pd.DataFrame(columns=['M',
                                  'NumberFeasibleMean', 'NumberFeasibleStd',
                                  'ObjRatioMean', 'ObjRatioStd',
                                  'RuntimeMean', 'RuntimeStd'])

tablesByM = {}
for table in os.listdir(base_dir):
    if not table.endswith(".csv"):
        continue
    parts = table.split(".")[0].split("_")
    N = parts[2]
    if N == "5":
        continue
    M = parts[3]
    df = pd.read_csv(os.path.join(base_dir, table))
    df['Table'] = table.split(".")[0]
    if M not in tablesByM:
        tablesByM[M] = []
    tablesByM[M].append(df)

for M, dfs in tablesByM.items():
    df = pd.concat(dfs, ignore_index=True)
    feasible_mask = (df['deter_feas'] == 1) & (df['prob_feas'] == 1)

    numberOfFeasibleForTable = df.groupby(['Hardness', 'Table']).apply(lambda x: feasible_mask[x.index].sum())
    numberOfFeasibleByHardness = numberOfFeasibleForTable.reset_index(name='numberOfFeasible')
    numFeasibleMean = round(numberOfFeasibleByHardness['numberOfFeasible'].mean(), 4)
    numFeasibleStd = round(numberOfFeasibleByHardness['numberOfFeasible'].std(), 4)

    objRatioMean = round(df[feasible_mask]['ObjRatio'].mean(), 4)
    objRatioStd = round(df[feasible_mask]['ObjRatio'].std(), 4)

    runtimeMean = round(df['Runtime'].mean(), 4)
    runtimeStd = round(df['Runtime'].std(), 4)

    result_df = pd.concat([result_df, pd.DataFrame([{
        'M': M,
        'NumberFeasibleMean': numFeasibleMean,
        'NumberFeasibleStd': numFeasibleStd,
        'ObjRatioMean': objRatioMean,
        'ObjRatioStd': objRatioStd,
        'RuntimeMean': runtimeMean,
        'RuntimeStd': runtimeStd,
    }])], ignore_index=True)

result_df['NumberFeasible'] = result_df['NumberFeasibleMean'].astype(str) + ' ± ' + result_df['NumberFeasibleStd'].astype(str)
result_df['ObjRatio'] = result_df['ObjRatioMean'].astype(str) + ' ± ' + result_df['ObjRatioStd'].astype(str)
result_df['Runtime'] = result_df['RuntimeMean'].astype(str) + ' ± ' + result_df['RuntimeStd'].astype(str)

feasible_pivot = result_df.pivot_table(index=[], columns='M', values='NumberFeasible', aggfunc='first').reset_index(drop=True)
feasible_pivot.columns.name = None
feasible_pivot.rename(columns={m: f'NumberFeasible_{m}' for m in ['10', '100', '10000']}, inplace=True)

objratio_pivot = result_df.pivot_table(index=[], columns='M', values='ObjRatio', aggfunc='first').reset_index(drop=True)
objratio_pivot.columns.name = None
objratio_pivot.rename(columns={m: f'ObjRatio_{m}' for m in ['10', '100', '10000']}, inplace=True)

runtime_pivot = result_df.pivot_table(index=[], columns='M', values='Runtime', aggfunc='first').reset_index(drop=True)
runtime_pivot.columns.name = None
runtime_pivot.rename(columns={m: f'Runtime_{m}' for m in ['10', '100', '10000']}, inplace=True)

display_df = pd.concat([feasible_pivot, objratio_pivot, runtime_pivot], axis=1)

print(display_df.to_string(index=False))

fig, ax = plt.subplots(figsize=(20, len(display_df) * 0.5 + 1))
ax.axis('off')
tbl = ax.table(cellText=display_df.values, colLabels=display_df.columns, loc='center', cellLoc='center')
tbl.auto_set_font_size(False)
tbl.set_fontsize(7)
tbl.auto_set_column_width(col=list(range(len(display_df.columns))))
plt.tight_layout()
plt.savefig(output_pdf_path, bbox_inches='tight')
