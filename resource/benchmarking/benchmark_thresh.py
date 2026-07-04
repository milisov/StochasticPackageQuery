import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

base_dir = "/home/fm2288/StochasticPackageQuery/test/HyperParameterBenchmark3RSRelaxedObjCons/"
output_pdf_path = os.path.join(base_dir, "threshold_search.pdf")

result_df = pd.DataFrame(columns=['Threshold', 'M',
                                  'NumberFeasibleMean', 'NumberFeasibleStd',
                                  'ObjRatioMean', 'ObjRatioStd',
                                  'RuntimeMean', 'RuntimeStd'])

for thresh_dir in os.listdir(base_dir):
    thresh_path = os.path.join(base_dir, thresh_dir)
    if not os.path.isdir(thresh_path):
        continue

    threshold = thresh_dir.split("_")[1]

    tablesByM = {}
    for table in os.listdir(thresh_path):
        parts = table.split(".")[0].split("_")
        N = parts[1]
        if N == "5":
            continue
        M = parts[2]
        df = pd.read_csv(os.path.join(thresh_path, table))
        df['Table'] = table.split(".")[0]
        if M not in tablesByM:
            tablesByM[M] = []
        tablesByM[M].append(df)

    for M, dfs in tablesByM.items():
        df = pd.concat(dfs, ignore_index=True)
        feasible_mask = (df['deter_feas'] == 1) & (df['prob_feas'] == 1)

        numberOfFeasibleForTable = df.groupby(['Hardness', 'Table']).apply(lambda x: feasible_mask[x.index].sum())
        print(numberOfFeasibleForTable)
        numberOfFeasibleByHardness = numberOfFeasibleForTable.reset_index(name='numberOfFeasible')
        numFeasibleMean = round(numberOfFeasibleByHardness['numberOfFeasible'].mean(), 4)
        numFeasibleStd = round(numberOfFeasibleByHardness['numberOfFeasible'].std(), 4)

        objRatioMean = round(df[feasible_mask]['ObjRatio'].mean(), 4)
        objRatioStd = round(df[feasible_mask]['ObjRatio'].std(), 4)

        runtimeMean = round(df['Runtime'].mean(), 4)
        runtimeStd = round(df['Runtime'].std(), 4)

        result_df = pd.concat([result_df, pd.DataFrame([{
            'Threshold': threshold,
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

feasible_pivot = result_df.pivot_table(index='Threshold', columns='M', values='NumberFeasible', aggfunc='first').reset_index()
feasible_pivot.columns.name = None
feasible_pivot.rename(columns={m: f'NumberFeasible_{m}' for m in ['10', '100', '10000']}, inplace=True)

objratio_pivot = result_df.pivot_table(index='Threshold', columns='M', values='ObjRatio', aggfunc='first').reset_index()
objratio_pivot.columns.name = None
objratio_pivot.rename(columns={m: f'ObjRatio_{m}' for m in ['10', '100', '10000']}, inplace=True)

runtime_pivot = result_df.pivot_table(index='Threshold', columns='M', values='Runtime', aggfunc='first').reset_index()
runtime_pivot.columns.name = None
runtime_pivot.rename(columns={m: f'Runtime_{m}' for m in ['10', '100', '10000']}, inplace=True)

display_df = feasible_pivot.merge(objratio_pivot, on='Threshold').merge(runtime_pivot, on='Threshold')

sum_feasible = result_df.groupby('Threshold')['NumberFeasibleMean'].sum().reset_index(name='SumNumberFeasibleMean')
sum_obj = result_df.groupby('Threshold')['ObjRatioMean'].sum().reset_index(name='SumObjRatioMean')
display_df = display_df.merge(sum_feasible, on='Threshold').merge(sum_obj, on='Threshold')
display_df['Threshold'] = display_df['Threshold'].astype(float)
display_df.sort_values(by=['SumNumberFeasibleMean','SumObjRatioMean'], ascending=[False, False], inplace=True)
# display_df.drop(columns=['SumNumberFeasibleMean', 'SumObjRatioMean'], inplace=True)

print(display_df.to_string(index=False))

fig, ax = plt.subplots(figsize=(20, len(display_df) * 0.5 + 1))
ax.axis('off')
tbl = ax.table(cellText=display_df.values, colLabels=display_df.columns, loc='center', cellLoc='center')
tbl.auto_set_font_size(False)
tbl.set_fontsize(7)
tbl.auto_set_column_width(col=list(range(len(display_df.columns))))
plt.tight_layout()
plt.savefig(output_pdf_path, bbox_inches='tight')
