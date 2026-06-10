import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

dir = "/home/fm2288/StochasticPackageQuery/test/HyperParameterBenchmark1/"

outputdir = "/home/fm2288/StochasticPackageQuery/test/HyperParameterBenchmark1/benchmark_summary_scenarios_cap.pdf"

result_df = pd.DataFrame(columns=['Tuples', 'Scenarios', 'Cap', 'M',
                                  'NumberFeasibleMean', 'NumberFeasibleStd',
                                  'ObjRatioMean', 'ObjRatioStd',
                                  'RuntimeMean', 'RuntimeStd'])


for hyperparam in os.listdir(dir):
    tuples = hyperparam.split("_")[1]
    scenarios = hyperparam.split("_")[3]
    cap = hyperparam.split("_")[5]

    tablesByM = {}
    for table in os.listdir(os.path.join(dir, hyperparam)):
        M = table.split(".")[0].split("_")[2]
        df = pd.read_csv(os.path.join(dir, hyperparam, table))
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
            'Tuples': tuples,
            'Scenarios': scenarios,
            'Cap': cap,
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

feasible_pivot = result_df.pivot_table(index=['Tuples', 'Scenarios', 'Cap'], columns='M', values='NumberFeasible', aggfunc='first').reset_index()
feasible_pivot.columns.name = None
feasible_pivot.rename(columns={m: f'NumberFeasible_{m}' for m in ['10', '100', '10000']}, inplace=True)

objratio_pivot = result_df.pivot_table(index=['Tuples', 'Scenarios', 'Cap'], columns='M', values='ObjRatio', aggfunc='first').reset_index()
objratio_pivot.columns.name = None
objratio_pivot.rename(columns={m: f'ObjRatio_{m}' for m in ['10', '100', '10000']}, inplace=True)

display_df = feasible_pivot.merge(objratio_pivot, on=['Tuples', 'Scenarios', 'Cap'])

sum_feasible = result_df.groupby(['Tuples', 'Scenarios', 'Cap'])['NumberFeasibleMean'].sum().reset_index(name='SumNumberFeasibleMean')
sum_obj = result_df.groupby(['Tuples', 'Scenarios', 'Cap'])['ObjRatioMean'].sum().reset_index(name='SumObjRatioMean')
display_df = display_df.merge(sum_feasible, on=['Tuples', 'Scenarios', 'Cap'])
display_df = display_df.merge(sum_obj, on=['Tuples', 'Scenarios', 'Cap'])
display_df.sort_values(by=['SumNumberFeasibleMean', 'SumObjRatioMean'], ascending=False, inplace=True)
# display_df.drop(columns=['SumNumberFeasibleMean'], inplace=True)

print(display_df)

fig, ax = plt.subplots(figsize=(20, len(display_df) * 0.5 + 1))
ax.axis('off')
tbl = ax.table(cellText=display_df.values, colLabels=display_df.columns, loc='center', cellLoc='center')
tbl.auto_set_font_size(False)
tbl.set_fontsize(7)
tbl.auto_set_column_width(col=list(range(len(display_df.columns))))
plt.tight_layout()
plt.savefig(outputdir, bbox_inches='tight')
