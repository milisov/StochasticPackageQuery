import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

dir = "/home/fm2288/StochasticPackageQuery/test/HyperParameterBenchmark2/"

output_csv_path = "/home/fm2288/StochasticPackageQuery/test/HyperParameterBenchmark2/benchmark_secondary_summary_tuples.csv"

output_pdf_path = "/home/fm2288/StochasticPackageQuery/test/HyperParameterBenchmark2/benchmark_secondary_summary_tuples.pdf"

result_df = pd.DataFrame(columns=['Tuples', 'Scenarios', 'Cap',
                                  'NumberFeasibleMean', 'NumberFeasibleStd',
                                  'ObjRatioMean', 'ObjRatioStd',
                                  'RuntimeMean', 'RuntimeStd'])


for hyperparam in os.listdir(dir):
    tuples = hyperparam.split("_")[1]
    scenarios = hyperparam.split("_")[3]
    cap = hyperparam.split("_")[5]

    all_dfs = []
    for table in os.listdir(os.path.join(dir, hyperparam)):
        df = pd.read_csv(os.path.join(dir, hyperparam, table))
        df['Table'] = table.split(".")[0]
        all_dfs.append(df)

    df = pd.concat(all_dfs, ignore_index=True)
    feasible_mask = (df['deter_feas'] == 1) & (df['prob_feas'] == 1)

    numberOfFeasibleForTable = df.groupby(['Hardness', 'Table']).apply(lambda x: feasible_mask[x.index].sum())
    numberOfFeasibleByHardness = numberOfFeasibleForTable.reset_index(name='numberOfFeasible')
    numFeasibleMean = round(numberOfFeasibleByHardness['numberOfFeasible'].mean(), 4)
    numFeasibleStd = round(numberOfFeasibleByHardness['numberOfFeasible'].std(), 4)

    has_solution_mask = (df['deter_feas'] != 0) | (df['prob_feas'] != 0)
    objRatioMean = round(df[has_solution_mask]['ObjRatio'].mean(), 4)
    objRatioStd = round(df[has_solution_mask]['ObjRatio'].std(), 4)

    runtimeMean = round(df['Runtime'].mean(), 4)
    runtimeStd = round(df['Runtime'].std(), 4)

    result_df = pd.concat([result_df, pd.DataFrame([{
        'Tuples': tuples,
        'Scenarios': scenarios,
        'Cap': cap,
        'NumberFeasibleMean': numFeasibleMean,
        'NumberFeasibleStd': numFeasibleStd,
        'ObjRatioMean': objRatioMean,
        'ObjRatioStd': objRatioStd,
        'RuntimeMean': runtimeMean,
        'RuntimeStd': runtimeStd,
    }])], ignore_index=True)

result_df.sort_values(by=['NumberFeasibleMean', 'ObjRatioMean'], ascending=False, inplace=True)
result_df.to_csv(output_csv_path, index=False)

display_df = result_df[['Tuples', 'Scenarios', 'Cap']].copy()
display_df['NumberFeasible'] = result_df['NumberFeasibleMean'].astype(str) + ' ± ' + result_df['NumberFeasibleStd'].astype(str)
display_df['ObjRatio'] = result_df['ObjRatioMean'].astype(str) + ' ± ' + result_df['ObjRatioStd'].astype(str)
display_df['Runtime'] = result_df['RuntimeMean'].astype(str) + ' ± ' + result_df['RuntimeStd'].astype(str)


print(display_df)

fig, ax = plt.subplots(figsize=(20, len(display_df) * 0.5 + 1))
ax.axis('off')
tbl = ax.table(cellText=display_df.values, colLabels=display_df.columns, loc='center', cellLoc='center')
tbl.auto_set_font_size(False)
tbl.set_fontsize(7)
tbl.auto_set_column_width(col=list(range(len(display_df.columns))))
plt.tight_layout()
plt.savefig(output_pdf_path, bbox_inches='tight')
