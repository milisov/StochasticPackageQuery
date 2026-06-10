import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('/home/fm2288/StochasticPackageQuery/test/HyperParameterBenchmark2/benchmark_secondary_summary_tuples.csv')
df.sort_values(by='RuntimeMean', inplace=True)

fig, ax = plt.subplots(figsize=(12, 6))

ax.plot(df['RuntimeMean'], df['ObjRatioMean'], marker='o', linestyle='-')

for _, row in df.iterrows():
    ax.annotate(str(row['NumberFeasibleMean']),
                xy=(row['RuntimeMean'], row['ObjRatioMean']),
                xytext=(0, 8), textcoords='offset points',
                ha='center', fontsize=8)
    ax.annotate(str(row['Tuples']),
                xy=(row['RuntimeMean'], row['ObjRatioMean']),
                xytext=(0, -12), textcoords='offset points',
                ha='center', fontsize=8)

ax.set_xscale('log')
ax.set_xlabel('Runtime Mean (log scale)')
ax.set_ylabel('ObjRatio Mean')
ax.set_title('ObjRatio vs Runtime (annotated with Mean Feasibility)')
plt.tight_layout()
plt.savefig('/home/fm2288/StochasticPackageQuery/test/HyperParameterBenchmark2/benchmark_secondary_plot.pdf', bbox_inches='tight')
plt.show()
