import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

base_dir = "/home/fm2288/StochasticPackageQuery/test/AblationStudy2/"
comparisonDirs = ["stage1fixedN", "stage1lcvar", "stage1random"]
tuple_prefixes = ["3", "4", "5"]

output_csv_path = os.path.join(base_dir, "tuple_ablation_results.csv")
output_pdf_path = os.path.join(base_dir, "tuple_ablation_plot.pdf")

# Vertical offsets per method so annotations don't stack on top of each other
annotation_offsets = [8, 18, 28]

all_results = []
fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharey=False)

for ax, tprefix in zip(axes, tuple_prefixes):
    files_to_concat = [f"stocks_{tprefix}_10.csv", f"stocks_{tprefix}_100.csv", f"stocks_{tprefix}_10000.csv"]

    for idx, comp_dir in enumerate(comparisonDirs):
        all_dfs = []
        for fname in files_to_concat:
            df = pd.read_csv(os.path.join(base_dir, comp_dir, fname))
            df['File'] = fname.split(".")[0]
            all_dfs.append(df)

        combined = pd.concat(all_dfs, ignore_index=True)
        feasible_mask = (combined['deter_feas'] == 1) & (combined['prob_feas'] == 1)
        feasible_df = combined[feasible_mask]

        # Total feasible count per hardness across all combined files
        total_feasible = feasible_df.groupby('Hardness').size()

        # Mean, Q1, Q3 of ObjRatio per hardness
        mean_obj_ratio = feasible_df.groupby('Hardness')['ObjRatio'].mean().round(4)
        q1_obj_ratio = feasible_df.groupby('Hardness')['ObjRatio'].quantile(0.25)
        q3_obj_ratio = feasible_df.groupby('Hardness')['ObjRatio'].quantile(0.75)

        for h in sorted(mean_obj_ratio.index):
            all_results.append({
                'Tuples': tprefix,
                'Method': comp_dir,
                'Hardness': h,
                'TotalFeasible': total_feasible.get(h, 0),
                'MeanObjRatio': mean_obj_ratio.get(h, np.nan),
            })

        line, = ax.plot(mean_obj_ratio.index, mean_obj_ratio.values, marker='o', label=comp_dir)
        color = line.get_color()

        # Q1–Q3 shaded band in the same color as the line
        ax.fill_between(
            mean_obj_ratio.index,
            q1_obj_ratio.reindex(mean_obj_ratio.index),
            q3_obj_ratio.reindex(mean_obj_ratio.index),
            alpha=0.15,
            color=color,
        )

        for h in mean_obj_ratio.index:
            ax.annotate(
                f'{total_feasible[h]}',
                xy=(h, mean_obj_ratio[h]),
                xytext=(0, annotation_offsets[idx]),
                textcoords='offset points',
                ha='center',
                fontsize=8,
                color=color,
            )

    ax.set_xlabel('Hardness')
    ax.set_ylabel('Mean ObjRatio (Feasible Solutions)')
    ax.set_title(f'Tuples = {tprefix}')
    ax.set_xticks(range(6))
    ax.legend()
    ax.grid(True, alpha=0.3)

fig.suptitle('Mean Objective Ratio by Hardness', fontsize=13)
plt.tight_layout()
plt.savefig(output_pdf_path, bbox_inches='tight')
plt.show()

result_df = pd.DataFrame(all_results)
result_df.to_csv(output_csv_path, index=False)
print(result_df)
