import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

ablation_base = "/home/fm2288/StochasticPackageQuery/test/AblationStudyTerminationCondition/"
rs_mad_base = "/home/fm2288/StochasticPackageQuery/test/RS_MAD_upto_1024/"

ablation_dirs = ["stage1fixedN", "stage1lcvar", "stage1random"]
tuple_prefixes = ["3", "4", "5"]
thresholds = [("0.010000", "0.01"), ("0.100000", "0.1")]

annotation_offsets = [8, 18, 28, 38]

script_dir = os.path.dirname(os.path.abspath(__file__))

for thresh_dir, thresh_label in thresholds:
    output_pdf_path = os.path.join(script_dir, f"tuple_ablation_threshold_{thresh_label}.pdf")

    fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharey=False)

    for ax, tprefix in zip(axes, tuple_prefixes):
        ablation_files = [f"stocks_{tprefix}_10.csv", f"stocks_{tprefix}_100.csv", f"stocks_{tprefix}_10000.csv"]
        rs_files = [f"RSstocks_{tprefix}_10.csv", f"RSstocks_{tprefix}_100.csv", f"RSstocks_{tprefix}_10000.csv"]

        methods = []
        for comp_dir in ablation_dirs:
            methods.append((comp_dir, os.path.join(ablation_base, comp_dir, f"threshold_{thresh_dir}/"), ablation_files))
        methods.append(("RS_MAD_upto_1024", os.path.join(rs_mad_base, f"threshold_{thresh_dir}/"), rs_files))

        for idx, (method_name, method_dir, fnames) in enumerate(methods):
            all_dfs = []
            for fname in fnames:
                fpath = os.path.join(method_dir, fname)
                if os.path.exists(fpath):
                    df = pd.read_csv(fpath)
                    df['File'] = fname.split(".")[0]
                    all_dfs.append(df)

            if not all_dfs:
                continue

            combined = pd.concat(all_dfs, ignore_index=True)
            feasible_mask = (combined['deter_feas'] == 1) & (combined['prob_feas'] == 1)
            feasible_df = combined[feasible_mask]

            total_feasible = feasible_df.groupby('Hardness').size()
            mean_obj_ratio = feasible_df.groupby('Hardness')['ObjRatio'].mean().round(4)
            q1_obj_ratio = feasible_df.groupby('Hardness')['ObjRatio'].quantile(0.25)
            q3_obj_ratio = feasible_df.groupby('Hardness')['ObjRatio'].quantile(0.75)

            if mean_obj_ratio.empty:
                continue

            line, = ax.plot(mean_obj_ratio.index, mean_obj_ratio.values, marker='o', label=method_name)
            color = line.get_color()

            ax.fill_between(
                mean_obj_ratio.index,
                q1_obj_ratio.reindex(mean_obj_ratio.index),
                q3_obj_ratio.reindex(mean_obj_ratio.index),
                alpha=0.15,
                color=color,
            )

            for h in mean_obj_ratio.index:
                ax.annotate(
                    f'{total_feasible.get(h, 0)}',
                    xy=(h, mean_obj_ratio[h]),
                    xytext=(0, annotation_offsets[idx]),
                    textcoords='offset points',
                    ha='center',
                    fontsize=8,
                    color=color,
                )

        ax.set_xlabel('Hardness')
        ax.set_ylabel('Mean ObjRatio (Feasible Solutions)')
        ax.set_title(f'N = {tprefix}')
        ax.set_xticks(range(6))
        ax.legend()
        ax.grid(True, alpha=0.3)

    fig.suptitle(f'Mean Objective Ratio by Hardness (threshold={thresh_label})', fontsize=13)
    plt.tight_layout()
    plt.savefig(output_pdf_path, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_pdf_path}")
