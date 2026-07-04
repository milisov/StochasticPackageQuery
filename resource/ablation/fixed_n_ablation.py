import pandas as pd
import os
import matplotlib.pyplot as plt

ablation_base = "/home/fm2288/StochasticPackageQuery/test/AblationStudyFixedN"
rs_mad_base   = "/home/fm2288/StochasticPackageQuery/test/RS_MAD"

N_values = sorted([int(d) for d in os.listdir(ablation_base) if d.isdigit()])
tables   = ["5_10", "5_100", "5_10000"]
methods  = [
    ("stage1lcvar",  "lcvar"),
    ("stage1random", "random"),
    ("RS_MAD",       "RS_MAD"),
]

annotation_y_offsets = [10, 22, -14]

script_dir = os.path.dirname(os.path.abspath(__file__))

fig, axes = plt.subplots(1, 3, figsize=(18, 6))

for ax, table in zip(axes, tables):
    for (method_key, method_label), y_off in zip(methods, annotation_y_offsets):
        xs, ys, ns = [], [], []
        for N in N_values:
            if method_key == "RS_MAD":
                fpath = os.path.join(rs_mad_base, str(N), f"RSstocks_{table}.csv")
            else:
                fpath = os.path.join(ablation_base, str(N), method_key, f"stocks_{table}.csv")
            if not os.path.exists(fpath):
                continue
            df = pd.read_csv(fpath)
            feas = df[(df["deter_feas"] == 1) & (df["prob_feas"] == 1)]
            if feas.empty:
                continue
            xs.append(df["Runtime"].mean())
            ys.append(feas["ObjRatio"].mean())
            ns.append(len(feas))

        if not xs:
            continue
        line, = ax.plot(xs, ys, marker="o", label=method_label)
        color = line.get_color()
        for x, y, n in zip(xs, ys, ns):
            ax.annotate(str(n), xy=(x, y),
                        xytext=(0, y_off), textcoords="offset points",
                        ha="center", fontsize=7, color=color)

    ax.set_xscale("log")
    ax.set_xlabel("Avg Runtime (ms, log scale)")
    ax.set_ylabel("Avg ObjRatio (feasible)")
    ax.set_title(f"stocks_{table}")
    ax.legend()

fig.tight_layout()
out = os.path.join(script_dir, "fixed_n_ablation.pdf")
fig.savefig(out)
print(f"Saved {out}")
