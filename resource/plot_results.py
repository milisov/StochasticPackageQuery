"""
Plot results for RCL (seeded) vs MAD (RS) algorithms.

Directory structure expected:
  mad_dir/   -> RSstocks_{i}_{s}.csv           (i in 3,4,5 ; s in 10,100,10000)
  rcl_dir/   -> RCLSEEDTIMEBUDGET_stocks_{i}_{s}_seeded.csv

Outputs (saved to output_dir):
  objective_plots.pdf   - 9 subplots (one per instance), objective vs hardness
  surplus_plots.pdf     - 9 subplots, surplus vs hardness  (99% CI)
"""

import os
import ast
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

MAD_DIR    = "/home/fm2288/StochasticPackageQuery/test/ExpImputedDataFinalAlgo/RobustSatisficing"
RCL_DIR    = "/home/fm2288/StochasticPackageQuery/stochastic-sketchrefine/results/TimeBudgetTenImputedData"
OUTPUT_DIR = "/home/fm2288/StochasticPackageQuery/test/Plots"

TUPLES    = [3, 4, 5]
SCENARIOS = [10, 100, 10000]
HARDNESS  = [0, 1, 2, 3, 4, 5]

FEAS_THRESH = 0.95
CI_Z        = 2.576          # z for 99 % CI

COLORS = {"RCL": "#3a7bd5", "MAD": "#e05252"}
LABELS = {"RCL": "RCL",     "MAD": "MAD"}

os.makedirs(OUTPUT_DIR, exist_ok=True)


# ══════════════════════════════════════════════════════════════
# 1.  DATA LOADING
# ══════════════════════════════════════════════════════════════

def parse_list_col(val):
    try:
        return ast.literal_eval(str(val).strip())
    except Exception:
        return [np.nan]


def _tidy(df):
    df.columns     = df.columns.str.strip()
    df["feas_val"] = df["feas"].apply(parse_list_col).apply(
                         lambda x: x[0] if x else np.nan)
    df["surplus_val"] = df["surplus"].apply(parse_list_col).apply(
                            lambda x: x[0] if x else np.nan)
    df["Objective"] = pd.to_numeric(df["Objective"], errors="coerce")
    df["Hardness"]  = pd.to_numeric(df["Hardness"],  errors="coerce").astype(int)
    return df


def load_all():
    data = {}
    for i in TUPLES:
        for s in SCENARIOS:
            entry = {}

            rcl_path = os.path.join(RCL_DIR, f"RCLSEEDTIMEBUDGET_stocks_{i}_{s}_seeded.csv")
            if os.path.exists(rcl_path):
                entry["RCL"] = _tidy(pd.read_csv(rcl_path))
            else:
                print(f"[WARN] RCL file not found: {rcl_path}")

            mad_path = os.path.join(MAD_DIR, f"RSstocks_{i}_{s}.csv")
            if os.path.exists(mad_path):
                entry["MAD"] = _tidy(pd.read_csv(mad_path))
            else:
                print(f"[WARN] MAD file not found: {mad_path}")

            data[(i, s)] = entry
    return data


def aggregate_objective(df):
    rows = []
    for h in HARDNESS:
        grp      = df[df["Hardness"] == h]
        feasible = grp[grp["feas_val"] >= FEAS_THRESH]["Objective"].dropna()
        n_feas   = len(feasible)
        mean_obj = feasible.mean() if n_feas > 0 else np.nan
        se       = feasible.std(ddof=1) / np.sqrt(n_feas) if n_feas > 1 else 0.0
        ci       = CI_Z * se
        rows.append({"Hardness": h, "n_feasible": n_feas, "mean_obj": mean_obj, "ci_half": ci})
    return pd.DataFrame(rows).set_index("Hardness")


def aggregate_surplus(df):
    rows = []
    for h in HARDNESS:
        vals = df[df["Hardness"] == h]["surplus_val"].dropna()
        n    = len(vals)
        m    = vals.mean() if n > 0 else np.nan
        se   = vals.std(ddof=1) / np.sqrt(n) if n > 1 else 0.0
        rows.append({"Hardness": h, "mean": m, "ci_half": CI_Z * se})
    return pd.DataFrame(rows).set_index("Hardness")


# ══════════════════════════════════════════════════════════════
# 3.  SINGLE-AXIS PLOT HELPERS
# ══════════════════════════════════════════════════════════════

_OBJ_OFFSET = {"RCL": 8, "MAD": -14}


def _style_ax(ax):
    ax.set_xticks(HARDNESS)
    ax.xaxis.set_major_locator(ticker.FixedLocator(HARDNESS))
    ax.set_xlabel("Hardness", fontsize=9)
    ax.legend(fontsize=8)
    ax.grid(axis="y", linestyle="--", alpha=0.4)


def plot_objective_ax(ax, entry, title):
    for alg in ("RCL", "MAD"):
        df = entry.get(alg)
        if df is None:
            continue
        agg = aggregate_objective(df)
        x  = np.array(HARDNESS)
        y  = agg["mean_obj"].reindex(HARDNESS).values
        n  = agg["n_feasible"].reindex(HARDNESS).fillna(0).astype(int).values
        ci = agg["ci_half"].reindex(HARDNESS).fillna(0).values

        ax.plot(x, y, "o-", color=COLORS[alg], label=LABELS[alg],
                linewidth=2, markersize=6, zorder=3)
        ax.fill_between(x, y - ci, y + ci,
                        color=COLORS[alg], alpha=0.15, zorder=2)

        for xi, yi, ni in zip(x, y, n):
            if not np.isnan(yi):
                ax.annotate(str(ni), xy=(xi, yi),
                            xytext=(4, _OBJ_OFFSET[alg]),
                            textcoords="offset points",
                            color=COLORS[alg], fontsize=8, fontweight="bold")

    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.set_ylabel("Objective (feasible seeds only)", fontsize=9)
    _style_ax(ax)


def plot_surplus_ax(ax, entry, title):
    for alg in ("RCL", "MAD"):
        df = entry.get(alg)
        if df is None:
            continue
        agg = aggregate_surplus(df)
        x   = np.array(HARDNESS)
        y   = agg["mean"].reindex(HARDNESS).values
        ci  = agg["ci_half"].reindex(HARDNESS).fillna(0).values

        ax.plot(x, y, "o-", color=COLORS[alg], label=LABELS[alg],
                linewidth=2, markersize=5, zorder=3)
        ax.fill_between(x, y - ci, y + ci,
                        color=COLORS[alg], alpha=0.15, zorder=2)

    ax.axhline(0, color="red", linestyle="--", linewidth=1.2, zorder=4)
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.set_ylabel("Surplus", fontsize=9)
    _style_ax(ax)


# ══════════════════════════════════════════════════════════════
# 4.  MAIN
# ══════════════════════════════════════════════════════════════

def main():
    print("Loading data …")
    data = load_all()

    instance_list = [(i, s) for i in TUPLES for s in SCENARIOS]   # 9 instances

    # ── Figure 1 : Objective plots ────────────────────────────
    print("Building objective plots …")
    fig1, axes1 = plt.subplots(3, 3, figsize=(15, 13))
    fig1.suptitle(
        "Objective by Hardness — feasible solutions only, 10-min budget (99% CI)",
        fontsize=13, fontweight="bold"
    )
    for idx, (i, s) in enumerate(instance_list):
        plot_objective_ax(axes1[idx // 3][idx % 3], data.get((i, s), {}), f"{i}_{s}")

    fig1.tight_layout(rect=[0, 0, 1, 0.95])
    obj_path = os.path.join(OUTPUT_DIR, "objective_plots.pdf")
    fig1.savefig(obj_path, bbox_inches="tight")
    print(f"  Saved → {obj_path}")

    # ── Figure 2 : Surplus plots ──────────────────────────────
    print("Building surplus plots …")
    fig2, axes2 = plt.subplots(3, 3, figsize=(15, 13))
    fig2.suptitle(
        "Surplus by Hardness — 10-min budget (99% CI)",
        fontsize=13, fontweight="bold"
    )
    for idx, (i, s) in enumerate(instance_list):
        plot_surplus_ax(axes2[idx // 3][idx % 3], data.get((i, s), {}), f"{i}_{s}")

    fig2.tight_layout(rect=[0, 0, 1, 0.95])
    sur_path = os.path.join(OUTPUT_DIR, "surplus_plots.pdf")
    fig2.savefig(sur_path, bbox_inches="tight")
    print(f"  Saved → {sur_path}")

    plt.close("all")
    print("\nDone! All plots saved to:", OUTPUT_DIR)


if __name__ == "__main__":
    main()