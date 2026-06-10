"""
Plot ablation results: MAD (ImputedGoodHyperParam) vs ablate_tuples and ablate_scenarios.

Directory structure expected:
  mad_dir/         -> RSstocks_{i}_{s}.csv
  ablation_dir/    -> stocks_{i}_{s}_scenarios.csv
                   -> stocks_{i}_{s}_tuples.csv

Outputs (saved to output_dir), one pair per time budget (2, 5, 10 min):
  ablation_objective_{t}min.pdf
  ablation_surplus_{t}min.pdf
"""

import os
import ast
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

MAD_DIR      = "/home/fm2288/StochasticPackageQuery/test/ImputedGoodHyperParam/RobustSatisficing"
ABLATION_DIR = "/home/fm2288/StochasticPackageQuery/test/AblationStudy"
OUTPUT_DIR   = "/home/fm2288/StochasticPackageQuery/test/Plots"

TUPLES    = [3, 4, 5]
SCENARIOS = [10, 100, 10000]
HARDNESS  = [0, 1, 2, 3, 4, 5]

FEAS_THRESH = 0.95
CI_Z        = 2.576          # z for 99 % CI

TIME_BUDGETS = [
    ("30s", 30  * 1000),
    ("1min",  1 * 60 * 1000),
    ("2min",  2 * 60 * 1000),
    ("5min",  5 * 60 * 1000),
    ("10min", 10 * 60 * 1000),
]

COLORS = {
    "MAD":               "#e05252",
    "Ablation-Scenarios": "#2ca02c",
    "Ablation-Tuples":    "#ff7f0e",
}
LABELS = {
    "MAD":               "MAD",
    "Ablation-Scenarios": "No Scenarios",
    "Ablation-Tuples":    "No Tuples",
}
_OBJ_OFFSETS = {
    "MAD":               -14,
    "Ablation-Scenarios":  8,
    "Ablation-Tuples":   -26,
}

ALG_ORDER = ("MAD", "Ablation-Scenarios", "Ablation-Tuples")

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
    df.columns = df.columns.str.strip()
    df["Objective"] = pd.to_numeric(df["Objective"], errors="coerce")
    df["Hardness"]  = pd.to_numeric(df["Hardness"],  errors="coerce").astype(int)
    return df


def load_all():
    data = {}
    for i in TUPLES:
        for s in SCENARIOS:
            entry = {}

            mad_path = os.path.join(MAD_DIR, f"RSstocks_{i}_{s}.csv")
            if os.path.exists(mad_path):
                entry["MAD"] = _tidy(pd.read_csv(mad_path))
            else:
                print(f"[WARN] MAD file not found: {mad_path}")

            for suffix, key in (
                ("scenarios", "Ablation-Scenarios"),
                ("tuples",    "Ablation-Tuples"),
            ):
                path = os.path.join(ABLATION_DIR, f"stocks_{i}_{s}_{suffix}.csv")
                if os.path.exists(path):
                    entry[key] = _tidy(pd.read_csv(path))
                else:
                    print(f"[WARN] Ablation file not found: {path}")

            data[(i, s)] = entry
    return data


# ══════════════════════════════════════════════════════════════
# 2.  TIME-BUDGET SOLUTION SELECTION
# ══════════════════════════════════════════════════════════════

def _parse_solution_entry(entry):
    """Parse one entry like '1197.86;[1];[0.05];16003.3;1;1;[1];[0.05];15893.4;500;50'."""
    parts = entry.split(';')
    # fields: 0=ts, 1=optim_feas, 2=optim_surplus, 3=optim_obj,
    #         4=deter_feas, 5=prob_feas, 6=valid_feas, 7=valid_surplus,
    #         8=valid_obj, 9=n_tuples, 10=n_scenarios
    return {
        'timestamp_ms':  float(parts[0]),
        'valid_surplus': float(parts[7].strip().strip('[]')),
        'valid_obj':     float(parts[8]),
    }


def parse_solutions(sol_str):
    """Parse solutions column string into a list of dicts."""
    sol_str = str(sol_str).strip()
    if not sol_str or sol_str == 'nan':
        return []
    inner = sol_str[1:-1]                   # strip outer [ ]
    raw_entries = inner.split(');(')
    results = []
    for raw in raw_entries:
        entry = raw.strip().strip('(').strip(')')
        try:
            results.append(_parse_solution_entry(entry))
        except Exception:
            pass
    return results


def _is_better(a, b):
    """Return True if solution a is better than solution b."""
    a_feas = a['valid_surplus'] >= 0
    b_feas = b['valid_surplus'] >= 0
    if a_feas and not b_feas:
        return True
    if not a_feas and b_feas:
        return False
    if a_feas and b_feas:
        return a['valid_obj'] > b['valid_obj']
    return a['valid_surplus'] > b['valid_surplus']


def best_solution_at_time(sol_str, time_limit_ms):
    """Return (valid_obj, valid_surplus) of best solution found within time_limit_ms."""
    candidates = [s for s in parse_solutions(sol_str) if s['timestamp_ms'] <= time_limit_ms]
    if not candidates:
        return np.nan, np.nan
    best = candidates[0]
    for s in candidates[1:]:
        if _is_better(s, best):
            best = s
    return best['valid_obj'], best['valid_surplus']


def apply_time_budget(df, time_limit_ms):
    """Return df with best_obj and best_surplus columns from the solutions list."""
    results = df['solutions'].apply(
        lambda s: pd.Series(best_solution_at_time(s, time_limit_ms),
                            index=['best_obj', 'best_surplus'])
    )
    return df.assign(best_obj=results['best_obj'], best_surplus=results['best_surplus'])


# ══════════════════════════════════════════════════════════════
# 3.  AGGREGATION HELPERS
# ══════════════════════════════════════════════════════════════

def aggregate_objective(df, time_limit_ms):
    df = apply_time_budget(df, time_limit_ms)
    rows = []
    for h in HARDNESS:
        grp      = df[df["Hardness"] == h]
        feasible = grp[grp["best_surplus"] >= 0]["best_obj"].dropna()
        n_feas   = len(feasible)
        mean_obj = feasible.mean() if n_feas > 0 else np.nan
        se       = feasible.std(ddof=1) / np.sqrt(n_feas) if n_feas > 1 else 0.0
        rows.append({"Hardness": h, "n_feasible": n_feas, "mean_obj": mean_obj, "ci_half": CI_Z * se})
    return pd.DataFrame(rows).set_index("Hardness")


def aggregate_surplus(df, time_limit_ms):
    df = apply_time_budget(df, time_limit_ms)
    rows = []
    for h in HARDNESS:
        vals = df[df["Hardness"] == h]["best_surplus"].dropna()
        n    = len(vals)
        m    = vals.mean() if n > 0 else np.nan
        se   = vals.std(ddof=1) / np.sqrt(n) if n > 1 else 0.0
        rows.append({"Hardness": h, "mean": m, "ci_half": CI_Z * se})
    return pd.DataFrame(rows).set_index("Hardness")


# ══════════════════════════════════════════════════════════════
# 4.  SINGLE-AXIS PLOT HELPERS
# ══════════════════════════════════════════════════════════════

def _style_ax(ax):
    ax.set_xticks(HARDNESS)
    ax.xaxis.set_major_locator(ticker.FixedLocator(HARDNESS))
    ax.set_xlabel("Hardness", fontsize=9)
    ax.legend(fontsize=8)
    ax.grid(axis="y", linestyle="--", alpha=0.4)


def plot_objective_ax(ax, entry, title, time_limit_ms):
    for alg in ALG_ORDER:
        df = entry.get(alg)
        if df is None:
            continue
        agg = aggregate_objective(df, time_limit_ms)
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
                            xytext=(4, _OBJ_OFFSETS[alg]),
                            textcoords="offset points",
                            color=COLORS[alg], fontsize=8, fontweight="bold")

    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.set_ylabel("Objective (feasible seeds only)", fontsize=9)
    _style_ax(ax)


def plot_surplus_ax(ax, entry, title, time_limit_ms):
    for alg in ALG_ORDER:
        df = entry.get(alg)
        if df is None:
            continue
        agg = aggregate_surplus(df, time_limit_ms)
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
# 5.  MAIN
# ══════════════════════════════════════════════════════════════

def main():
    print("Loading data …")
    data = load_all()
    instance_list = [(i, s) for i in TUPLES for s in SCENARIOS]   # 9 instances

    for t_label, t_ms in TIME_BUDGETS:
        print(f"\n── {t_label} budget ──")

        fig1, axes1 = plt.subplots(3, 3, figsize=(15, 13))
        fig1.suptitle(
            f"Ablation — Objective by Hardness ({t_label} budget, feasible only, 99% CI)",
            fontsize=13, fontweight="bold"
        )
        for idx, (i, s) in enumerate(instance_list):
            plot_objective_ax(axes1[idx // 3][idx % 3],
                              data.get((i, s), {}), f"{i}_{s}", t_ms)
        fig1.tight_layout(rect=[0, 0, 1, 0.95])
        obj_path = os.path.join(OUTPUT_DIR, f"ablation_objective_{t_label}.pdf")
        fig1.savefig(obj_path, bbox_inches="tight")
        print(f"  Saved → {obj_path}")

        fig2, axes2 = plt.subplots(3, 3, figsize=(15, 13))
        fig2.suptitle(
            f"Ablation — Surplus by Hardness ({t_label} budget, 99% CI)",
            fontsize=13, fontweight="bold"
        )
        for idx, (i, s) in enumerate(instance_list):
            plot_surplus_ax(axes2[idx // 3][idx % 3],
                            data.get((i, s), {}), f"{i}_{s}", t_ms)
        fig2.tight_layout(rect=[0, 0, 1, 0.95])
        sur_path = os.path.join(OUTPUT_DIR, f"ablation_surplus_{t_label}.pdf")
        fig2.savefig(sur_path, bbox_inches="tight")
        print(f"  Saved → {sur_path}")

        plt.close("all")

    print("\nDone! All plots saved to:", OUTPUT_DIR)


if __name__ == "__main__":
    main()
