from multiprocessing import Pool
import warnings
import time
import os
import numpy as np
import re
import argparse
import json
import csv
import copy
import matplotlib.pyplot as plt
import pandas as pd

from CVaRification.RCLSolve import RCLSolve
from DbInfo.PortfolioInfo import PortfolioInfo
from DbInfo.TpchInfo import TpchInfo
from Naive.Naive import Naive
from SummarySearch.SummarySearch import SummarySearch
from OfflinePreprocessing.DistPartition import DistPartition
from PgConnection.PgConnection import PgConnection
from DbInfo.DbInfo import DbInfo
# from QueryHardness.HardnessEvaluator import HardnessEvaluator
# from QueryHardness.RCLSolveBasedHardness import RCLSolveBasedHardness 
from ScenarioGenerator.PorfolioScenarioGenerator.GainScenarioGenerator import GainScenarioGenerator
from ScenarioGenerator.TpchScenarioGenerators.PriceScenarioGenerator import PriceScenarioGenerator
from SeedManager.SeedManager import SeedManager
from OfflinePreprocessing.MonotonicDequeUnitTest import MonotonicDequeUnitTest
from OfflinePreprocessing.OptimalPartitioningUnitTest import OptimalPartitioningUnitTest
from StochasticPackageQuery.Parser.Parser import Parser
from Utils.RelationalOperators import RelationalOperators
from Utils.Stochasticity import Stochasticity
from ValueGenerator.ValueGenerator import ValueGenerator
from Validator.Validator import Validator
from SketchRefine.Sketch import Sketch
from SketchRefine.SketchRefine import SketchRefine

# Directory Constants
WORKLOAD_BASE_STD = '/home/fm2288/StochasticPackageQuery/test/Queries'
WORKLOAD_BASE_SEED = '/home/fm2288/StochasticPackageQuery/test/QueriesSeeds'
TIMEOUT_SECONDS = 10 * 60

# --- Helper Functions ---

def get_and_plot_profit_stddev(table_name="stocks_5_validate", column_name="profit"):
    """
    Queries the database to get the standard deviation per id,
    then plots a histogram of stddev values across all ids.

    Args:
        table_name: The table to query (default: stocks_5_validate)
        column_name: The column to compute stddev for (default: profit)

    Returns:
        List of (id, stddev) tuples.
    """
    # Get standard deviation per id
    sql = f"""
        SELECT id, STDDEV(val) as std
        FROM {table_name}, unnest({column_name}) AS val
        GROUP BY id
        ORDER BY id;
    """
    print("I am running the query")
    PgConnection.Execute(sql)
    results = PgConnection.Fetch()

    ids = [row[0] for row in results]
    stddevs = [float(row[1]) if row[1] else 0.0 for row in results]

    print("I am done running the query")

    # Box plot with log scale
    plt.figure(figsize=(10, 6))
    plt.boxplot(stddevs, vert=True)
    plt.yscale('log')
    plt.ylabel('Standard Deviation (log scale)')
    plt.title(f'{column_name} StdDev Distribution in {table_name}')
    plt.tight_layout()
    plt.savefig(f'{table_name}_{column_name}_stddev_boxplot.png')
    plt.show()

    print(f"Number of IDs: {len(ids)}")
    print(f"Mean StdDev: {np.mean(stddevs):.4f}")
    print(f"Min: {min(stddevs):.4f}, Max: {max(stddevs):.4f}")

    # Percentiles
    percentiles = [90, 95, 97, 99, 99.5, 99.7, 99.8, 99.9]
    print("\nPercentiles:")
    for p in percentiles:
        val = np.percentile(stddevs, p)
        print(f"  {p}th: {val:.4f}")

    return list(zip(ids, stddevs))

def format_solution_list(solutions):
    """Formats solution tuples into a semi-colon separated string for CSV."""
    if not solutions:
        return "[]"
    formatted = []
    for sol in solutions:
        rt = f"{sol[0]:.5f}"
        feas = str(sol[1]).replace(',', ';')
        surp = str(sol[2]).replace(',', ';')
        obj = f"{sol[3]:.5f}"
        oos_deter_feas = int(sol[4])
        oos_prob_feas = int(sol[5]) 
        oos_feas = str(sol[6]).replace(',', ';')
        oos_surp = str(sol[7]).replace(',', ';')
        oos_obj = f"{sol[8]:.5f}"
        no_of_scenarios = sol[9]
        formatted.append(f"({rt};{feas};{surp};{obj};{oos_deter_feas};{oos_prob_feas};{oos_feas};{oos_surp};{oos_obj};{no_of_scenarios})")
    return "[" + ",".join(formatted) + "]"

def validate_and_prepare_row(package_dict, query, hardness, runtime, solutions_history, seed=None, algorithm=None):
    """
    Validates the result and returns the row dictionary.
    Does NOT write to file directly (to avoid race conditions).
    """
    
    # Strip _validate suffix if present
    if "validate" in query.get_relation():
        query.set_relation(query.get_relation().replace("_validate", ""))

    valid_query = copy.deepcopy(query)

    baseName = query.get_relation().split("_seeded_")[0]  # e.g., stocks_3_100
    parts = baseName.split("_")  # ['stocks', '3', '100']
    N = int(parts[1])
    validateTableName = f"{parts[0]}_{parts[1]}_validate"  # stocks_3_validate
    valid_query.set_relation(validateTableName)

    # Create a fresh Validator instance (thread/process safe)
    validator = Validator(
        query=valid_query,
        dbInfo=PortfolioInfo,
        no_of_validation_scenarios=10000
    )

    deter_feasible, prob_feasible, feasibilities, surpluses, objective_value = \
        validator.validate_package(package_dict)

    # Compute ObjRatio against DETER baseline
    obj_ratio = None
    if algorithm is None: 
        deter_csv_path = os.path.join("results", "DETER", "DETER.csv")
        if os.path.exists(deter_csv_path):
            deter_df = pd.read_csv(deter_csv_path)
            match = deter_df[(deter_df['N'] == N) & (deter_df['Hardness'] == hardness)]
            if not match.empty:
                deter_obj = match['Objective'].mean()
                if deter_obj != 0:
                    obj_ratio = round(objective_value / deter_obj, 5)

    print(f"Validation Results -> Hardness: {hardness}, Seed: {seed}, Obj: {objective_value}, ObjRatio: {obj_ratio}, DeterFeas: {deter_feasible}, ProbFeas: {prob_feasible}")

    row = None
    if algorithm and algorithm.upper() == 'DETER':
        row = {
            'Hardness': hardness,
            'Objective': objective_value,
            'deter_feas': int(deter_feasible),
            'prob_feas': int(prob_feasible),
            'feas': str(np.round(feasibilities, 4)).replace(',', ';'),
            'surplus': str(np.round(surpluses, 4)).replace(',', ';'),
        }
    else:
        # Prepare Row
        row = {
            'Hardness': hardness,
            'Objective': objective_value,
            'ObjRatio': obj_ratio,
            'deter_feas': int(deter_feasible),
            'prob_feas': int(prob_feasible),
            'feas': str(np.round(feasibilities, 4)).replace(',', ';'),
            'surplus': str(np.round(surpluses, 4)).replace(',', ';'),
            'Runtime': runtime,
            'solutions': format_solution_list(solutions_history)
        }

    # Add Seed if it exists
    if seed is not None and algorithm is None:
        row['Seed'] = seed

    return row

def process_single_task(args):
    """
    This is the Worker Function that runs on a separate core.
    It performs the Solve + Validation and returns the result Row.
    """
    full_path, query_name, h, seed, algorithm, M, N, results_base = args
    
    try:
        print(f"Starting Processing: {os.path.basename(full_path)} (H={h}, Seed={seed})")
        
        with open(full_path, 'r') as f:
            query = Parser().parse(f.readlines())
            
            SeedManager.reinitialize_seed() 
            start = time.time()

            if algorithm.upper() == 'HARDNESS':
                return None
            else:
                # Solver Logic
                solver = None
                # Base parameters
                params = {
                    'query': query, 'linear_relaxation': False, 'dbInfo': PortfolioInfo,
                    'init_no_of_scenarios': min(M, 100), 'no_of_validation_scenarios': M,
                    'approximation_bound': 0.0
                }

                # Algorithm Mapping
                is_rcl_type = algorithm.upper() in ['RCL', 'RCLSEED', 'RCLSEEDTIMEBUDGET']
                is_ss_type = algorithm.upper() in ['SS', 'SSSEED', 'DETER']
                is_sketchrefine_type = algorithm.upper() in ['SKETCHREFINESEED']

                if is_rcl_type:
                    if algorithm.upper() == 'RCLSEEDTIMEBUDGET':
                        params.update({'sampling_tolerance': 0.2, 'bisection_threshold': 0.1, 'start_time': start, 'timeout': TIMEOUT_SECONDS, 'timeBudgeted': True})
                    else:
                        params.update({'sampling_tolerance': 0.2, 'bisection_threshold': 0.1})
                    solver = RCLSolve(**params)
                elif algorithm.upper() == 'NAIVE':
                    solver = Naive(**params)
                elif algorithm.upper() == 'DETER':
                    baseName = query.get_relation().split("_seeded_")[0]
                    parts = baseName.split("_")
                    query.set_relation(f"{parts[0]}_{parts[1]}_validate")
                    params.update({'init_no_of_scenarios': 10**4, 'no_of_validation_scenarios': 10**4, 'init_no_of_summaries': 1})
                    solver = SummarySearch(**params)
                elif is_ss_type:
                    params['init_no_of_summaries'] = 1
                    solver = SummarySearch(**params)
                elif is_sketchrefine_type:
                    solver = SketchRefine(query=query, dbInfo=PortfolioInfo)
                else:
                    return None

                if is_rcl_type:
                    solver.solve(can_add_scenarios=True)
                elif algorithm.upper() == 'DETER':
                    solver.solve(start_time=start, timeout=TIMEOUT_SECONDS, getDeterministic=True)
                elif is_sketchrefine_type:
                    solver.solve()
                else:
                    solver.solve(start_time=start, timeout=TIMEOUT_SECONDS)

                # Log Metrics
                if hasattr(solver, 'get_metrics'):
                    metrics = solver.get_metrics()
                    metrics.set_query(query_name)
                    metrics.set_hardness(h)
                    
                    rt_ms = metrics.get_runtime() * 1000
                    
                    # # Save JSON log (Filenames are unique, so this is safe in parallel)
                    # json_name = f"{query_name}_{algorithm}_{h}.json"
                    # metrics.log_to_json(os.path.join(results_base, json_name), rt_ms)
                    
                    algo = algorithm.upper() if algorithm.upper() == "DETER" else None
                    # Return the row data to the main process
                    row = validate_and_prepare_row(
                        metrics.get_package(), query, h, rt_ms,
                        metrics.get_solutions(), seed, algo
                    )
                    if algorithm.upper() == 'DETER' and row is not None:
                        row['N'] = N
                    return row
    except Exception as e:
        print(f"Error processing {full_path}: {e}")
        return None


def run_experiment(workload_directory, algorithm, M, N, relation_name):
    print(f"--- Running Parallel Experiment 90 Cores) ---")
    print(f"Algorithm: {algorithm}")
    print(f"Workload Base: {workload_directory}")

    is_seeded = algorithm.upper() in ['RCLSEED', 'SSSEED', 'RCLSEEDTIMEBUDGET', 'SKETCHREFINESEED', 'DETER']

    # Define paths for Results
    if algorithm.upper() == "DETER":
        results_dir = os.path.join("results", "DETER")
        csv_path = os.path.join(results_dir, f"{algorithm.upper()}.csv")
    else:
        results_dir = os.path.join("results", "ImputedDataTolerance")
        csv_path = os.path.join(results_dir, f"{algorithm.upper()}_{relation_name}.csv")
    os.makedirs(results_dir, exist_ok=True)

    # --- 1. Prepare List of Tasks ---
    tasks = []
    if is_seeded:
        print("Detected Seeded Experiment. Generating task list...")
        for h in range(0, 6): # Hardness 0 to 5
            for seed in range(1, 11): # Seed 1 to 10
                filename = f"stocks_{N}_{M}_seeded_{seed}.spaql"
                full_path = os.path.join(workload_directory, str(h), filename)
                query_name = f"stocks_{N}_{M}_seeded_{seed}"
                
                if os.path.exists(full_path):
                    tasks.append((full_path, query_name, h, seed, algorithm, M, N, results_dir))
    else:
        print("Detected Standard Algorithm. Generating task list...")
        if os.path.exists(workload_directory):
            for f in sorted(os.listdir(workload_directory)):
                match = re.search(r"(stocks_\d_\d)_(.*\d+).spaql", f)
                if match:
                    query_name = match.group(1)
                    h = int(match.group(2))
                    if h >= 0:
                        full_path = os.path.join(workload_directory, f)
                        tasks.append((full_path, query_name, h, None, algorithm, M, N, results_dir))
    
    print(f"Total Tasks Found: {len(tasks)}")

    # --- 2. Initialize CSV Writer ---
    with open(csv_path, 'a', newline='') as csvfile:
        if algorithm.upper() == 'HARDNESS':
            fieldnames = ['hardness', 'computed_hardness']
        elif algorithm.upper() == 'DETER':
            fieldnames = ['N', 'Hardness', 'Objective', 'deter_feas', 'prob_feas', 'feas', 'surplus']
        else:
            fieldnames = ['Hardness', 'Objective', 'ObjRatio', 'deter_feas', 'prob_feas', 'feas', 'surplus', 'Runtime', 'solutions']
            if is_seeded:
                fieldnames.insert(1, 'Seed')
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        csvfile.flush() # <--- FORCE WRITE TO DISK IMMEDIATELY
        print(f"Header written to {csv_path}")

        # --- 3. Parallel Execution ---
        NUM_PROCESSES = 15
        if N == 3:
            NUM_PROCESSES = 30
        elif N == 5 and M == 10000:
            NUM_PROCESSES = 6
        print(f"Launching Pool with {NUM_PROCESSES} processes...")

        with Pool(processes=NUM_PROCESSES) as pool:
            # imap_unordered yields results as soon as they finish
            for result in pool.imap_unordered(process_single_task, tasks):
                if result:
                    writer.writerow(result)
                    csvfile.flush() # Ensure data is written immediately
                    
    print("Experiment Complete.")

def run_all_rs_seed():
    Ns = [3, 4, 5]
    Ms = [10, 100, 10000]
    for N in Ns:
        for M in Ms:
            rel_name = f"stocks_{N}_{M}_seeded"
            workload_dir = os.path.join(WORKLOAD_BASE_SEED, rel_name, "RCL")
            if os.path.exists(workload_dir):
                run_experiment(workload_dir, 'RCLSeedTimeBudget', M, N, rel_name)
            else:
                print(f"Warning: Workload directory not found: {workload_dir}")

if __name__ == '__main__':
    warnings.filterwarnings('ignore')
    parser = argparse.ArgumentParser()
    parser.add_argument('algorithm', type=str, choices=['RCL', 'SS', 'Naive', 'DETER', 'HARDNESS', 'RCLSeed', 'SSSeed', 'RCLSeedTimeBudget', 'SketchRefineSeed', 'stddev', 'AllRSSeed'])
    parser.add_argument('N', type=int, nargs='?', default=5)
    parser.add_argument('M', type=int, nargs='?', default=100)
    parser.add_argument('--seeded', action='store_true', help='Run partitioning for all 10 seeds (only for partition command)')
    args = parser.parse_args()

    # --- Handle stddev command ---
    if args.algorithm.lower() == 'stddev':
        get_and_plot_profit_stddev(f"stocks_{args.N}_validate", "profit")
    elif args.algorithm == 'AllRSSeed':
        run_all_rs_seed()
    else:
        # --- Directory Logic ---
        if args.algorithm.upper() in ['RCLSEED', 'SSSEED', 'RCLSEEDTIMEBUDGET', 'SKETCHREFINESEED', "DETER"]:
            rel_name = f"stocks_{args.N}_{args.M}_seeded"
            workload_dir = os.path.join(WORKLOAD_BASE_SEED, rel_name, "RCL")
        else:
            rel_name = f"stocks_{args.N}_{args.M}"
            workload_dir = os.path.join(WORKLOAD_BASE_STD, rel_name, "RCL")
        if not os.path.exists(workload_dir):
            print(f"Error: Workload directory not found: {workload_dir}")
        else:
            run_experiment(workload_dir, args.algorithm, args.M, args.N, rel_name)