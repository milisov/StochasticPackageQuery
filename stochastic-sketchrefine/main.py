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

from CVaRification.CVaRification import CVaRification
from CVaRification.StaircaseCVaRification import StaircaseCVaRification
from CVaRification.RCLSolve import RCLSolve
from DbInfo.PortfolioInfo import PortfolioInfo
from DbInfo.TpchInfo import TpchInfo
from Naive.Naive import Naive
from SummarySearch.SummarySearch import SummarySearch
from OfflinePreprocessing.DistPartition import DistPartition
from PgConnection.PgConnection import PgConnection
from QueryHardness.HardnessEvaluator import HardnessEvaluator
from QueryHardness.RCLSolveBasedHardness import RCLSolveBasedHardness 
from ScenarioGenerator.PorfolioScenarioGenerator.GainScenarioGenerator import GainScenarioGenerator
from ScenarioGenerator.TpchScenarioGenerators.PriceScenarioGenerator import PriceScenarioGenerator
from SeedManager.SeedManager import SeedManager
from OfflinePreprocessing.MonotonicDequeUnitTest import MonotonicDequeUnitTest
from OfflinePreprocessing.OptimalPartitioningUnitTest import OptimalPartitioningUnitTest
from StochasticPackageQuery.Parser.Parser import Parser
from Utils.Stochasticity import Stochasticity
from UnitTestRunner import UnitTestRunner
from ValueGenerator.ValueGenerator import ValueGenerator
from Validator.Validator import Validator

# Directory Constants
WORKLOAD_BASE_STD = '/home/fm2288/StochasticPackageQuery/test/Queries'
WORKLOAD_BASE_SEED = '/home/fm2288/StochasticPackageQuery/test/QueriesSeeds'
TIMEOUT_SECONDS = 10 * 60

# --- Helper Functions ---

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
        oos_feas = str(sol[4]).replace(',', ';')
        oos_surp = str(sol[5]).replace(',', ';')
        oos_obj = f"{sol[6]:.5f}"
        no_of_scenarios = sol[7]
        formatted.append(f"({rt};{feas};{surp};{obj};{oos_feas};{oos_surp};{oos_obj};{no_of_scenarios})")
    return "[" + ",".join(formatted) + "]"

def validate_and_prepare_row(package_dict, query, hardness, runtime, solutions_history, seed=None):
    """
    Validates the result and returns the row dictionary.
    Does NOT write to file directly (to avoid race conditions).
    """
    
    # Strip _validate suffix if present
    if "validate" in query.get_relation():
        query.set_relation(query.get_relation().replace("_validate", ""))

    valid_query = copy.deepcopy(query)
    valid_query.set_relation(query.get_relation() + "_validate")
    
    # Create a fresh Validator instance (thread/process safe)
    validator = Validator(
        query=valid_query,
        dbInfo=PortfolioInfo,
        no_of_validation_scenarios=10000
    )

    feasibilities = []
    surpluses = []
    objective_value = 0.0

    if package_dict:
        for constraint in query.get_constraints():
            if constraint.is_risk_constraint():
                feas = validator.get_var_constraint_satisfaction(package_dict, constraint)
                feasibilities.append(feas)
                p = constraint.get_probability_threshold()
                surpluses.append(feas - p)
        objective_value = validator.get_validation_objective_value(package_dict)

    print(f"Validation Results -> Hardness: {hardness}, Seed: {seed}, Obj: {objective_value}")

    # Prepare Row
    row = {
        'Hardness': hardness,
        'Objective': objective_value,
        'feas': str(np.round(feasibilities, 4)).replace(',', ';'),
        'surplus': str(np.round(surpluses, 4)).replace(',', ';'),
        'Runtime': runtime,
        'solutions': format_solution_list(solutions_history)
    }
    
    # Add Seed if it exists
    if seed is not None:
        row['Seed'] = seed

    return row

def process_single_task(args):
    """
    This is the Worker Function that runs on a separate core.
    It performs the Solve + Validation and returns the result Row.
    """
    full_path, query_name, h, seed, algorithm, M, results_base = args
    
    try:
        print(f"Starting Processing: {os.path.basename(full_path)} (H={h}, Seed={seed})")
        
        with open(full_path, 'r') as f:
            query = Parser().parse(f.readlines())
            
            SeedManager.reinitialize_seed() 
            start = time.time()

            if algorithm.upper() == 'HARDNESS':
                # Hardness Evaluator Logic
                evaluator = RCLSolveBasedHardness(
                    query=query, linear_relaxation=False, dbInfo=PortfolioInfo,
                    init_no_of_scenarios=100, no_of_validation_scenarios=10**M,
                    approximation_bound=0.05, sampling_tolerance=1.00, bisection_threshold=0.1
                )
                evaluator.solve()
                comp_h = evaluator.get_model_probability()
                return {'hardness': h, 'computed_hardness': comp_h}
            
            else:
                # Solver Logic
                solver = None
                # Base parameters
                params = {
                    'query': query, 'linear_relaxation': False, 'dbInfo': PortfolioInfo,
                    'init_no_of_scenarios': 100, 'no_of_validation_scenarios': 10**M,
                    'approximation_bound': 0.0
                }

                # Algorithm Mapping
                is_rcl_type = algorithm.upper() in ['RCL', 'RCLSEED', 'RCLSEEDTIMEBUDGET']
                is_ss_type = algorithm.upper() in ['SS', 'SSSEED', 'DETER']

                if is_rcl_type:
                    if algorithm.upper() == 'RCLSEEDTIMEBUDGET':
                        params.update({'sampling_tolerance': 0.2, 'bisection_threshold': 0.1, 'start_time': start, 'timeout': TIMEOUT_SECONDS, 'timeBudgeted': True})
                    else:
                        params.update({'sampling_tolerance': 0.2, 'bisection_threshold': 0.1})
                    solver = RCLSolve(**params)
                elif algorithm.upper() == 'NAIVE':
                    solver = Naive(**params)
                elif algorithm.upper() == 'DETER':
                    query.set_relation(query.get_relation() + "_validate")
                    params.update({'init_no_of_scenarios': 10**4, 'no_of_validation_scenarios': 10**4, 'init_no_of_summaries': 1})
                    solver = SummarySearch(**params)
                elif is_ss_type:
                    params['init_no_of_summaries'] = 1
                    solver = SummarySearch(**params)
                else:
                    return None

                if is_rcl_type:
                    solver.solve(can_add_scenarios=True)
                elif algorithm.upper() == 'DETER':
                    solver.solve(start_time=start, timeout=TIMEOUT_SECONDS, getDeterministic=True)
                else:
                    solver.solve(start_time=start, timeout=TIMEOUT_SECONDS)

                # Log Metrics
                if hasattr(solver, 'get_metrics'):
                    metrics = solver.get_metrics()
                    metrics.set_query(query_name)
                    metrics.set_hardness(h)
                    
                    rt_ms = metrics.get_runtime() * 1000
                    
                    # Save JSON log (Filenames are unique, so this is safe in parallel)
                    json_name = f"{query_name}_{algorithm}_{h}.json"
                    metrics.log_to_json(os.path.join(results_base, json_name), rt_ms)
                    
                    # Return the row data to the main process
                    return validate_and_prepare_row(
                        metrics.get_package(), query, h, rt_ms, 
                        metrics.get_solutions(), seed
                    )
    except Exception as e:
        print(f"Error processing {full_path}: {e}")
        return None


def run_experiment(workload_directory, algorithm, M, N, relation_name):
    print(f"--- Running Parallel Experiment 90 Cores) ---")
    print(f"Algorithm: {algorithm}")
    print(f"Workload Base: {workload_directory}")
    
    # Define paths for Results
    results_base = os.path.join("results", algorithm.upper(), relation_name)
    csv_dir = os.path.join(results_base, "Results")
    os.makedirs(csv_dir, exist_ok=True)
    
    csv_path = os.path.join(csv_dir, f"{algorithm.upper()}_{N}_{M}.csv")
    is_seeded = algorithm.upper() in ['RCLSEED', 'SSSEED', 'RCLSEEDTIMEBUDGET']

    # --- 1. Prepare List of Tasks ---
    tasks = []
    if is_seeded:
        print("Detected Seeded Experiment. Generating task list...")
        for h in range(0, 9): # Hardness 0 to 8
            for seed in range(1, 11): # Seed 1 to 10
                filename = f"stocks_{N}_{M}_seeded_{seed}.spaql"
                full_path = os.path.join(workload_directory, str(h), filename)
                query_name = f"stocks_{N}_{M}_seeded_{seed}"
                
                if os.path.exists(full_path):
                    # Store tuple of args needed for process_single_task
                    tasks.append((full_path, query_name, h, seed, algorithm, M, results_base))
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
                        tasks.append((full_path, query_name, h, None, algorithm, M, results_base))
    
    print(f"Total Tasks Found: {len(tasks)}")

    # --- 2. Initialize CSV Writer ---
    with open(csv_path, 'a', newline='') as csvfile:
        if algorithm.upper() == 'HARDNESS':
            fieldnames = ['hardness', 'computed_hardness']
        else:
            fieldnames = ['Hardness', 'Objective', 'feas', 'surplus', 'Runtime', 'solutions']
            if is_seeded:
                fieldnames.insert(1, 'Seed')

        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        csvfile.flush() # <--- FORCE WRITE TO DISK IMMEDIATELY
        print(f"Header written to {csv_path}")

        # --- 3. Parallel Execution ---
        NUM_PROCESSES = 15
        print(f"Launching Pool with {NUM_PROCESSES} processes...")

        with Pool(processes=NUM_PROCESSES) as pool:
            # imap_unordered yields results as soon as they finish
            for result in pool.imap_unordered(process_single_task, tasks):
                if result:
                    writer.writerow(result)
                    csvfile.flush() # Ensure data is written immediately
                    
    print("Experiment Complete.")

if __name__ == '__main__':
    warnings.filterwarnings('ignore')
    parser = argparse.ArgumentParser()
    parser.add_argument('N', type=int)
    parser.add_argument('M', type=int)
    parser.add_argument('algorithm', type=str, choices=['RCL', 'SS', 'Naive', 'DETER', 'HARDNESS', 'RCLSeed', 'SSSeed', 'RCLSeedTimeBudget'])
    args = parser.parse_args()

    # --- Directory Logic ---
    if args.algorithm.upper() in ['RCLSEED', 'SSSEED', 'RCLSEEDTIMEBUDGET']:
        rel_name = f"stocks_{args.N}_{args.M}_seeded"
        workload_dir = os.path.join(WORKLOAD_BASE_SEED, rel_name, "RCL")
    else:
        rel_name = f"stocks_{args.N}_{args.M}"
        workload_dir = os.path.join(WORKLOAD_BASE_STD, rel_name, "RCL")

    if not os.path.exists(workload_dir):
        print(f"Error: Workload directory not found: {workload_dir}")
    else:
        run_experiment(workload_dir, args.algorithm, args.M, args.N, rel_name)