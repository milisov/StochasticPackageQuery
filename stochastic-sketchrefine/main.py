from scenario_generation_demo import gen_scenarios
from CVaRification.CVaRification import CVaRification
from CVaRification.StaircaseCVaRification import StaircaseCVaRification
from CVaRification.RCLSolve import RCLSolve
from DbInfo.PortfolioInfo import PortfolioInfo
from DbInfo.TpchInfo import TpchInfo
from Naive.Naive import Naive
from SummarySearch.SummarySearch import SummarySearch
from OfflinePreprocessing.DistPartition import DistPartition
from PgConnection.PgConnection import PgConnection
# Assuming RCLSolveBasedHardness exists in a module, added an import for it.
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
import warnings
import time
import os
import numpy as np
import re
import argparse
import json
import csv

relation = "stocks_{}_{}"
relation_validation = "stocks_3_4_validate"
timeout_seconds = 10 * 60

json_solutions_directory = "/home/fm2288/StochasticPackageQuery/test/ExperimentsAUC/RobustSatisficing/Solutions/"
workload_base_directory = '/home/fm2288/StochasticPackageQuery/test/Queries'

def validate_and_store_results(package_dict: dict, query, hardness: int, runtime: float, writer: csv.DictWriter):
    """
    Validates a single solution and writes the results to a CSV file using a provided writer.

    Args:
        package_dict (dict): The solution package.
        query: The parsed query object.
        hardness (int): The hardness level of the query.
        runtime (float): The execution time of the algorithm in ms.
        writer (csv.DictWriter): The CSV writer object to use for storing results.
    """

    if "validate" in query.get_relation():
        query.set_relation(query.get_relation().replace("_validate", ""))

    validator = Validator(
        query=query,
        dbInfo=PortfolioInfo,
        no_of_validation_scenarios=10000
    )

    feasibility = 0
    for constraint in query.get_constraints():
        if constraint.is_risk_constraint():
            feasibility = validator.get_var_constraint_satisfaction(package_dict, constraint, True)

    objective_value = validator.get_validation_objective_value(package_dict, True)
    
    print(f"Validation Results - h={hardness}, feasibility: {feasibility}, Objective: {objective_value}")

    # Write the results to the pre-opened CSV file
    writer.writerow({
        'hardness': hardness,
        'objective': objective_value,
        'feas': feasibility,
        'runtime': runtime
    })

def run_experiment(workload_directory, algorithm, M, N, relation_name):
    """
    Runs experiments for a given algorithm and scenario configuration.

    Args:
        workload_directory (str): The directory containing the .spaql query files.
        algorithm (str): The algorithm to use.
        M (int): The M parameter for the experiment.
        N (int): The N parameter for the experiment.
        relation_name (str): The full name of the relation (e.g., 'stocks_3_4').
    """
    print(f"--- Running Experiment ---")
    print(f"Algorithm: {algorithm}, Relation: {relation_name}")
    print(f"Workload Directory: {workload_directory}")
    print("--------------------------")

    # --- Set up NEW results directory structure ---
    relation_results_dir = os.path.join("results", algorithm.upper(), relation_name)
    csv_results_dir = os.path.join(relation_results_dir, "Results")
    os.makedirs(csv_results_dir, exist_ok=True)
    
    print(f"JSON solution files will be stored in: {relation_results_dir}")
    print(f"CSV summary file will be stored in: {csv_results_dir}")

    # --- Set up CSV file and writer ---
    csv_filename = f"{algorithm.upper()}_{N}_{M}.csv"
    csv_filepath = os.path.join(csv_results_dir, csv_filename)
    
    with open(csv_filepath, 'a', newline='') as csvfile:
        # Conditionally set fieldnames based on the algorithm
        if algorithm.upper() == 'HARDNESS':
            fieldnames = ['hardness', 'computed_hardness']
        else:
            fieldnames = ['hardness', 'objective', 'feas', 'runtime']
            
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        print(f"Created CSV file for this run: {csv_filepath}")

        # --- Iterate through queries ---
        for file in sorted(os.listdir(workload_directory)):
            queryNameHardness = re.search(r"(stocks_\d_\d)_(.*\d+).spaql", file)
            if not queryNameHardness:
                continue
            
            queryName = queryNameHardness.group(1)
            h = int(queryNameHardness.group(2))
            if h < 8:
                continue
            
            full_query_path = os.path.join(workload_directory, file)
            print(f"\nProcessing query: {file} (Hardness: {h})")

            with open(full_query_path, 'r') as f:
                query = Parser().parse(f.readlines())
                SeedManager.reinitialize_seed()

                # Handle HARDNESS algorithm separately
                if algorithm.upper() == "HARDNESS":
                    hardness_evaluator = RCLSolveBasedHardness(
                        query=query, linear_relaxation=False,
                        dbInfo=PortfolioInfo,
                        init_no_of_scenarios=100,
                        no_of_validation_scenarios=10**M,
                        approximation_bound=0.05,
                        sampling_tolerance=1.00,
                        bisection_threshold=0.1
                    )
                    hardness_evaluator.solve()
                    computed_hardness = hardness_evaluator.get_model_probability()
                    print(f'Original Hardness = {h}, Computed Hardness = {computed_hardness}')
                    
                    # Write the hardness results to the CSV file
                    writer.writerow({
                        'hardness': h,
                        'computed_hardness': computed_hardness
                    })
                    
                # Logic for all other algorithms
                else:
                    solver = None
                    algo_params = {
                        'query': query,
                        'linear_relaxation': False,
                        'dbInfo': PortfolioInfo,
                        'init_no_of_scenarios': 100,
                        'no_of_validation_scenarios': 10**M,
                        'approximation_bound': 0.05
                    }

                    start_time = time.time()
                    if algorithm.upper() == 'RCL':
                        rcl_params = algo_params.copy()
                        rcl_params.update({'sampling_tolerance': 0.2, 'bisection_threshold': 0.1})
                        solver = RCLSolve(**rcl_params)
                    elif algorithm.upper() == 'SS':
                        ss_params = algo_params.copy()
                        ss_params['init_no_of_summaries'] = 1
                        solver = SummarySearch(**ss_params)
                    elif algorithm.upper() == 'NAIVE':
                        solver = Naive(**algo_params)
                    elif algorithm.upper() == 'DETER':
                        query.set_relation(query.get_relation() + "_validate")
                        ss_params = algo_params.copy()
                        ss_params.update({'init_no_of_scenarios': 10**4, 'no_of_validation_scenarios': 10**4})
                        ss_params['init_no_of_summaries'] = 1
                        solver = SummarySearch(**ss_params)
                    else:
                        raise ValueError(f"Unknown algorithm specified: {algorithm}")

                    if algorithm.upper() == 'RCL':
                        solver.solve(can_add_scenarios=True, start_time=start_time, timeout=timeout_seconds)
                    elif algorithm.upper() == 'DETER':
                        solver.solve(start_time=start_time, timeout=timeout_seconds, getDeterministic=True)
                    else:
                        solver.solve(start_time=start_time, timeout=timeout_seconds)
                        
                    end_time = time.time()
                    runtime_in_ms = (end_time - start_time) * 1000

                    print(f"Query {file} finished in {runtime_in_ms:.2f} ms.")

                    if hasattr(solver, 'get_metrics') and algorithm.upper() != "HARDNESS":
                        metrics = solver.get_metrics()
                        metrics.set_query(queryName)
                        metrics.set_hardness(h)

                        json_filename = f"{queryName}_{algorithm}_{h}.json"
                        full_json_path = os.path.join(relation_results_dir, json_filename)
                        metrics.log_to_json(full_json_path, runtime_in_ms)
                        print(f"Metrics logged to {full_json_path}")

                        solution_dict = metrics.get_package()
                        
                        validate_and_store_results(
                            package_dict=solution_dict,
                            query=query,
                            hardness=h,
                            runtime=runtime_in_ms,
                            writer=writer
                        )
                    else:
                        print(f"Note: '{algorithm}' solver has no 'get_metrics' method. No log or validation generated.")
            print("done")

if __name__ == '__main__':
    warnings.filterwarnings('ignore')

    parser = argparse.ArgumentParser(description="Run Stochastic Package Query solvers.")
    parser.add_argument('N', type=int, help='Parameter N for relation name.')
    parser.add_argument('M', type=int, help='Parameter M for relation name (and scenario count).')
    # Added 'HARDNESS' to the list of available choices
    parser.add_argument('algorithm', type=str, choices=['RCL', 'SS', 'Naive', 'DETER', 'HARDNESS'],
                        help='The algorithm to use for solving the queries.')

    args = parser.parse_args()

    relation_name = "stocks_{}_{}".format(args.N, args.M)
    workload_directory = os.path.join(workload_base_directory, relation_name, "RCL")

    print(f"Targeting workload directory: {workload_directory}")
    
    if not os.path.isdir(workload_directory):
        print(f"Error: Workload directory not found: {workload_directory}")
    else:
        run_experiment(
            workload_directory=workload_directory,
            algorithm=args.algorithm,
            M=args.M,
            N=args.N,
            relation_name=relation_name
        )