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
from QueryHardness.HardnessEvaluator import HardnessEvaluator
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


if __name__ == '__main__':
    # warnings.filterwarnings('ignore')
    '''
    partitioner = DistPartition(
        relation='Stock_Investments_Half',
        dbInfo=PortfolioInfo
    )

    id_lists = ValueGenerator(
        relation='Stock_Investments_Half',
        base_predicate='',
        attribute='id'
    ).get_values()
    
    ids = []

    for tuple in id_lists:
        ids.append(tuple[0])
    
    partitioner.partition(ids)
    print('Number of partitions:', partitioner.get_no_of_partitions())
    '''
    id_lists = ValueGenerator(
        relation='stocks_3_2',
        base_predicate='',
        attribute='profit'
    ).get_values()
    
    ids = []

    for tuple in id_lists:
        ids.append(tuple[0])
    iter = 0
    workload_directory = '/home/fm2288/StochasticPackageQuery/test/Queries/stocks_4_2'
    for file in os.listdir(workload_directory):
        queryNameHardness = re.search(r"(stocks_\d_\d)_(.*).spaql", file)
        queryName = queryNameHardness.group(1)
        hardness = float(queryNameHardness.group(2))
        with open(
            workload_directory + '/' + file, 'r') as f:
            query = Parser().parse(f.readlines())
            SeedManager.reinitialize_seed()
            start_time = time.time()
            rclSolve = RCLSolve(
                query=query, linear_relaxation=False,
                dbInfo=PortfolioInfo,
                init_no_of_scenarios=100,
                no_of_validation_scenarios=100,
                approximation_bound=0.02,
                sampling_tolerance=0.2,
                bisection_threshold=0.01,
            )
        rclSolve.solve(can_add_scenarios = False)
        end_time = time.time()

        runtime = end_time - start_time
        runtime_in_ms = runtime * 1000

        rclMetrics = rclSolve.get_metrics()
        rclMetrics.set_query(queryName)
        rclMetrics.set_hardness(hardness)
        filename = queryName + ".json"
        rclMetrics.log_to_json(filename, runtime_in_ms)