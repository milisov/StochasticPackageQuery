import numpy as np
import os
from numpy.random import SFC64, SeedSequence, Generator
from concurrent.futures import ThreadPoolExecutor
from PgConnection.PgConnection import PgConnection
from ScenarioGenerator.ScenarioGenerator import ScenarioGenerator
from Utils.Relation_Prefixes import Relation_Prefixes

SUBSTEPS = 10

def _process_ticker(ticker, group, seed, no_of_scenarios):
    _, last_row = group[-1]
    _, _, last_price, last_volatility, last_volatility_coeff, last_drift = last_row
    last_volatility_coeff = np.sqrt(last_volatility_coeff)

    sell_after_map = {int(row[1] * 2): idx for idx, row in group}
    sorted_dates = sorted(sell_after_map.keys())

    hashed_value = (seed + _hash(ticker)) % (10 ** 8)
    rng = Generator(SFC64(SeedSequence(hashed_value)))

    DRIFT_DAMPEN = 0.2
    exp_vol = last_volatility * last_volatility_coeff
    drift_coeff = (last_drift - 0.5 * exp_vol ** 2) * DRIFT_DAMPEN

    intervals = [d - prev_d for prev_d, d in zip([0] + sorted_dates[:-1], sorted_dates)]
    sub_step_sizes = [interval / SUBSTEPS for interval in intervals]

    expanded_dt = np.array(
        [sub_dt for sub_dt in sub_step_sizes for _ in range(SUBSTEPS)],
        dtype=np.float32
    )
    sqrt_dt = np.sqrt(expanded_dt)
    total_steps = len(expanded_dt)

    # Generate 20% more scenarios upfront
    extended_scenarios = int(no_of_scenarios * 1.2)
    Z = rng.standard_normal(size=(extended_scenarios, total_steps)).astype(np.float32)

    log_increments = drift_coeff * expanded_dt + exp_vol * sqrt_dt * Z
    cum_log_prices = np.cumsum(log_increments, axis=1)
    price_paths = last_price * np.exp(cum_log_prices)  # (extended_scenarios, total_steps)

    # Trim bottom and top 10% by final price, consistently across all timestamps
    final_prices = price_paths[:, -1]
    extended_scenarios = int(no_of_scenarios * 1.2)
    low  = (extended_scenarios - no_of_scenarios) // 2
    high = low + no_of_scenarios  # guarantees exactly no_of_scenarios paths
    sorted_indices = np.argsort(final_prices)[low:high]
    trimmed_paths = price_paths[sorted_indices]          # (remaining, total_steps)

    extraction_indices = [SUBSTEPS * (i + 1) - 1 for i in range(len(sorted_dates))]

    rng.shuffle(trimmed_paths)
    return {
        sell_after_map[date]: (trimmed_paths[:, extraction_indices[i]] - last_price).tolist()
        for i, date in enumerate(sorted_dates)
    }


def _hash(s: str) -> int:
    hashed_value = 0
    for char in s:
        hashed_value = hashed_value * 7727 + ord(char)
        hashed_value %= 2593697387
    return hashed_value


class GainScenarioGenerator(ScenarioGenerator):

    def __init__(self, relation: str, base_predicate: str = ''):
        self.__relation = relation
        self.__base_predicate = base_predicate or '1=1'

    def __get_info_partition(self, no_of_scenarios: int = None):
        partition_relation = self.__relation.split('_')[0] + '_' + self.__relation.split('_')[1] + '_partition'
        if no_of_scenarios is not None:
            sql_query = (
                f'SELECT profit[1:{no_of_scenarios}] FROM {partition_relation} '
                f'WHERE {self.__base_predicate} ORDER BY id;'
            )
        else:
            sql_query = (
                f'SELECT profit FROM {partition_relation} '
                f'WHERE {self.__base_predicate} ORDER BY id;'
            )
        PgConnection.Execute(sql_query)
        return PgConnection.Fetch()

    def __get_info(self):
        sql_query = (
            f'SELECT profit FROM {self.__relation} '
            f'WHERE {self.__base_predicate} ORDER BY id;'
        )
        PgConnection.Execute(sql_query)
        return PgConnection.Fetch()

    def generate_scenarios(self, seed: int, no_of_scenarios: int, no_of_duplicates = None) -> list[list[float]]:

        if no_of_duplicates is not None:
            # print(no_of_duplicates, no_of_scenarios)
            no_of_scenarios = no_of_scenarios * no_of_duplicates
            rows = self.__get_info_partition(no_of_scenarios)
        else:
            rows = self.__get_info()

        gains = []
        for row in rows:
            scenario_array = row[0]
            if len(scenario_array) > no_of_scenarios:
                gains.append(list(scenario_array[:no_of_scenarios]))
            else:
                gains.append(list(scenario_array))

        return gains

    def generate_scenarios_from_partition(
        self, seed: int, no_of_scenarios: int, partition_id: int
    ) -> list[list[float]]:
        self.__relation = (
            f'{self.__relation} AS r INNER JOIN '
            f'{Relation_Prefixes.PARTITION_RELATION_PREFIX}{self.__relation} '
            f'AS p ON r.id = p.tuple_id'
        )
        self.__base_predicate = (
            f'{self.__base_predicate} AND ' if self.__base_predicate != '1=1' else ''
        ) + f'p.partition_id = {partition_id}'

        return self.generate_scenarios(seed, no_of_scenarios)