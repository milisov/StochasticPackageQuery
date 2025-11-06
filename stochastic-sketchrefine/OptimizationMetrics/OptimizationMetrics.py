import time
import json
import os

class OptimizationMetrics:

    def __init__(
        self,
        algorithm_name: str,
        linear_relaxation: bool,
    ) -> None:
        self.__algorithm_name = \
            algorithm_name
        self.__linear_relaxation = \
            linear_relaxation
        self.__runtime = 0
        self.__objective_value = 0
        self.__optimizer_starting_time = None
        self.__optimizer_runtime = 0
        self.__number_of_optimization_calls = 0
        self.__number_of_scenarios_needed = 0
        self.__starting_time = None
        self.__package = []
        self.__query = None
        self.__hardness = None
        self.__p = None
    
    def get_optimizer_runtime(self):
        return self.__optimizer_runtime

    def start_execution(self):
       self.__starting_time = time.time()

    def get_package(self):
        return self.__package

    def end_execution(
        self, objective_value, no_of_scenarios
    ):
        assert self.__starting_time is not None
        self.__runtime = \
            time.time() - self.__starting_time
        self.__objective_value = \
            objective_value
        print("Metrics Objective Value = ", objective_value)
        self.__number_of_scenarios_needed = \
            no_of_scenarios


    def start_optimizer(self):
        self.__optimizer_starting_time = \
            time.time()


    def end_optimizer(self):
        assert self.__optimizer_starting_time\
            is not None
        self.__optimizer_runtime += time.time()\
            - self.__optimizer_starting_time
        print("Current Optimizer Runtime:",self.__optimizer_runtime)
        self.__number_of_optimization_calls += 1

    def set_package(self, res_package):
        print("Metrics Res Package = ", res_package)
        self.__package = res_package

    def set_query(self, q):
        self.__query = q

    def set_hardness(self, h):
        self.__hardness = h
    
    def set_p(self, p):
        self.__p = p

    def log(self):
        print('Algorithm:',
              self.__algorithm_name)
        print('Linear Relaxation:',
              self.__linear_relaxation)
        print('Runtime:',
              self.__runtime)
        print('Objective Value:',
              self.__objective_value)
        print('Number of optimization calls:', 
              self.__number_of_optimization_calls)
        print('Total Optimizer Runtime:',
              self.__optimizer_runtime)
        print('Number of scenarios needed:',
              self.__number_of_scenarios_needed)

    def log_to_json(self, filename="results.json", runtime_in_ms = None):
        print(runtime_in_ms)
        data = {
            "Algorithm": self.__algorithm_name,
            "Query": self.__query,
            "Hardness": self.__hardness,
            "p": self.__p,
            "Linear Relaxation": self.__linear_relaxation,
            "Runtime": self.__runtime,
            "Package": self.__package,
            "Objective Value": self.__objective_value,
            "Number of Optimization Calls": self.__number_of_optimization_calls,
            "Total Optimizer Runtime": self.__optimizer_runtime,
            "Number of Scenarios Needed": self.__number_of_scenarios_needed,
            "Runtime": runtime_in_ms
        }

        # Check if file exists and read existing data
        if os.path.exists(filename):
            with open(filename, "r") as json_file:
                try:
                    existing_data = json.load(json_file)  # Load existing JSON
                    if not isinstance(existing_data, list):
                        existing_data = []  # Ensure it's a list
                except json.JSONDecodeError:
                    existing_data = []  # Handle case where file is empty or invalid
        else:
            existing_data = []

        # Append new data
        existing_data.append(data)

        # Write back to file in valid JSON format
        with open(filename, "w") as json_file:
            json.dump(existing_data, json_file, indent=4)