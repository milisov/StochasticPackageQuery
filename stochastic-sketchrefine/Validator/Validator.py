from DbInfo import DbInfo
from Hyperparameters.Hyperparameters import Hyperparameters
from StochasticPackageQuery.Constraints.ExpectedSumConstraint.ExpectedSumConstraint import ExpectedSumConstraint
from StochasticPackageQuery.Constraints.VaRConstraint.VaRConstraint import VaRConstraint
from StochasticPackageQuery.Constraints.CVaRConstraint.CVaRConstraint import CVaRConstraint
from StochasticPackageQuery.Query import Query
from Utils.RelationalOperators import RelationalOperators
from Utils.TailType import TailType
from Utils.ObjectiveType import ObjectiveType
from ValueGenerator.ValueGenerator import ValueGenerator
from Utils.Stochasticity import Stochasticity
import math
import numpy as np


class Validator:

    def __init__(self,query: Query,
                 dbInfo: DbInfo,
                 no_of_validation_scenarios: int,
                 precomputed_scenarios=None,
                 precomputed_values=None):
        self.__query = query
        self.__dbInfo = dbInfo
        self.__no_of_validation_scenarios = \
            no_of_validation_scenarios
        self.__all_scenarios = precomputed_scenarios
        self.__all_values = precomputed_values


    def __get_stochastic_attributes(self):
        attributes = set()
        for constraint in self.__query.get_constraints():
            if constraint.is_expected_sum_constraint():
                attributes.add(
                    constraint.get_attribute_name())
            if constraint.is_risk_constraint():
                attributes.add(
                    constraint.get_attribute_name())
        if self.__query.get_objective().get_stochasticity() \
            == Stochasticity.STOCHASTIC:
            attributes.add(
                self.__query.get_objective().\
                    get_attribute_name()
            )
        return attributes
    

    def __get_scenarios_and_ids(self, package_dict: dict,
                              attribute: str):
        base_predicate = ''
        ids_with_multiplicities = []
        for id in package_dict:
            ids_with_multiplicities.append(
                (id, package_dict[id]))
            if len(base_predicate) > 0:
                base_predicate += " or "
            base_predicate += " id=" + str(id)
        ids_with_multiplicities.sort()
        scenarios = []
        if len(package_dict) > 0:
            if self.__all_scenarios is None:
                sc = ValueGenerator(
                        relation = self.__query.get_relation(),
                        base_predicate = base_predicate,
                        attribute = attribute
                        ).get_values()
                for s in sc:
                    scenarios.append(s[0])
            else:
                for id in package_dict:
                    i = id - 1 # get only the positive scenarios
                    scenarios.append(self.__all_scenarios[attribute][i])            


        return scenarios, ids_with_multiplicities


    def __get_deterministic_values_and_ids(self, package_dict: dict,
                                           attribute: str):
        ids_with_multiplicities = []
        for id_val in package_dict:
            ids_with_multiplicities.append((id_val, package_dict[id_val]))
        ids_with_multiplicities.sort()
        values = []
        if len(package_dict) > 0:
            if self.__all_values is not None and attribute in self.__all_values:
                for id_val, _ in ids_with_multiplicities:
                    values.append(self.__all_values[attribute][id_val - 1])
            else:
                base_predicate = ' or '.join(f'id={id_val}' for id_val, _ in ids_with_multiplicities)
                raw = ValueGenerator(
                    relation=self.__query.get_relation(),
                    base_predicate=base_predicate,
                    attribute=attribute
                ).get_values()
                for v in raw:
                    values.append(v[0])
        return values, ids_with_multiplicities


    def get_validation_objective_value(self, package_dict) -> float:
        if package_dict is None:
            if self.__query.get_objective().get_objective_type() == \
                ObjectiveType.MAXIMIZATION:
                return -math.inf
            else:
                return math.inf
        
        if len(package_dict) == 0:
            return 0
            
        attribute = \
            self.__query.get_objective().get_attribute_name()
        scenarios, ids_with_multiplicities = \
            self.__get_scenarios_and_ids(
                package_dict, attribute)
        idx = 0
        objective_value = 0
        # print(ids_with_multiplicities)
        for tuple_values in scenarios:
            _, multiplicity = ids_with_multiplicities[idx]
            idx += 1
            #print('Validation Average:', np.average(tuple_values))
            #print('Validation multiplicity:', multiplicity)
            objective_value += np.average(tuple_values)*multiplicity
        
        return objective_value
    

    def get_package_size_constraint_feasibility(
            self, package_dict, constraint, EPS=0.05) -> bool:
        if package_dict is None:
            return True
        size = float(sum(package_dict.values()))
        limit = float(constraint.get_package_size_limit())
        op = constraint.get_inequality_sign()
        if op == RelationalOperators.LESS_THAN_OR_EQUAL_TO:
            if size <= limit + EPS:
                return True
            else:
                return False
        elif op == RelationalOperators.GREATER_THAN_OR_EQUAL_TO:
            if size >= limit - EPS:
                return True
            else:
                return False
        else:
            if abs(size - limit) <= EPS:
                return True
            else:
                return False


    def get_deterministic_constraint_feasibility(
            self, package_dict, constraint, EPS=0.05) -> bool:
        if package_dict is None:
            return True
        det_sum = self.get_deterministic_sum_value(package_dict, constraint)
        limit = constraint.get_sum_limit()
        op = constraint.get_inequality_sign()
        if op == RelationalOperators.LESS_THAN_OR_EQUAL_TO:
            if det_sum <= limit + EPS:
                return True
            else:
                return False
        elif op == RelationalOperators.GREATER_THAN_OR_EQUAL_TO:
            if det_sum >= limit - EPS:
                return True
            else:
                return False
        else:
            if abs(det_sum - limit) <= EPS:
                return True
            else:
                return False


    def get_deterministic_sum_value(self, package_dict, constraint) -> float:
        if package_dict is None or len(package_dict) == 0:
            return 0.0
        attribute = constraint.get_attribute_name()
        values, ids_with_multiplicities = \
            self.__get_deterministic_values_and_ids(package_dict, attribute)
        det_sum = 0.0
        for i, (_, multiplicity) in enumerate(ids_with_multiplicities):
            det_sum += values[i] * multiplicity
        return det_sum


    def get_expected_sum_value(
            self, package_dict,
            expected_sum_constraint: ExpectedSumConstraint) -> float:
        if package_dict is None:
            return 0.0
        attribute = expected_sum_constraint.get_attribute_name()
        scenarios, ids_with_multiplicities = \
            self.__get_scenarios_and_ids(
                package_dict, attribute)
        idx = 0
        expected_sum = 0.0
        for scenario in scenarios:
            _, multiplicity = ids_with_multiplicities[idx]
            idx += 1
            expected_sum += np.sum(scenario) * multiplicity
        expected_sum /= self.__no_of_validation_scenarios
        return expected_sum


    def get_expected_sum_constraint_feasibility(
            self, package_dict,
            expected_sum_constraint: ExpectedSumConstraint,
            EPS=0.05) -> bool:
        if package_dict is None:
            return True

        expected_sum = self.get_expected_sum_value(
            package_dict, expected_sum_constraint)
        limit = expected_sum_constraint.get_sum_limit()
        op = expected_sum_constraint.get_inequality_sign()

        if op == RelationalOperators.GREATER_THAN_OR_EQUAL_TO:
            if expected_sum >= limit - EPS:
                return True
            else:
                return False

        if op == RelationalOperators.EQUALS:
            if abs(expected_sum - limit) <= EPS:
                return True
            else:
                return False

        if expected_sum <= limit + EPS:
            return True
        else:
            return False
    

    def get_var_among_validation_scenarios(
        self, package_dict,
        var_constraint: VaRConstraint
    ) -> float:
        # print("Package dict =", package_dict)
        attribute = var_constraint.get_attribute_name()
        scenarios, ids_with_multiplicities = \
            self.__get_scenarios_and_ids(
                package_dict, attribute
            )

        # print("Scenarios =", scenarios)
        # print("Ids with multiplicities =", ids_with_multiplicities)
        scenario_scores = []
        for _ in range(self.__no_of_validation_scenarios):
            idx = 0
            scenario_score = 0
            for __, multiplicity in ids_with_multiplicities:
                scenario_score += \
                    scenarios[idx][_] * multiplicity
                idx += 1
            scenario_scores.append(
                scenario_score
            )
        scenario_scores.sort()
        # print(scenario_scores)
        scenarios_to_consider = \
            int(np.floor((var_constraint.get_probability_threshold()*\
                      self.__no_of_validation_scenarios)))
        if var_constraint.get_inequality_sign() == RelationalOperators.GREATER_THAN_OR_EQUAL_TO:
            return scenario_scores[self.__no_of_validation_scenarios - scenarios_to_consider - 1]
        else:
            return scenario_scores[scenarios_to_consider]


    def get_var_constraint_satisfaction(
            self, package_dict,
            var_constraint: VaRConstraint
    ) -> float:
        EPS = 0.05
        if package_dict is None:
            return 1.00
        attribute = var_constraint.get_attribute_name()
        scenarios, ids_with_multiplicities = \
            self.__get_scenarios_and_ids(
                package_dict, attribute
            )
        satisfying_scenarios = 0
        for _ in range(self.__no_of_validation_scenarios):
            idx = 0
            scenario_score = 0
            for __, multiplicity in ids_with_multiplicities:
                scenario_score += \
                    scenarios[idx][_] * multiplicity
                idx += 1
            if var_constraint.get_inequality_sign() == \
                RelationalOperators.GREATER_THAN_OR_EQUAL_TO:
                if scenario_score >= var_constraint.get_sum_limit() - EPS:
                    #print("scenario score = ", scenario_score, "v =",var_constraint.get_sum_limit())
                    satisfying_scenarios += 1
            elif var_constraint.get_inequality_sign() == \
                RelationalOperators.LESS_THAN_OR_EQUAL_TO:
                if scenario_score <= var_constraint.get_sum_limit() + EPS:
                    satisfying_scenarios += 1
        print("Satisfying scenarios =", satisfying_scenarios, "out of", self.__no_of_validation_scenarios)
        return satisfying_scenarios / self.__no_of_validation_scenarios
    

    def get_cvar_constraint_satisfaction(
        self, package_dict,
        cvar_constraint: CVaRConstraint
    ) -> float:
        if package_dict is None:
            if cvar_constraint.get_inequality_sign() == \
                RelationalOperators.GREATER_THAN_OR_EQUAL_TO:
                return math.inf
            else:
                return -math.inf
        attribute = cvar_constraint.get_attribute_name()
        
        scenarios, ids_with_multiplicities = \
            self.__get_scenarios_and_ids(
                package_dict, attribute
            )
        scenario_scores = []
        for scenario_no in range(self.__no_of_validation_scenarios):
            idx = 0
            scenario_score = 0
            for __, multiplicity in ids_with_multiplicities:
                scenario_score += \
                    scenarios[idx][scenario_no] * multiplicity
                idx += 1
            scenario_scores.append(scenario_score)
        
        if cvar_constraint.get_tail_type() == TailType.HIGHEST:
            scenario_scores.sort(reverse=True)
        else:
            scenario_scores.sort()

        no_of_scenarios_to_consider = \
            int(
                np.floor(
                    self.__no_of_validation_scenarios*\
                    cvar_constraint.get_percentage_of_scenarios()\
                        /100
                )
            )

        return np.average(scenario_scores[
            0: no_of_scenarios_to_consider])    


    def get_var_constraint_feasibility(
        self, package_dict,
        var_constraint: VaRConstraint
    ) -> bool:
        if package_dict is None:
            return True
        probability = \
            self.get_var_constraint_satisfaction(
                package_dict, var_constraint
            )
        print("Probability =",probability, "p = ",var_constraint.get_probability_threshold())
        return (
            probability >= \
            var_constraint.get_probability_threshold()
        )

    def get_cvar_constraint_feasibility(
        self, package_dict,
        cvar_constraint: CVaRConstraint
    ) -> bool:
        print("CVaR feasibility")
        if package_dict is None:
            return True
        cvar = \
            self.get_cvar_constraint_satisfaction(
                package_dict, cvar_constraint,
            )
        
        if cvar_constraint.get_inequality_sign() == \
            RelationalOperators.LESS_THAN_OR_EQUAL_TO:
            return (cvar <= cvar_constraint.get_sum_limit())
        
        elif cvar_constraint.get_inequality_sign() == \
            RelationalOperators.GREATER_THAN_OR_EQUAL_TO:
            return (cvar >= cvar_constraint.get_sum_limit())
        
        return cvar == cvar_constraint.get_sum_limit()
    

    def validate_package(self, package_dict):
        deter_feasible = True
        prob_feasible = True
        feasibilities = []
        surpluses = []

        if package_dict:
            for constraint in self.__query.get_constraints():
                if constraint.is_package_size_constraint():
                    if not self.get_package_size_constraint_feasibility(package_dict, constraint):
                        deter_feasible = False
                elif constraint.is_deterministic_constraint():
                    if not self.get_deterministic_constraint_feasibility(package_dict, constraint):
                        deter_feasible = False
                elif constraint.is_expected_sum_constraint():
                    if not self.get_expected_sum_constraint_feasibility(package_dict, constraint):
                        deter_feasible = False
                elif constraint.is_var_constraint():
                    feas = self.get_var_constraint_satisfaction(package_dict, constraint)
                    p = constraint.get_probability_threshold()
                    feasibilities.append(feas)
                    surpluses.append(feas - p)
                    if feas < p:
                        prob_feasible = False
                elif constraint.is_cvar_constraint():
                    if not self.get_cvar_constraint_feasibility(package_dict, constraint):
                        prob_feasible = False

        objective_value = self.get_validation_objective_value(package_dict) if package_dict else 0.0
        return deter_feasible, prob_feasible, feasibilities, surpluses, objective_value


    def is_package_validation_feasible(
        self, package_dict,
    ):
        print("VALIDATING")
        if package_dict is None:
            print("NO SOLUTION FOR THAT ALPHA")
            return True
        
        for constraint in self.__query.get_constraints():
            if constraint.is_expected_sum_constraint():
                if not self.get_expected_sum_constraint_feasibility(
                    package_dict, constraint
                ):
                    return False
            if constraint.is_var_constraint():
                if not self.get_var_constraint_feasibility(
                    package_dict, constraint
                ):
                    return False
                
            if constraint.is_cvar_constraint():
                if not self.get_cvar_constraint_feasibility(
                    package_dict, constraint
                ):
                    return False
        return True


    def is_package_1_pm_epsilon_approximate(
        self, package_dict: dict,
        epsilon: float, upper_bound: float
    ):
        if package_dict is None:
            return False

        objective = \
            self.__query.get_objective()
        
        objective_value = \
            self.get_validation_objective_value(
                package_dict
            )
        if objective.get_objective_type() == \
            ObjectiveType.MAXIMIZATION:
            print("Objective Value", objective_value, "(1 - epsilon) * upper_bound", (1 - epsilon) * upper_bound)
            return objective_value >= \
                (1 - epsilon) * upper_bound
        
        return objective_value <= \
            (1 + epsilon) * upper_bound

    