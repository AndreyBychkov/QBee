import sys
import math
import sympy as sp
import pandas as pd
import multiprocessing as mp

from .structures import *
from .substitution_heuristics import get_heuristics, get_heuristic_sorter
from tqdm.autonotebook import tqdm
from copy import deepcopy
from queue import Queue, Empty
from collections import deque
from typing import List, Optional, Callable
from .util import reset_progress_bar, refresh_and_close_progress_bars


def quadratize(system: EquationSystem, mode: str = "optimal", auxiliary_eq_type: str = "differential", heuristics: str = "default",
               method_optimal: str = "iddfs", initial_max_depth: int = 1, limit_depth: Optional[int] = None, debug: Optional[str] = None,
               log_file=None) -> QuadratizationResult:
    """
    Transforms the system into quadratic form using variable substitution technique.

    :param system: polynomial system of DAE
    :param mode: type of search.
    :param auxiliary_eq_type: auxiliary equation form.
    :param heuristics: next substitution choice method.
    :param method_optimal: graph search algorithm in the optimal mode.
    :param initial_max_depth: for some methods checks all systems where the number of substitutions does not exceed this number.
                              Put here your assumption of how long a chain of optimal substitutions might be.
    :param limit_depth: maximum number of substitutions. Raise error if no quadratic system is found within limit depth.
    :param debug: printing mode while quadratization is performed.
    :param log_file: output file for evaluation logging. Must be in 'csv' format.
    :returns: quadratic polynomial system

    Mode
    -----------------
    **optimal**
        find optimal transformation. The most time-consuming mode;
    **heuristic**
        find sub-optimal transformation. Works much faster than 'optimal'. You can choose heuristics in 'heuristics' parameter;

    Auxiliary equations type
    -----------------
    **algebraic**
        adds auxiliary equations in form y = f(x, y)
    **differential**
         adds auxiliary equations in form y' = f(x, y)

    Heuristics
    -----------------
    **random**
        choose next possible substitution in random way;
    **frequent-first**
        choose most frequent possible substitution as the next one;
    **free-variables-count**
        choose a substitution with the least free variables;
    **auxiliary-equation-degree**
        choose a monomial substitution with the least generated auxiliary equation degree;
    **auxiliary-equation-quadratic-discrepancy**
        choose a monomial substitution which generated auxiliary equation is closest to quadratic form;
    **summary-monomial-degree**
        choose a monomial substitution with maximal reduction of system's degree;

    Method Optimal
    -----------------
    **bfs**
        Breadth-First Search
    **iddfs**
        Iterative Deepening Depth-First Search

    Debug
    ---------------
    **None** or **silent**
        prints nothing;
    **info**
        prints substitution for each iteration;
    **debug**
        prints equations in system with substitution for each iteration;

    """
    if not system.is_polynomial():
        raise RuntimeError("System is not polynomialized. Polynomialize it first.")
    if mode == 'optimal':
        return _quadratize_optimal(system, auxiliary_eq_type, heuristics, method_optimal, initial_max_depth, limit_depth, debug, log_file)
    elif mode == 'heuristic':
        return _quadratize_heuristic(system, auxiliary_eq_type, heuristics, debug, log_file)
    else:
        raise ValueError("mode must be 'optimal' or 'heuristic'")


def _quadratize_heuristic(system: EquationSystem, auxiliary_eq_type: str, heuristics: str, debug: Optional[str] = None,
                          log_file: Optional[str] = None) -> QuadratizationResult:
    log_rows_list = list()
    new_system = deepcopy(system)
    statistics = EvaluationStatistics(0, 0, heuristics)
    substitutions = list()

    while not new_system.is_quadratic():
        iter_fun = _heuristic_iter_choose(auxiliary_eq_type)
        hash_before, hash_after, substitution = iter_fun(new_system, heuristics)

        new_system._debug_system_print(debug)
        substitutions.append(substitution)
        statistics.steps += 1
        if log_file:
            _quad_log_append(log_rows_list, hash_before, hash_after, substitution)

    if log_file:
        log_df = pd.DataFrame(log_rows_list)
        log_df.to_csv(log_file, index=False)

    if not (debug is None or debug == 'silent'):
        print('-' * 100)

    statistics.depth = len(new_system.equations) - len(system.equations)
    return QuadratizationResult(new_system, statistics, tuple(substitutions))


def _heuristic_iter_choose(auxiliary_eq_type: str) -> Callable:
    if auxiliary_eq_type == 'differential':
        return _heuristic_differential_iter
    elif auxiliary_eq_type == 'algebraic':
        return _heuristic_algebraic_iter
    else:
        raise ValueError("auxiliary_eq_type must be 'algebraic' or 'differential'")


def _heuristic_differential_iter(system: EquationSystem, method: str):
    hash_before = system.equations_hash

    substitution = get_heuristics(method)(system)
    new_variable, new_variable_dot = system.variables.create_variable_with_derivative()
    system.replace_monomial(substitution, new_variable)
    system._substitution_equations.append(sp.Eq(new_variable, substitution))

    system._equations.append(sp.Eq(new_variable_dot, system._calculate_Lie_derivative(substitution)).expand())
    hash_after = system.equations_hash

    return hash_before, hash_after, substitution


def _heuristic_algebraic_iter(system: EquationSystem, method: str):
    raise NotImplementedError()


def _quadratize_optimal(system: EquationSystem, auxiliary_eq_type: str, heuristics: str = "default", method="iddfs",
                        initial_max_depth: int = 1, limit_depth: Optional[int] = None, debug=None,
                        log_file: Optional[str] = None) -> QuadratizationResult:
    disable_pbar = True if (debug is None or debug == 'silent') else False

    if limit_depth is None:
        limit_depth = sys.maxsize

    log_rows_list = list()
    result = _optimal_method_choose(system, auxiliary_eq_type, heuristics, method, initial_max_depth, limit_depth, disable_pbar, log_rows_list)

    if log_file:
        log_rows_list.append(
            {'from': system.equations_hash, 'name': result.system.equations_hash, 'substitution': list(result.substitutions)})
        log_df = pd.DataFrame(log_rows_list)
        log_df.to_csv(log_file, index=False)
    return result


def _optimal_method_choose(system: EquationSystem, auxiliary_eq_type, heuristics, method: str, initial_max_depth, limit_depth, disable_pbar,
                           log_rows_list) -> QuadratizationResult:
    if method == "bfs":
        return _optimal_bfs(system, auxiliary_eq_type, limit_depth, disable_pbar, log_rows_list)
    elif method == "iddfs":
        return _optimal_iddfs(system, auxiliary_eq_type, heuristics, initial_max_depth, limit_depth, disable_pbar, log_rows_list)
    else:
        raise ValueError("Optimal method must be 'bfs' of 'iddfs'")


def _optimal_bfs(system: EquationSystem, auxiliary_eq_type: str, limit_depth, disable_pbar: tqdm,
                 log_rows_list: Optional[List]) -> QuadratizationResult:
    processed_systems_pbar = tqdm(unit="node", desc="Systems processed: ", position=0, disable=disable_pbar)
    queue_pbar = tqdm(unit="node", desc="Nodes in queue: ", position=1, disable=disable_pbar)

    substitution_chains = set()

    system_queue = Queue()
    system_queue.put((system, list()), block=True)
    initial_eq_number = len(system.equations)
    statistics = EvaluationStatistics(depth=0, steps=0, method_name='BFS')

    quad_reached = False
    while not quad_reached:
        if system_queue.empty():
            raise RuntimeError("Limit depth passed. No quadratic system is found.")
        curr_system, curr_substitutions = system_queue.get_nowait()
        curr_depth = len(curr_system.equations) - initial_eq_number

        queue_pbar.update(-1)
        statistics.steps += 1
        processed_systems_pbar.update(1)
        processed_systems_pbar.postfix = f"Current depth level: {curr_depth} / {limit_depth}"

        if curr_system.is_quadratic():
            statistics.depth = curr_depth
            refresh_and_close_progress_bars(processed_systems_pbar, queue_pbar)
            return QuadratizationResult(curr_system, statistics, tuple(curr_substitutions))

        if curr_depth == limit_depth:
            continue

        possible_substitutions = curr_system.get_possible_substitutions()
        for substitution in map(sp.Poly.as_expr, possible_substitutions):
            supplemented_substitutions = curr_substitutions + [substitution]
            supplemented_substitutions_set = frozenset(supplemented_substitutions)
            if supplemented_substitutions_set in substitution_chains:
                continue

            new_system = _make_new_system(curr_system, auxiliary_eq_type, substitution)
            system_queue.put((new_system, supplemented_substitutions))
            substitution_chains.add(supplemented_substitutions_set)
            queue_pbar.update(1)

            if log_rows_list is not None:
                _quad_log_append(log_rows_list, curr_system.equations_hash, new_system.equations_hash, substitution)


def _optimal_bfs_parallel(system: EquationSystem, auxiliary_eq_type: str):
    system_queue = Queue()
    system_queue.put_nowait(system)
    initial_eq_number = len(system.equations)

    quad_reached = False
    while not quad_reached:
        curr_system = system_queue.get()

        if curr_system.is_quadratic():
            return curr_system

        possible_substitutions = curr_system.get_possible_substitutions()
        for substitution in map(sp.Poly.as_expr, possible_substitutions):
            new_system = deepcopy(curr_system)
            new_variable = new_system.free.create_variable()
            equation_add_fun = new_system._auxiliary_equation_type_choose(auxiliary_eq_type)
            equation_add_fun(new_variable, substitution)

            system_queue.put(new_system)


def _optimal_bfs_parallel_worker(system: EquationSystem, system_queue: Queue, auxiliary_eq_type: str):
    while True:
        try:
            curr_system = system_queue.get_nowait()
        except Empty:
            break
        else:
            pass


def _optimal_iddfs(system: EquationSystem, auxiliary_eq_type: str, heuristics: str, initial_max_depth: int, limit_depth: int, disable_pbar: tqdm,
                   log_rows_list: Optional[List]) -> QuadratizationResult:
    processed_systems_pbar = tqdm(unit="node", desc="Systems processed: ", position=0, disable=disable_pbar)
    stack_pbar = tqdm(unit="node", desc="Nodes in queue: ", position=1, disable=disable_pbar)
    high_depth_stack_pbar = tqdm(unit="node", desc="Nodes in higher depth queue: ", position=2, disable=disable_pbar)

    heuristic_sorter = get_heuristic_sorter(heuristics)
    system_stack = deque()
    system_high_depth_stack = deque()

    substitution_chains = set()

    curr_depth = 0
    curr_max_depth = initial_max_depth
    system_stack.append((system, curr_depth, list()))
    stack_pbar.update(1)

    statistics = EvaluationStatistics(depth=0, steps=0, method_name='ID-DFS')

    quad_reached = False
    while not quad_reached:
        if len(system_stack) == 0:
            if len(system_high_depth_stack) == 0:
                raise RuntimeError("Limit depth passed. No quadratic system is found.")
            system_stack = system_high_depth_stack
            system_high_depth_stack = deque()
            curr_max_depth += int(math.ceil(math.log(curr_depth + 1)))

            reset_progress_bar(stack_pbar, len(system_stack))
            reset_progress_bar(high_depth_stack_pbar, 0)

        curr_system, curr_depth, curr_substitutions = system_stack.popleft()

        stack_pbar.update(-1)
        stack_pbar.refresh()
        statistics.steps += 1
        processed_systems_pbar.update(1)
        processed_systems_pbar.postfix = f"Current max depth level: {curr_max_depth} / {limit_depth}"

        if curr_system.is_quadratic():
            statistics.depth = curr_depth
            refresh_and_close_progress_bars(processed_systems_pbar, stack_pbar, high_depth_stack_pbar)
            return QuadratizationResult(curr_system, statistics, tuple(curr_substitutions))

        if curr_depth == limit_depth:
            continue

        possible_substitutions = heuristic_sorter(curr_system)
        for substitution in map(sp.Poly.as_expr, possible_substitutions):
            supplemented_substitutions = curr_substitutions + [substitution]
            supplemented_substitutions_set = frozenset(supplemented_substitutions)
            if supplemented_substitutions_set in substitution_chains:
                continue

            new_system = _make_new_system(curr_system, auxiliary_eq_type, substitution)
            substitution_chains.add(supplemented_substitutions_set)

            if curr_depth < curr_max_depth:
                system_stack.append((new_system, curr_depth + 1, supplemented_substitutions))
                stack_pbar.update(1)
            else:
                system_high_depth_stack.append((new_system, curr_depth + 1, supplemented_substitutions))
                high_depth_stack_pbar.update(1)

            if log_rows_list is not None:
                _quad_log_append(log_rows_list, curr_system.equations_hash, new_system.equations_hash, substitution)


def _make_new_system(system: EquationSystem, auxiliary_eq_type, substitution) -> EquationSystem:
    new_system = deepcopy(system)
    new_variable = new_system.variables.create_variable()
    equation_add_fun = new_system.auxiliary_equation_type_choose(auxiliary_eq_type)
    equation_add_fun(new_variable, substitution)
    return new_system


def _quad_log_append(row_list: List, hash_before, hash_after, substitution):
    row_list.append({'from': hash_before, 'name': hash_after, 'substitution': substitution})
