import sys
import math
import sympy as sp
import pandas as pd
import multiprocessing as mp

from structures import *
from replacement_heuristics import get_heuristics, get_heuristic_sorter
from tqdm import tqdm
from copy import deepcopy
from queue import Queue, Empty
from collections import deque
from typing import List, Optional, Callable
from statistics import *


def quadratic_linearize(system: EquationSystem, mode: str = "optimal", auxiliary_eq_type: str = "differential", heuristics: str = "default",
                        method_optimal: str = "iddfs", initial_max_depth: int = 1, limit_depth: Optional[int] = None, debug: Optional[str] = None,
                        log_file=None) -> EquationSystem:
    """
    Transforms the system into quadratic-linear form using variable replacement technique.

    :param system: polynomialized system of DAE
    :param mode: type of search.
    :param auxiliary_eq_type: auxiliary equation form.
    :param heuristics: next replacement choice method.
    :param method_optimal: graph search algorithm in the optimal mode.
    :param initial_max_depth: for some methods checks all systems where the number of replacements does not exceed this number.
                              Put here your assumption of how long a chain of optimal replacements might be.
    :param limit_depth: maximum number of replacements. Raise error if no quadratic-linear system is found within limit depth.
    :param debug: printing mode while quadratic linearization is performed.
    :param log_file: output file for evaluation logging. Must be in 'csv' format.
    :returns: quadratic-linearized system

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
        choose next possible replacement in random way;
    **frequent-first**
        choose most frequent possible replacement as the next one;
    **free-variables-count**
        choose a replacement with the least free variables;
    **auxiliary-equation-degree**
        choose a monomial replacement with the least generated auxiliary equation degree;
    **auxiliary-equation-ql-discrepancy**
        choose a monomial replacement which generated auxiliary equation is closest to quadratic form;
    **summary-monomial-degree**
        choose a monomial replacement with maximal reduction of system's degree;

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
        prints replacement for each iteration;
    **debug**
        prints equations in system with replacement for each iteration;

    """
    if not system.is_polynomial():
        raise RuntimeError("System is not polynomialized. Polynomialize it first.")
    if mode == 'optimal':
        return _quadratic_linearize_optimal(system, auxiliary_eq_type, heuristics, method_optimal, initial_max_depth, limit_depth, debug, log_file)
    elif mode == 'heuristic':
        return _quadratic_linearize_heuristic(system, auxiliary_eq_type, heuristics, debug, log_file)
    else:
        raise ValueError("mode must be 'optimal' or 'heuristic'")


def _quadratic_linearize_heuristic(system: EquationSystem, auxiliary_eq_type: str, heuristics: str, debug: Optional[str] = None,
                                   log_file: Optional[str] = None) -> EquationSystem:
    log_rows_list = list()
    new_system = deepcopy(system)
    new_system.statistics = EvaluationStatistics(0, 0, heuristics)

    while not new_system.is_quadratic_linear():
        iter_fun = _heuristic_iter_choose(auxiliary_eq_type)
        hash_before, hash_after, replacement = iter_fun(new_system, heuristics)

        new_system._debug_system_print(debug)
        new_system.statistics.steps += 1
        if log_file:
            _ql_log_append(log_rows_list, hash_before, hash_after, replacement)

    if log_file:
        log_df = pd.DataFrame(log_rows_list)
        log_df.to_csv(log_file, index=False)

    if not (debug is None or debug == 'silent'):
        print('-' * 100)

    new_system.statistics.depth = len(new_system.equations) - len(system.equations)
    return new_system


def _heuristic_iter_choose(auxiliary_eq_type: str) -> Callable:
    if auxiliary_eq_type == 'differential':
        return _heuristic_differential_iter
    elif auxiliary_eq_type == 'algebraic':
        return _heuristic_algebraic_iter
    else:
        raise ValueError("auxiliary_eq_type must be 'algebraic' or 'differential'")


def _heuristic_differential_iter(system: EquationSystem, method: str):
    hash_before = system.equations_hash

    replacement = get_heuristics(method)(system)
    new_symbol, new_symbol_dot = system.variables.create_symbol_with_derivative()
    system.replace_monomial(replacement, new_symbol)
    system._replacement_equations.append(sp.Eq(new_symbol, replacement))

    system._equations.append(sp.Eq(new_symbol_dot, system._calculate_Lie_derivative(replacement)).expand())
    hash_after = system.equations_hash

    return hash_before, hash_after, replacement


def _heuristic_algebraic_iter(system: EquationSystem, method: str):
    raise NotImplementedError()


def _quadratic_linearize_optimal(system: EquationSystem, auxiliary_eq_type: str, heuristics: str = "default", method="iddfs",
                                 initial_max_depth: int = 1, limit_depth: Optional[int] = None, debug=None,
                                 log_file: Optional[str] = None) -> EquationSystem:
    disable_pbar = True if (debug is None or debug == 'silent') else False
    progress_bar = tqdm(total=1, unit='node', desc="System nodes queued: ", disable=disable_pbar)

    if limit_depth is None:
        limit_depth = sys.maxsize

    log_rows_list = list()
    solution = _optimal_method_choose(system, auxiliary_eq_type, heuristics, method, initial_max_depth, limit_depth, progress_bar, log_rows_list)

    if log_file:
        log_df = pd.DataFrame(log_rows_list)
        log_df.to_csv(log_file, index=False)
    return solution


def _optimal_method_choose(system: EquationSystem, auxiliary_eq_type, heuristics, method: str, initial_max_depth, limit_depth, progress_bar,
                           log_rows_list):
    if method == "bfs":
        return _optimal_bfs(system, auxiliary_eq_type, limit_depth, progress_bar, log_rows_list)
    elif method == "iddfs":
        return _optimal_iddfs(system, auxiliary_eq_type, heuristics, initial_max_depth, limit_depth, progress_bar, log_rows_list)
    else:
        raise ValueError("Optimal method must be 'bfs' of 'iddfs'")


def _optimal_bfs(system: EquationSystem, auxiliary_eq_type: str, limit_depth, progress_bar: tqdm, log_rows_list: Optional[List]) -> EquationSystem:
    system_queue = Queue()
    system_queue.put(system, block=True)
    initial_eq_number = len(system.equations)
    statistics = EvaluationStatistics(depth=0, steps=0, method_name='BFS')

    ql_reached = False
    while not ql_reached:
        if system_queue.empty():
            raise RuntimeError("Limit depth passed. No quadratic-linear system is found.")
        curr_system = system_queue.get_nowait()
        curr_depth = len(curr_system.equations) - initial_eq_number

        progress_bar.update(-1)
        progress_bar.postfix = f"Current depth level: {curr_depth}"
        statistics.steps += 1
        progress_bar.total -= 1

        if curr_system.is_quadratic_linear():
            statistics.depth = curr_depth
            curr_system.statistics = statistics
            progress_bar.close()
            return curr_system

        if curr_depth == limit_depth:
            continue

        possible_replacements = curr_system.get_possible_replacements()
        progress_bar.total += len(possible_replacements)
        for replacement in map(sp.Poly.as_expr, possible_replacements):
            new_system = deepcopy(curr_system)
            new_symbol = new_system.variables.create_symbol()
            equation_add_fun = new_system.auxiliary_equation_type_choose(auxiliary_eq_type)
            equation_add_fun(new_symbol, replacement)

            system_queue.put(new_system)
            progress_bar.update(1)
            if log_rows_list is not None:
                _ql_log_append(log_rows_list, curr_system.equations_hash, new_system.equations_hash, replacement)


def _optimal_bfs_parallel(system: EquationSystem, auxiliary_eq_type: str):
    system_queue = Queue()
    system_queue.put_nowait(system)
    initial_eq_number = len(system.equations)

    ql_reached = False
    while not ql_reached:
        curr_system = system_queue.get()

        if curr_system.is_quadratic_linear():
            return curr_system

        possible_replacements = curr_system.get_possible_replacements()
        for replacement in map(sp.Poly.as_expr, possible_replacements):
            new_system = deepcopy(curr_system)
            new_symbol = new_system.variables.create_symbol()
            equation_add_fun = new_system._auxiliary_equation_type_choose(auxiliary_eq_type)
            equation_add_fun(new_symbol, replacement)

            system_queue.put(new_system)


def _optimal_bfs_parallel_worker(system: EquationSystem, system_queue: Queue, auxiliary_eq_type: str):
    while True:
        try:
            curr_system = system_queue.get_nowait()
        except Empty:
            break
        else:
            pass


def _optimal_iddfs(system: EquationSystem, auxiliary_eq_type: str, heuristics: str, initial_max_depth: int, limit_depth: int, progress_bar: tqdm,
                   log_rows_list: Optional[List]) -> EquationSystem:
    heuristic_sorter = get_heuristic_sorter(heuristics)
    system_stack = deque()
    system_high_depth_stack = deque()

    initial_eq_number = len(system.equations)
    curr_depth = 0
    curr_max_depth = initial_max_depth
    system_stack.append((system, curr_depth))

    statistics = EvaluationStatistics(depth=0, steps=0, method_name='ID-DFS')

    ql_reached = False
    while not ql_reached:
        if len(system_stack) == 0:
            if len(system_high_depth_stack) == 0:
                raise RuntimeError("Limit depth passed. No quadratic-linear system is found.")
            system_stack = system_high_depth_stack
            system_high_depth_stack = deque()
            curr_max_depth += int(math.ceil(math.log(curr_depth + 1)))

        curr_system, curr_depth = system_stack.pop()

        progress_bar.update(-1)
        progress_bar.postfix = f"Current max depth level: {curr_max_depth}"
        statistics.steps += 1
        progress_bar.total -= 1

        if curr_system.is_quadratic_linear():
            progress_bar.close()
            statistics.depth = curr_depth
            curr_system.statistics = statistics
            return curr_system

        if curr_depth == limit_depth:
            continue

        possible_replacements = heuristic_sorter(curr_system)
        progress_bar.total += len(possible_replacements)
        for replacement in map(sp.Poly.as_expr, possible_replacements[::-1]):
            new_system = deepcopy(curr_system)
            new_symbol = new_system.variables.create_symbol()
            equation_add_fun = new_system.auxiliary_equation_type_choose(auxiliary_eq_type)
            equation_add_fun(new_symbol, replacement)

            progress_bar.update(1)
            if curr_depth < curr_max_depth:
                system_stack.append((new_system, curr_depth + 1))
            else:
                system_high_depth_stack.append((new_system, curr_depth + 1))

            if log_rows_list is not None:
                _ql_log_append(log_rows_list, curr_system.equations_hash, new_system.equations_hash, replacement)




def _ql_log_append(row_list: List, hash_before, hash_after, replacement):
    row_list.append({'from': hash_before, 'name': hash_after, 'replacement': replacement})
