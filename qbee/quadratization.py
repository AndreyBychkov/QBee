import sys
import math
import sympy as sp
import pandas as pd
import multiprocessing as mp

from .structures import *
from .substitution_heuristics import get_heuristics, get_heuristic_sorter, get_heuristics_substitution_sorted
from tqdm.autonotebook import tqdm
from copy import deepcopy
from queue import Queue, Empty
from collections import deque
from typing import List, Optional, Callable
from .util import *


def quadratize(system: EquationSystem, search_algorithm: str = "ID-DLS", auxiliary_eq_type: str = "differential", heuristics: str = "default",
               initial_max_depth: int = 1, limit_depth: Optional[int] = None, debug: Optional[str] = None, log_file=None) -> QuadratizationResult:
    """
    Transforms the system into quadratic form using variable substitution technique.

    :param system: polynomial system of DAE
    :param search_algorithm: graph search algorithm for quadratization.
    :param auxiliary_eq_type: auxiliary equation form.
    :param heuristics: next substitution choice method.
    :param initial_max_depth: some search algorithms evaluate every systems where the number of substitutions does not exceed this number.
                              Put here your assumption of how long a chain of optimal substitutions might be.
    :param limit_depth: maximum number of substitutions. Raise error if no quadratic system is found within limit depth.
    :param debug: printing mode while quadratization is performed.
    :param log_file: output file for evaluation logging. Must be in 'csv' format.
    :returns: quadratic polynomial system

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
    **BFS**
        Heuristic Breadth-First Search
    **Best-First**
        Best-First Search: Heuristic Depth-First Search
    **ID-DLS**
        Iterative Deepening Depth-Limited Search with heuristics

    Debug
    ---------------
    **None** or **silent**
        display nothing;
    **info**
        display progress bars;

    """
    if not system.is_polynomial():
        raise RuntimeError("System is not polynomialized. Polynomialize it first.")

    disable_pbar = True if (debug is None or debug == 'silent') else False

    if limit_depth is None:
        limit_depth = sys.maxsize

    log_rows_list = list() if log_file else None
    algo = _quadratization_algorithms[search_algorithm]
    result = algo(system, auxiliary_eq_type, heuristics, initial_max_depth, limit_depth, disable_pbar, log_rows_list)

    if log_file:
        log_rows_list.append(
            {'from': system.equations_hash, 'name': result.system.equations_hash, 'substitution': list(result.substitutions)})
        log_df = pd.DataFrame(log_rows_list)
        log_df.to_csv(log_file, index=False)
    return result


def _best_first(system: EquationSystem, auxiliary_eq_type: str, heuristics: str, initial_max_depth: int, limit_depth: int,
                disable_pbar: bool, log_rows_list: Optional[List]) -> QuadratizationResult:
    # TODO(implement)
    pass


def _bfs(system: EquationSystem, auxiliary_eq_type: str, heuristics: str, initial_max_depth: int, limit_depth: int, disable_pbar: bool,
         log_rows_list: Optional[List]) -> QuadratizationResult:
    processed_systems_pbar = tqdm(unit="node", desc="Systems processed: ", position=0, disable=disable_pbar)
    queue_pbar = tqdm(unit="node", desc="Nodes in queue: ", position=1, disable=disable_pbar)

    if heuristics == "default":
        heuristics = "none"
    heuristic_sorter = get_heuristic_sorter(heuristics)

    system_queue = Queue()
    system_queue.put((system, list()), block=True)
    initial_eq_number = len(system.equations)
    statistics = EvaluationStatistics(depth=0, steps=0, method_name='BFS')

    quad_reached = False
    while not quad_reached:
        if system_queue.empty():
            raise RuntimeError("Limit depth passed. No quadratic system is found.")

        prev_system, substitution_chain = system_queue.get_nowait()
        if substitution_chain:
            last_substitution = substitution_chain[-1]
            curr_system = _make_new_system(prev_system, auxiliary_eq_type, last_substitution)

            if log_rows_list is not None:
                _log_append(log_rows_list, prev_system.equations_hash, curr_system.equations_hash, last_substitution)
        else:
            curr_system = prev_system

        curr_depth = len(curr_system.equations) - initial_eq_number

        queue_pbar.update(-1)
        statistics.steps += 1
        processed_systems_pbar.update(1)
        processed_systems_pbar.postfix = f"Current depth level: {curr_depth} / {limit_depth}"

        if curr_system.is_quadratic():
            statistics.depth = curr_depth
            refresh_and_close_progress_bars(processed_systems_pbar, queue_pbar)
            return QuadratizationResult(curr_system, statistics, tuple(substitution_chain))

        if curr_depth == limit_depth:
            continue

        possible_substitutions = heuristic_sorter(curr_system)
        for substitution in map(sp.Poly.as_expr, possible_substitutions):
            system_queue.put((curr_system, substitution_chain + [substitution]))
            queue_pbar.update(1)


def _iddls(system: EquationSystem, auxiliary_eq_type: str, heuristics: str, initial_max_depth: int, limit_depth: int, disable_pbar: bool,
           log_rows_list: Optional[List]) -> QuadratizationResult:
    processed_systems_pbar = tqdm(unit="node", desc="Systems processed: ", position=0, disable=disable_pbar)
    stack_pbar = tqdm(unit="node", desc="Nodes in queue: ", position=1, disable=disable_pbar)
    high_depth_stack_pbar = tqdm(unit="node", desc="Nodes in higher depth queue: ", position=2, disable=disable_pbar)

    heuristic_sorter = get_heuristic_sorter(heuristics)
    system_stack = deque()
    system_high_depth_stack = deque()

    curr_depth = 0
    curr_max_depth = initial_max_depth
    system_stack.append((system, curr_depth, list()))
    stack_pbar.update(1)

    statistics = EvaluationStatistics(depth=0, steps=0, method_name='ID-DLS')

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

        prev_system, curr_depth, substitution_chain = system_stack.popleft()
        if substitution_chain:
            last_substitution = substitution_chain[-1]
            curr_system = _make_new_system(prev_system, auxiliary_eq_type, last_substitution)

            if log_rows_list is not None:
                _log_append(log_rows_list, prev_system.equations_hash, curr_system.equations_hash, last_substitution)
        else:
            curr_system = prev_system

        stack_pbar.update(-1)
        stack_pbar.refresh()
        statistics.steps += 1
        processed_systems_pbar.update(1)
        processed_systems_pbar.postfix = f"Current max depth level: {curr_max_depth} / {limit_depth}"

        if curr_system.is_quadratic():
            statistics.depth = curr_depth
            refresh_and_close_progress_bars(processed_systems_pbar, stack_pbar, high_depth_stack_pbar)
            return QuadratizationResult(curr_system, statistics, tuple(substitution_chain))

        if curr_depth == limit_depth:
            continue

        possible_substitutions = heuristic_sorter(curr_system)[::-1]
        for substitution in map(sp.Poly.as_expr, possible_substitutions):
            if curr_depth < curr_max_depth:
                system_stack.appendleft((curr_system, curr_depth + 1, substitution_chain + [substitution]))
                stack_pbar.update(1)
            else:
                system_high_depth_stack.appendleft((curr_system, curr_depth + 1, substitution_chain + [substitution]))
                high_depth_stack_pbar.update(1)


def _mmdr(system: EquationSystem, auxiliary_eq_type: str, heuristics: str, initial_max_depth: int, limit_depth: int, disable_pbar: bool,
          log_rows_list: Optional[List]) -> QuadratizationResult:
    """
    Minimal Monomial Decomposition Reduction

    First choose a monomial with the least number of decompositions. Then select a substitution using heuristics.
    """
    processed_systems_pbar = tqdm(unit="node", desc="Systems processed: ", position=0, disable=disable_pbar)
    queue_pbar = tqdm(unit="node", desc="Nodes in queue: ", position=1, disable=disable_pbar)

    if heuristics == "default":
        heuristics = "none"
    heuristic_sorter = get_heuristic_sorter(heuristics)

    system_queue = Queue()
    system_queue.put((system, list()), block=True)
    initial_eq_number = len(system.equations)
    statistics = EvaluationStatistics(depth=0, steps=0, method_name='MMDR')

    quad_reached = False
    while not quad_reached:
        if system_queue.empty():
            raise RuntimeError("Limit depth passed. No quadratic system is found.")

        prev_system, substitution_chain = system_queue.get_nowait()
        if substitution_chain:
            last_substitution = substitution_chain[-1]
            curr_system = _make_new_system(prev_system, auxiliary_eq_type, last_substitution)

            if log_rows_list is not None:
                _log_append(log_rows_list, prev_system.equations_hash, curr_system.equations_hash, last_substitution)
        else:
            curr_system = prev_system

        curr_depth = len(curr_system.equations) - initial_eq_number

        queue_pbar.update(-1)
        statistics.steps += 1
        processed_systems_pbar.update(1)
        processed_systems_pbar.postfix = f"Current depth level: {curr_depth} / {limit_depth}"

        if curr_system.is_quadratic():
            statistics.depth = curr_depth
            refresh_and_close_progress_bars(processed_systems_pbar, queue_pbar)
            return QuadratizationResult(curr_system, statistics, tuple(substitution_chain))

        if curr_depth == limit_depth:
            continue

        minimal_decomposition = sorted(curr_system.get_monomial_decompositions(), key=lambda d: len(d))[0]
        substitutions = heuristic_sorter(curr_system, get_possible_substitutions_from_decompositions([minimal_decomposition, ]))
        for substitution in map(sp.Poly.as_expr, substitutions):
            system_queue.put_nowait((curr_system, curr_depth + 1, substitution_chain + [substitution]))
            queue_pbar.update(1)


def _id_mmdr(system: EquationSystem, auxiliary_eq_type: str, heuristics: str, initial_max_depth: int, limit_depth: int, disable_pbar: bool,
             log_rows_list: Optional[List]) -> QuadratizationResult:
    """
    Iterative Deepening Minimal Monomial Decomposition Reduction

    First choose a monomial with the least number of decompositions. Then select a substitution using heuristics.
    """
    processed_systems_pbar = tqdm(unit="node", desc="Systems processed: ", position=0, disable=disable_pbar)
    stack_pbar = tqdm(unit="node", desc="Nodes in queue: ", position=1, disable=disable_pbar)

    heuristic_sorter = get_heuristics_substitution_sorted(heuristics)

    curr_depth = 0
    system_stack = deque()
    system_stack.append((system, curr_depth, list()))

    statistics = EvaluationStatistics(depth=0, steps=0, method_name='ID-MMDR')

    quad_reached = False
    while not quad_reached:
        prev_system, curr_depth, substitution_chain = system_stack.pop()
        if substitution_chain:
            last_substitution = substitution_chain[-1]
            curr_system = _make_new_system(prev_system, auxiliary_eq_type, last_substitution)

            if log_rows_list is not None:
                _log_append(log_rows_list, prev_system.equations_hash, curr_system.equations_hash, last_substitution)
        else:
            curr_system = prev_system

        stack_pbar.update(-1)
        stack_pbar.refresh()
        statistics.steps += 1
        processed_systems_pbar.update(1)
        processed_systems_pbar.postfix = f"Current depth level: {curr_depth} / {limit_depth}"

        if curr_system.is_quadratic():
            statistics.depth = curr_depth
            refresh_and_close_progress_bars(processed_systems_pbar, stack_pbar)
            return QuadratizationResult(curr_system, statistics, tuple(substitution_chain))

        if curr_depth == limit_depth:
            continue

        used_substitutions = set()
        decompositions = curr_system.get_monomial_decompositions()
        for dec in sorted(decompositions, key=lambda d: len(d), reverse=True):
            substitutions = heuristic_sorter(curr_system, get_possible_substitutions_from_decompositions([dec, ]))
            for substitution in map(sp.Poly.as_expr, substitutions):
                if substitution not in used_substitutions:
                    system_stack.append((curr_system, curr_depth + 1, substitution_chain + [substitution]))
                    stack_pbar.update(1)
                    used_substitutions.add(substitution)


def _new_bfs(system: PolynomialSystem, auxiliary_eq_type="differential", limit_depth=None) -> QuadratizationResult:
    processed_systems_pbar = tqdm(unit="node", desc="Systems processed: ", position=0)
    queue_pbar = tqdm(unit="node", desc="Nodes in queue: ", position=1)

    system_queue = Queue()
    system_queue.put((system, list()), block=True)
    initial_eq_number = len(system.equations)
    statistics = EvaluationStatistics(depth=0, steps=0, method_name='New BFS')

    quad_reached = False
    while not quad_reached:
        if system_queue.empty():
            raise RuntimeError("Limit depth passed. No quadratic system is found.")

        prev_system, substitution_chain = system_queue.get_nowait()
        if substitution_chain:
            last_substitution = substitution_chain[-1]
            curr_system = _make_new_poly_system(prev_system, auxiliary_eq_type, last_substitution)
        else:
            curr_system = prev_system

        curr_depth = len(curr_system.equations) - initial_eq_number

        queue_pbar.update(-1)
        statistics.steps += 1
        processed_systems_pbar.update(1)
        processed_systems_pbar.postfix = f"Current depth level: {curr_depth} / {limit_depth}"

        if curr_system.is_quadratic([poly_to_monomial(sub) for sub in substitution_chain]):
            statistics.depth = curr_depth
            refresh_and_close_progress_bars(processed_systems_pbar, queue_pbar)
            return QuadratizationResult(curr_system, statistics, tuple(substitution_chain))

        if curr_depth == limit_depth:
            continue

        possible_substitutions = curr_system.get_possible_substitutions()
        for substitution in possible_substitutions:
            system_queue.put((curr_system, substitution_chain + [substitution]))
            queue_pbar.update(1)


def _make_new_system(system: EquationSystem, auxiliary_eq_type, substitution) -> EquationSystem:
    new_system = deepcopy(system)
    new_variable = new_system.variables.create_variable()
    equation_add_fun = new_system.auxiliary_equation_type_choose(auxiliary_eq_type)
    equation_add_fun(new_variable, substitution)
    return new_system


def _make_new_poly_system(system: PolynomialSystem, auxiliary_eq_type, substitution) -> PolynomialSystem:
    new_system = system.clone()
    new_variable = new_system.variables.create_variable()
    equation_add_fun = new_system.auxiliary_equation_inserter(auxiliary_eq_type)
    equation_add_fun(new_variable, substitution)
    return new_system


def _log_append(row_list: List, hash_before, hash_after, substitution):
    row_list.append({'from': hash_before, 'name': hash_after, 'substitution': substitution})


_quadratization_algorithms = \
    {
        "BFS"       : _bfs,
        "Best-First": _best_first,
        "ID-DLS"    : _iddls,
        "ID-MMDR"   : _id_mmdr,
    }
