from qbee import *
from examples import *
from copy import deepcopy
from functools import partial
from typing import List, Optional, Set, Dict
from qbee.quadratization import QuadratizationResult
from queue import Queue
from sympy.core.function import AppliedUndef

no_quad_pruning = partial(pruning_by_nodes_without_quadratization_found, nodes_processed=1000)


def pde(start_system: List[Tuple[sp.Symbol, sp.Expr]], input_orders: Optional[Dict] = None) -> Optional[
    QuadratizationResult]:
    queue = Queue()
    queue.put(start_system)
    if input_orders is None:
        inputs = select_inputs(start_system)
        input_orders = {i: 0 for i in inputs}
    while not queue.empty():
        system = queue.get_nowait()
        inputs_pde = select_pde_inputs(system)
        input_orders_with_pde = {i: 0 for i in inputs_pde}
        input_orders_with_pde.update(input_orders)
        print_common(inputs_pde)
        quad_res = polynomialize_and_quadratize(system, input_orders_with_pde,
                                                pruning_functions=[no_quad_pruning, *default_pruning_rules])
        if quad_res:
            return quad_res
        for i in inputs_pde:
            new_sys = deepcopy(system)
            ex, dx = rm_last_diff(i)
            try:
                new_sys.append((i, get_rhs(system, ex).diff(dx)))
                queue.put(new_sys)
            except AttributeError as e:
                pass
    return None


def get_rhs(system, sym):
    for lhs, rhs in system:
        if lhs == sym:
            return rhs
    return None


def rm_last_diff(der):
    if len(der.variables) == 1:
        return der.expr, der.variables[0]
    elif len(der.variables) > 1:
        return sp.Derivative(der.expr, *der.variables[:-1]), der.variables[-1]


def select_inputs(system):
    lhs, rhs = zip(*system)
    funcs = set(reduce(lambda l, r: l | r, [eq.atoms(AppliedUndef) for eq in rhs]))
    lhs_args = set(sp.flatten([eq.args for eq in lhs if not isinstance(eq, sp.Derivative)]))
    return set(filter(lambda f: (f not in lhs) and (f not in lhs_args), funcs))


def select_pde_inputs(system):
    lhs, rhs = zip(*system)
    return set(filter(lambda v: v not in lhs, reduce(lambda l, r: l | r, [eq.atoms(sp.Derivative) for eq in rhs])))


def scalar_pde_with_parameter_and_input():
    x, c = functions("x, c")
    u = multivariable_functions("u", [x])
    p = parameters("p")
    system = [
        (u, p * u.diff(x) ** 3 + c),
    ]

    res = pde(system)
    assert res is not None
    print(res)


if __name__ == '__main__':
    scalar_pde_with_parameter_and_input()
