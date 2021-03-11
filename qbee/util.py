import pickle
import pandas as pd
import sympy as sp
import numpy as np
from sympy.polys.rings import PolyElement
from scipy.spatial import ConvexHull
from time import time
from tqdm import tqdm
from typing import Iterable, Collection, Union, Tuple, List
from sympy.polys.monomials import monomial_deg
from itertools import combinations, product
from functools import reduce, wraps
from operator import add
from copy import deepcopy


def parametrized(dec):
    def layer(*args, **kwargs):
        def repl(f):
            return dec(f, *args, **kwargs)

        return repl

    return layer


@parametrized
def timed(func, enabled=True):
    def wrapper(*args, **kwargs):
        start_time = time()
        res = func(*args, **kwargs)
        end_time = time()
        print()
        print(f"Elapsed time: {np.round(end_time - start_time, 3)}s.")
        return res

    if enabled:
        return wrapper
    else:
        return func


__log_pb_evaluated = False
__log_pb = None
__log_records = list()


@parametrized
def progress_bar(func, is_stop, enabled=True):
    def wrapper(*args, **kwargs):
        global __log_pb_evaluated
        global __log_pb
        if not __log_pb_evaluated:
            __log_pb_evaluated = True
            __log_pb = tqdm(desc="Nodes processed", unit=" nodes")
        if is_stop:
            __log_pb.close()
            __log_pb = None
            __log_pb_evaluated = False
        else:
            __log_pb.update(1)
            __log_pb.refresh()
        return func(*args, **kwargs)

    if enabled:
        return wrapper
    else:
        return func


@parametrized
def logged(method, enabled, log_file, is_stop=False):
    def wrapper(self, *args, **kwargs):
        global __log_records
        res = method(self, *args, **kwargs)
        if is_stop:
            df = pd.DataFrame(__log_records, columns=['from', 'to']).applymap(str)
            df.to_feather(log_file)
            __log_records = list()
        else:
            part_res, *_ = args
            for r in list(res):
                __log_records.append([part_res, r])
        return res

    if enabled:
        return wrapper
    else:
        return method


@parametrized
def dump_results(method, enabled, log_file):
    def wrapper(self, *args, **kwargs):
        res = method(self, *args, **kwargs)
        pickle.dump(res, open(log_file, "wb"))
        return res

    if enabled:
        return wrapper
    else:
        return method


@parametrized
def memoize_first(func, max_size):
    memo = {}

    def wrapper(*args):
        if args in memo:
            return memo[args]
        else:
            rv = func(*args)
            if len(memo) <= max_size:
                memo[args] = rv
            return rv

    return wrapper


def derivatives(names) -> Union[sp.Symbol, Tuple[sp.Symbol]]:
    """
    Add dot to input symbols.

    :param names: input symbols. Can be represented in different ways.
    :return: Input symbols wrapper with dots.

    Example:
        .. math:: \dot{x}, \dot{y} = derivatives([x, y])

    Names variants
    -----------------
    **symbol**
        derivatives('x'), derivatives(x: Symbol)
    **string with delimiters**
        derivatives('x, y, z'), derivatives('x y z')
    **iterable of symbols**
        derivatives([x, y, z]), derivatives(['x', 'y', 'z'])
    """
    if not isinstance(names, Iterable) and isinstance(names, sp.Symbol):
        return (make_derivative_symbol(names),)

    if isinstance(names, Iterable) and reduce(lambda a, b: a and b,
                                              map(lambda elem: isinstance(elem, sp.Symbol), names)):
        return tuple(map(make_derivative_symbol, names))

    symbols_output = sp.symbols(names)
    if isinstance(symbols_output, sp.Symbol):
        return make_derivative_symbol(symbols_output)
    else:
        return tuple(map(make_derivative_symbol, sp.symbols(names)))


def make_derivative_symbol(symbol) -> sp.Symbol:
    """
    Makes symbol with the dot from input symbol.

    If symbols with index must be in form 's_{index}'.

    Example:
        .. math:: \dot{x} = make\_derivative\_symbol(x)

    """
    str_symbol = str(symbol)
    if '_' in str_symbol:
        name, index, *other = str(symbol).split('_')
        return sp.Symbol(rf'\dot {name}' + '_' + '%s' % index)
    else:
        return sp.Symbol(rf'\dot {str_symbol}')


def monom2str(monom: tuple, gens):
    return sp.latex(monomial_to_poly(sp.Monomial(monom, gens)).as_expr())


def reset_progress_bar(pbar: tqdm, value):
    pbar.n = pbar.last_print_n = value
    pbar.start_t = pbar.last_print_t = time()
    pbar.refresh()


def refresh_and_close_progress_bars(*pbars: tqdm):
    for bar in pbars:
        bar.refresh()
        bar.close()


def get_decompositions(monomial):
    if len(monomial) == 0:
        return {(tuple(), tuple())}
    result = set()
    prev_result = get_decompositions(tuple(monomial[:-1]))
    for r in prev_result:
        for i in range(monomial[-1] + 1):
            a, b = tuple(list(r[0]) + [i]), tuple(list(r[1]) + [monomial[-1] - i])
            result.add((min(a, b), max(a, b)))
    return result


def symbol_from_derivative(derivative: sp.Symbol) -> sp.Symbol:
    return sp.Symbol(str(derivative).replace(r"\dot ", '', 1))


def poly_to_monomial(poly: sp.Poly) -> sp.Monomial:
    return sp.Monomial(poly.monoms()[0], poly.gens)


def monomial_to_poly(monom: sp.Monomial) -> sp.Poly:
    return sp.Poly(sp.prod([gen ** e for gen, e in zip(monom.gens, monom.exponents)]), monom.gens)


def mlist_to_poly(mlist: Collection[sp.Monomial], gens) -> sp.Poly:
    return sp.Poly(reduce(add, mlist), gens)


def get_hull_vertices(hull: ConvexHull):
    return list(filter(lambda x: sum(x) > 1e-8, hull.points[hull.vertices].tolist()))


def dominated(monom, monom_set):
    """
    Returns true iff the monom is coordinate-wise <=
    than one of the monoms in the set
    """
    for m in monom_set:
        if all([monom[i] <= m[i] for i in range(len(m))]):
            return True
    return False


def top_priority():
    """ Set the priority of the process to above-normal."""
    import os
    if os.name == 'posix':
        os.nice(19)
    else:
        import win32api, win32process, win32con

        pid = win32api.GetCurrentProcessId()
        handle = win32api.OpenProcess(win32con.PROCESS_ALL_ACCESS, True, pid)
        win32process.SetPriorityClass(handle, win32process.REALTIME_PRIORITY_CLASS)


def apply_quadratization(polynomials: List[PolyElement], quadratization: List[Tuple], new_var='w'):
    result = list(polynomials)
    for monom in quadratization:
        result.append(calc_Lie_derivative(polynomials, monom))
    subs = generalized_variables_dict(polynomials[0].ring.symbols, quadratization, new_var)
    result = list(map(PolyElement.as_expr, result))
    for i, poly in enumerate(result):
        result[i] = poly.subs(subs, simultaneous=True)
    return result


def generalized_variables_dict(orig_vars: List[sp.Symbol], quadratization: List[Tuple], new_vars_name):
    orig_var_list = list(zip(orig_vars, orig_vars))

    new_vars = sp.symbols([f'{new_vars_name}_{i + 1}' for i in range(len(quadratization))])
    quad_syms = [monom2PolyElem(m, orig_vars) for m in quadratization]
    quad_var_list = list(zip(quad_syms, new_vars))

    res = dict(orig_var_list + quad_var_list)
    for (left_k, left_v), (right_k, right_v) in product(orig_var_list + quad_var_list, repeat=2):
        res[left_k * right_k] = left_v * right_v
    return res


def calc_Lie_derivative(polys: List[PolyElement], new_var: Tuple) -> PolyElement:
    res = 0 * polys[0]  # Do it in order to attach the ring to 0
    gens = polys[0].ring.gens
    for i in range(len(new_var)):
        if new_var[i] > 0:
            monom = deepcopy(new_var)
            poly = monom2PolyElem(monom, gens).diff(gens[i]) * polys[i]
            res += poly
    return res


def monom2PolyElem(monom: tuple, gens):
    return sp.prod([gen ** e for gen, e in zip(gens, monom)])


def apply_quad_monom(monom, quad: List):
    pass
