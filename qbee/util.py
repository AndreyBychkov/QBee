import pickle
import pandas as pd
import sympy as sp
import numpy as np
from sympy.polys.rings import PolyElement
from time import time
from tqdm import tqdm
from typing import Iterable, Union, Tuple, List, Dict
from itertools import product, chain, combinations
from functools import reduce
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
        return sp.Symbol(rf"{name}_{index}'")
    else:
        return sp.Symbol(rf"{symbol}'")


def monom2str(monom: tuple, gens):
    return sp.latex(monomial_to_poly(sp.Monomial(monom, gens)).as_expr())


def symbol_from_derivative(derivative: sp.Symbol) -> sp.Symbol:
    return sp.Symbol(str(derivative).replace("\'", '', 1))


def monomial_to_poly(monom: sp.Monomial) -> sp.Poly:
    return sp.Poly(sp.prod([gen ** e for gen, e in zip(monom.gens, monom.exponents)]), monom.gens)


def tuple_to_monom(m: tuple, gens) -> sp.Expr:
    return sp.prod([g ** e for g, e in zip(gens, m)])


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


def apply_quadratization(polynomials: List[PolyElement], quadratization: List[Tuple], params, new_var_name='z_'):
    gens = [g for g in polynomials[0].ring.gens if str(g) not in map(str, params)]
    result = list(polynomials)
    for monom in quadratization:
        result.append(calc_Lie_derivative(polynomials, monom, gens))
    subs = generalized_variables_dict(gens, list(map(PolyElement.as_expr, quadratization)), new_var_name)
    result = list(map(PolyElement.as_expr, result))
    for i, poly in enumerate(result):
        ppoly = sp.Poly(poly, *[m.as_expr() for m in gens if m.as_expr() not in params])
        res_lst = list()
        for monom, coef in ppoly.terms():
            subs_monom = monom2PolyElem(monom, ppoly.gens).as_expr().xreplace(subs)
            res_lst.append(coef * subs_monom)
        result[i] = sp.Add(*res_lst)
    return result


def generalized_variables_dict(orig_vars: List[sp.Symbol], quadratization: List[PolyElement], new_vars_name):
    orig_vars = list(map(lambda m: m.as_expr(), orig_vars))
    orig_var_list = list(zip(orig_vars, orig_vars))

    new_vars = sp.symbols([new_vars_name + "{%d}" % i for i in range(len(quadratization))])
    quad_var_list = list(zip(map(lambda m: m.as_expr(), quadratization), new_vars))

    res = dict(orig_var_list + quad_var_list)
    for (left_k, left_v), (right_k, right_v) in product(orig_var_list + quad_var_list, repeat=2):
        key = left_k * right_k
        val = left_v * right_v
        if key != val:
            res[key] = val
    return res


def calc_Lie_derivative(polys: List[PolyElement], new_var: PolyElement, gens) -> PolyElement:
    res = polys[0].ring(0)
    new_var_tup = new_var.monoms()[0]
    for i in range(len(new_var_tup)):
        if new_var_tup[i] > 0:
            dx = new_var.ring(sp.Symbol(str(gens[i])))
            poly = new_var.diff(dx).set_ring(polys[0].ring) * polys[i]
            res += poly
    return res


def monom2PolyElem(monom: tuple, gens):
    return sp.prod([gen ** e for gen, e in zip(gens, monom)])


def generate_derivatives(inputs: Dict[PolyElement, int]) -> List[List[str]]:
    return [
        [sp.Symbol(str(var) + '\'' * n) for n in range(0, max_order + 1)]
        for var, max_order in inputs.items()
    ]
