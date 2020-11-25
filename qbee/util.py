import pandas as pd
import sympy as sp
import numpy as np
from time import time
from tqdm import tqdm
from typing import Iterable, Collection, Union, Tuple
from sympy.polys.monomials import monomial_deg
from itertools import combinations
from functools import reduce, wraps
from operator import add


def parametrized(dec):
    def layer(*args, **kwargs):
        def repl(f):
            return dec(f, *args, **kwargs)

        return repl

    return layer


def timed(func):
    def wrapper(*args, **kwargs):
        start_time = time()
        res = func(*args, **kwargs)
        end_time = time()
        print()
        print(f"Elapsed time: {np.round(end_time - start_time, 3)}s.")
        return res

    return wrapper


__log_pb_evaluated = False
__log_pb = None
__log_records = list()


@parametrized
def progress_bar(func, is_stop):
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

    return wrapper


@parametrized
def logged(method, enabled, log_file, is_stop=False):
    def wrapper(self, *args, **kwargs):
        global __log_records
        res = method(self, *args, **kwargs)
        if is_stop:
            pd.DataFrame(__log_records, columns=['from', 'to']).to_csv(log_file, index=None)
            __log_records = list()
        else:
            part_res, *_ = args
            __log_records.append([part_res, list(res)])
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


def can_substitutions_quadratize(monom: sp.Monomial, subs: Iterable[sp.Monomial]) -> bool:
    gens = tuple(monom.gens)
    var_subs = set(map(lambda var: sp.Monomial(var, gens), gens))
    var_subs.add(sp.Monomial((0,) * len(gens), gens))  # zero degree monomial
    expanded_subs = var_subs.union(subs)
    return _can_quad_0(monom) or _can_quad_2(monom, expanded_subs) or _can_quad_1(monom, expanded_subs)


def _can_quad_2(monom: sp.Monomial, subs: Iterable[sp.Monomial]) -> bool:
    for sub in subs:
        if monom == sub ** 2:
            return True
    return False


def _can_quad_1(monom: sp.Monomial, subs: Iterable[sp.Monomial]) -> bool:
    for left, right in combinations(subs, 2):
        if monom == left * right:
            return True
    return False


def _can_quad_0(monom: sp.Monomial) -> bool:
    return monomial_deg(monom) == 0
