import pickle
import pandas as pd
import sympy as sp
import configparser
from sympy.polys.rings import PolyElement
from tqdm.autonotebook import tqdm
from typing import Iterable, Union, Tuple, List, Dict, Optional
from functools import reduce, wraps
from .printer import str_qbee

INDEPENDENT_VARIABLE = sp.Symbol("_t")


def functions(names, *, laurent=True, **kwargs):
    """Dependent and input variables in differential equations. The syntax is identical to `sympy.symbols`_.

    Examples
        >>> x, y, z = functions('x, y, z')

    .. _sympy.symbols: https://docs.sympy.org/latest/modules/core.html?highlight=symbols#sympy.core.symbol.symbols
    """
    t = INDEPENDENT_VARIABLE
    funcs = sp.symbols(names, cls=sp.Function, is_laurent=laurent, **kwargs)
    if isinstance(funcs, Iterable):
        return [f(t) for f in funcs]
    return funcs(t)


def multivariable_functions(names, args, *, laurent=True, **kwargs):
    funcs = sp.symbols(names, cls=sp.Function, is_laurent=laurent, **kwargs)
    if isinstance(funcs, Iterable):
        return [f(*args) for f in funcs]
    return funcs(*args)


class Parameter(sp.Symbol):
    pass


def parameters(names, **kwargs):
    """Parameter variables in differential equations. The syntax is identical to `sympy.symbols`_.

        Examples
            >>> p, k = parameters('p, k')

        .. _sympy.symbols: https://docs.sympy.org/latest/modules/core.html?highlight=symbols#sympy.core.symbol.symbols
        """
    return sp.symbols(names, cls=Parameter, **kwargs)


class Config:
    def __init__(self, progress_bar_enabled=True,
                 logging_enabled=False,
                 logging_file="log.feather",
                 quad_systems_file="quad_systems.pkl"):
        self.progress_bar_enabled = progress_bar_enabled
        self.logging_enabled = logging_enabled
        self.logging_file = logging_file
        self.quad_systems_file = quad_systems_file

    def update_from_ini_file(self, path, section="DEFAULT"):
        parser = configparser.ConfigParser()
        with open(path, 'r') as f:
            parser.read_string(f.read())  # Because ConfigParser.read_file somehow does not work
            if pb_enabled := parser.get(section, "progress_bar_enabled"):
                self.progress_bar_enabled = eval(pb_enabled)
            if log_enabled := parser.get(section, "logging_enabled"):
                self.logging_enabled = eval(log_enabled)
            if log_file := parser.get(section, "logging_file"):
                self.logging_file = log_file
            if quad_file := parser.get(section, "quad_systems_file"):
                self.quad_systems_file = quad_file


QBEE_CONFIG = Config()


def parametrized(dec):
    def layer(*args, **kwargs):
        def repl(f):
            return dec(f, *args, **kwargs)

        return repl

    return layer


__log_pb_evaluated = False
__log_pb = None
__log_records = list()


@parametrized
def progress_bar(func, is_stop):
    postfix_str = "Current best order = {}"

    @wraps(func)
    def wrapper(*args, **kwargs):
        if QBEE_CONFIG.progress_bar_enabled:
            global __log_pb_evaluated
            global __log_pb
            if not __log_pb_evaluated:
                __log_pb_evaluated = True
                __log_pb = tqdm(desc="Nodes processed", unit=f" nodes", leave=False)
                __log_pb.set_postfix_str(postfix_str.format(args[-1]))
            if is_stop:
                __log_pb.close()
                __log_pb = None
                __log_pb_evaluated = False
            else:
                __log_pb.update(1)
                __log_pb.set_postfix_str(postfix_str.format(args[-1]))
                __log_pb.refresh()
        return func(*args, **kwargs)

    return wrapper


@parametrized
def logged(method, is_stop=False):
    @wraps(method)
    def wrapper(self, *args, **kwargs):
        res = method(self, *args, **kwargs)
        if QBEE_CONFIG.logging_enabled:
            global __log_records
            if is_stop:
                df = pd.DataFrame(__log_records, columns=['from', 'to']).applymap(str)
                df.to_feather(QBEE_CONFIG.logging_file)
                __log_records = list()
            else:
                part_res, *_ = args
                for r in list(res):
                    __log_records.append([part_res, r])
        return res

    return wrapper


def dump_results(method):
    @wraps(method)
    def wrapper(self, *args, **kwargs):
        res = method(self, *args, **kwargs)
        if QBEE_CONFIG.logging_enabled:
            pickle.dump(res, open(QBEE_CONFIG.quad_systems_file, "wb"))
        return res

    return wrapper


@parametrized
def memoize_first(func, max_size):
    memo = {}

    @wraps(func)
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
        for i in range(abs(monomial[-1]) + 1):
            i = i if monomial[-1] >= 0 else -i
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
    return str_qbee(sp.Monomial(monom, gens).as_expr())


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


def Lie_derivative_on_dicts(monom, vector_field):
    result = dict()
    for i in range(len(monom)):
        if monom[i] != 0:
            for m, c in vector_field[i].items():
                new_monom = [k + l for k, l in zip(monom, m)]
                new_monom[i] -= 1
                new_c = c * monom[i]
                new_monom = tuple(new_monom)
                result[new_monom] = result.get(new_monom, 0) + new_c
    return result


def apply_quadratization(polynomials: List[PolyElement], quadratization: List[Tuple],
                         new_var_name='z_', start_new_vars_with=0):
    n = len(polynomials)
    system_dicts = [p.to_dict() for p in polynomials]
    vector_field = {i: p for i, p in enumerate(system_dicts)}
    for monom in quadratization:
        system_dicts.append(Lie_derivative_on_dicts(monom, vector_field))

    squares = dict()
    generalized_vars = [tuple([1 if i == j else 0 for i in range(n)]) for j in range(n)] + list(quadratization)
    squares[tuple([0] * n)] = tuple([0] * len(generalized_vars))
    for i, v in enumerate(generalized_vars):
        squares[v] = tuple([1 if j == i else 0 for j in range(len(generalized_vars))])

    for i, v in enumerate(generalized_vars):
        for j, u in enumerate(generalized_vars):
            prod = tuple([a + b for a, b in zip(u, v)])
            new_monom = [0] * len(generalized_vars)
            new_monom[i] += 1
            new_monom[j] += 1
            new_monom = tuple(new_monom)
            if prod not in squares:
                squares[prod] = new_monom
            else:
                # prioritizing new variables over the old ones
                squares[prod] = min(squares[prod], new_monom)

    result_as_dicts = [{squares[m]: c for m, c in p.items()} for p in system_dicts]
    quad_varnames = [new_var_name + '{' + str(start_new_vars_with + i) + '}' for i in range(len(quadratization))]
    new_varnames = [str(g) for g in polynomials[0].ring.gens] + quad_varnames
    new_ring = sp.ring(new_varnames, polynomials[0].ring.domain)[0]
    result = [new_ring(d) for d in result_as_dicts]
    return result, [g.as_expr() for g in new_ring.gens], [g.as_expr() for g in new_ring.gens if str(g) in quad_varnames]


def monom2PolyElem(monom: tuple, gens):
    return sp.prod([gen ** e for gen, e in zip(gens, monom)])


def generate_derivatives(inputs: Dict[PolyElement, int]) -> List[List[str]]:
    return [
        [sp.Symbol(str(var) + '\'' * n) for n in range(0, max_order + 1)]
        for var, max_order in inputs.items()
    ]


def remove_repeating(seq: Iterable) -> List:
    return list(dict.fromkeys(seq))


def key_from_value(d: dict, value) -> Optional:
    for k, v in d.items():
        if v == value:
            return k
    return None
