import sympy as sp
from typing import Union, Set


def find_non_polynomial(expr: sp.Expr, in_vars: Set[sp.Symbol], mode="backward") -> Union[sp.Expr, None]:
    """
    Finds non-polynomial subexpression in 'expr'.

    :param in_vars: variables that can produce nonlinearities
    :param mode: if 'forward' search goes from trunk to leaves;
                 if 'backward' search goes from leaves to trunk.

    :returns: subexpression if found; else None
    """
    if mode == "backward":
        return _find_non_polynomial_backward(expr, in_vars)
    if mode == "forward":
        return _find_non_polynomial_forward(expr, in_vars)


def _find_non_polynomial_forward(expr: sp.Expr, in_vars: Set[sp.Symbol]) -> Union[sp.Expr, None]:
    if not is_polynomial_function(expr, in_vars) and not expr.is_Symbol and not expr.is_Number:
        return expr

    results = map(lambda e: _find_non_polynomial_forward(e, in_vars), expr.args)
    results = filter(lambda e: e is not None, results)
    results = list(results)

    if results:
        return results[0]
    return None


def _find_non_polynomial_backward(expr: sp.Expr, in_vars: Set[sp.Symbol]) -> Union[sp.Expr, None]:
    if expr.args:
        results = map(lambda e: _find_non_polynomial_backward(e, in_vars), expr.args)
        results = filter(lambda e: e is not None, results)
        results = list(results)
        if results:
            return results[0]

    if not is_polynomial_function(expr) and not expr.is_Symbol and not expr.is_Number and \
            _has_common_symbol(expr.free_symbols, in_vars):
        return expr
    return None


def is_polynomial_function(expr: sp.Expr):
    return True if expr.is_Add or \
                   expr.is_Mul or \
                   _is_positive_integer_power(expr) \
        else False


def _is_positive_integer_power(expr: sp.Expr):
    return expr.is_Pow and expr.exp > 0 and expr.exp.is_integer


def _has_common_symbol(lhs: Set[sp.Symbol], *rhs: Set[sp.Symbol]):
    return any([len(lhs.intersection(r)) > 0 for r in rhs])
