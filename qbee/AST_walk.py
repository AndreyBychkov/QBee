import sympy as sp
from typing import Union


def find_non_polynomial(expr: sp.Expr, mode="backward") -> Union[sp.Expr, None]:
    """
    Finds non-polynomial subexpression in 'expr'.

    :param mode: if 'forward' search goes from trunk to leaves;
                 if 'backward' search goes from leaves to trunk.

    :returns: subexpression if found; else None
    """
    if mode == "backward":
        return _find_non_polynomial_backward(expr)
    if mode == "forward":
        return _find_non_polynomial_forward(expr)


def _find_non_polynomial_forward(expr: sp.Expr) -> Union[sp.Expr, None]:
    if not is_polynomial_function(expr) and not expr.is_Symbol and not expr.is_Number:
        return expr

    results = map(_find_non_polynomial_forward, expr.args)
    results = filter(lambda e: e is not None, results)
    results = list(results)

    if results:
        return results[0]
    return None


def _find_non_polynomial_backward(expr: sp.Expr) -> Union[sp.Expr, None]:
    if expr.args:
        results = map(_find_non_polynomial_backward, expr.args)
        results = filter(lambda e: e is not None, results)
        results = list(results)
        if results:
            return results[0]

    if not is_polynomial_function(expr) and not expr.is_Symbol and not expr.is_Number:
        return expr
    return None


def is_polynomial_function(expr: sp.Expr):
    return True if expr.is_Add or \
                   expr.is_Mul or \
                   _is_positive_numeric_power(expr) \
        else False


def _is_positive_numeric_power(expr: sp.Expr):
    if not expr.is_Pow:
        return False

    power = expr.args[1]
    return power.is_Number and power > 0
