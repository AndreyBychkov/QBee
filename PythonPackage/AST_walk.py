import sympy as sp


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


def find_non_polynomial(expr: sp.Expr):
    if not is_polynomial_function(expr) and not expr.is_Symbol and not expr.is_Number:
        return expr

    results = map(find_non_polynomial, expr.args)
    results = filter(lambda e: e is not None, results)
    results = list(results)

    if results:
        return results[0]
    return None
