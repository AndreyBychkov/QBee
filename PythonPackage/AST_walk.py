import sympy as sp


def is_polynomial_function(expr: sp.Expr):
    return True if expr.is_Add or \
                   expr.is_Mul or \
                   (expr.is_Pow and expr.args[1].is_Number) \
        else False


def find_non_polynomial(expr: sp.Expr):
    if not is_polynomial_function(expr) and not expr.is_Symbol and not expr.is_Number:
        return expr

    results = map(find_non_polynomial, expr.args)
    results = filter(lambda e: e is not None, results)
    results = list(results)

    if results:
        return results[0]
    return None
