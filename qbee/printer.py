import sympy as sp


class QBeePrinter(sp.printing.StrPrinter):
    def _print_Derivative(self, expr):
        function, *vars = expr.args
        return "{}_{}".format(
            self._print(sp.Symbol(function.func.__name__)),
            ''.join([self._print(i[0]) * i[1] for i in vars]))

    def _print_Function(self, expr: sp.Function):
        # TODO
        # from qbee import INDEPENDENT_VARIABLE
        if sp.Symbol("_t") in expr.args:
            return expr.func.__name__
        return super()._print_Function(expr)

    def _print_Relational(self, expr):
        if isinstance(expr, sp.Eq):
            return f"{self._print(expr.lhs)} = {self._print(expr.rhs)}"
        return super()._print_Relational(expr)

    def _print_Symbol(self, expr):
        return str(expr).replace('{', '').replace('}', '')


def print_qbee(expr):
    print(str_qbee(expr))


def str_qbee(expr):
    return QBeePrinter().doprint(expr)
