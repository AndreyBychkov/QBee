import sympy as sp


class QBeePrinter(sp.printing.StrPrinter):
    def _print_Derivative(self, expr):
        function, *vars = expr.args
        return "{}_{}".format(
            self._print(sp.Symbol(function.func.__name__)),
            ''.join([self._print(i[0]) * i[1] for i in vars]))

    def _print_Function(self, expr: sp.Function):
        return expr.func.__name__

    def _print_Relational(self, expr):
        if isinstance(expr, sp.Eq):
            return f"{expr.lhs} = {expr.rhs}"
        return super()._print_Relational(expr)


def print_qbee(expr):
    print(str_qbee(expr))


def str_qbee(expr):
    return QBeePrinter().doprint(expr)
