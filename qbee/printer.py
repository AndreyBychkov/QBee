import sympy as sp


class CommonPrinter(sp.printing.StrPrinter):
    def _print_Derivative(self, expr):
        function, *vars = expr.args
        return "{}_{}".format(
            self._print(sp.Symbol(function.func.__name__)),
            ''.join([self._print(i[0]) * i[1] for i in vars]))

    def _print_Function(self, expr: sp.Function):
        return expr.func.__name__


def print_common(expr):
    print(str_common(expr))

def str_common(expr):
    return CommonPrinter().doprint(expr)
