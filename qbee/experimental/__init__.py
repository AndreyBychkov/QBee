import sympy as sp


class Variable(sp.Symbol):
    pass


class Input(sp.Symbol):
    pass


class Parameter(sp.Symbol):
    pass


def variables(names, **kwargs):
    return sp.symbols(names, cls=Variable, **kwargs)


def inputs(names, **kwargs):
    return sp.symbols(names, cls=Input, **kwargs)


def parameters(names, **kwargs):
    return sp.symbols(names, cls=Parameter, **kwargs)
