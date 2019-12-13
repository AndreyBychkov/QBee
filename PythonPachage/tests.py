import sympy as sp

from AST_walk import find_non_polynomial
from SymbolsHolder import SymbolsHolder


def test_symbols_holder():
    symbols = sp.symbols(["a", "b", "c"])
    holder = SymbolsHolder(symbols)

    holder.create_symbol()
    holder.create_symbol()

    symbols_names = set(map(lambda x: str(x), holder.get_symbols()))
    expected_names = {"a", "b", "c", "y_0", "y_1"}
    assert not symbols_names.symmetric_difference(expected_names)


def test_AST_walk():
    x, y = sp.symbols(["x", "y"])
    sin = sp.Function("sin")
    expr = x ** 2 - y + (x + y) * sin(x)

    found = find_non_polynomial(expr)
    assert found == sin(x)
