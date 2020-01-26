import sympy as sp

from SymbolsHolder import SymbolsHolder


def test_symbols_creation():
    symbols = sp.symbols(["a", "b", "c"])
    holder = SymbolsHolder(symbols)

    holder.create_symbol()
    holder.create_symbol()

    symbols_names = set(map(lambda x: str(x), holder.get_symbols()))
    expected_names = {"a", "b", "c", "y_{0}", "y_{1}"}
    assert not symbols_names.symmetric_difference(expected_names)


def test_long_indexes():
    true_res = sp.Symbol('y_{0}')

    holder = SymbolsHolder([])
    created_res = holder.create_symbol()

    assert true_res == created_res


def test_derivative_creation():
    possible_dot_expression = sp.Symbol(r'\dot y_{0}'), sp.Symbol(r'\dot{y}_{0}')

    holder = SymbolsHolder([])
    _, dot_created = holder.create_symbol_with_derivative()

    assert dot_created in possible_dot_expression
