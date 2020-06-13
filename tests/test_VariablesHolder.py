import pytest
import sympy as sp

from qbee.structures import VariablesHolder


def test_symbols_creation():
    symbols = sp.symbols(["a", "b", "c"])
    holder = VariablesHolder(symbols)

    holder.create_variable()
    holder.create_variable()

    symbols_names = set(map(lambda x: str(x), holder.free))
    expected_names = {"a", "b", "c", "y_{0}", "y_{1}"}
    assert not symbols_names.symmetric_difference(expected_names)


def test_long_indexes():
    true_res = sp.Symbol('y_{0}')

    holder = VariablesHolder([])
    created_res = holder.create_variable()

    assert true_res == created_res


def test_derivative_creation():
    possible_dot_expression = sp.Symbol(r'\dot y_{0}'), sp.Symbol(r'\dot{y}_{0}')

    holder = VariablesHolder([])
    _, dot_created = holder.create_variable_with_derivative()

    assert dot_created in possible_dot_expression
