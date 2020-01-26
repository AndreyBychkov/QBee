import sympy as sp
from collections.abc import Iterable
from typing import Tuple, List


def make_derivative_symbol(symbol) -> sp.Symbol:
    """
    Makes symbol with the dot from input symbol.

    Example:
        .. math:: \dot{x} = make\_derivative\_symbol(x)

    """
    return sp.Symbol(rf"\dot {symbol}")


class SymbolsHolder:
    """
    Class that manages symbol storage.
    """

    def __init__(self, symbols: Iterable):
        self._original_symbols = list(symbols)
        self._added_symbols = list()
        self._base_name = "y_"

    def create_symbol(self) -> sp.Symbol:
        """
        Creates a new symbol and stores it within itself.

        Example:
            .. math:: y_1 = holder.create\_symbol()

        :return: Created symbol
        """
        if not self._added_symbols:
            new_symbol = sp.Symbol(self._base_name + "0")
        else:
            last = str(self.get_last_added())
            new_index = int(last.split("_")[1]) + 1
            new_symbol = sp.Symbol(last[:-1] + str(new_index))

        self._added_symbols.append(new_symbol)
        return new_symbol

    def create_symbol_with_derivative(self) -> Tuple[sp.Symbol, sp.Symbol]:
        """
        Creates new a symbol with its derivative and stores them within itself.

        Example:
            .. math:: y_1, \dot{y}_1 = holder.create\_symbol\_with\_derivative()

        :returns: Created symbol with derivative
        """
        new_symbol = self.create_symbol()
        new_symbol_der = make_derivative_symbol(new_symbol)
        return new_symbol, new_symbol_der

    def add_symbols(self, symbols: Iterable):
        for new_symbol in symbols:
            self._added_symbols.append(new_symbol)

    def get_last_added(self) -> sp.Symbol:
        return self._added_symbols[-1]

    def get_symbols(self) -> List[sp.Symbol]:
        return self._original_symbols + self._added_symbols
