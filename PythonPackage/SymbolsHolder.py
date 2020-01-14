import sympy as sp
from collections.abc import Iterable


class SymbolsHolder:

    def __init__(self, symbols: Iterable):
        self._original_symbols = list(symbols)
        self._added_symbols = list()
        self._base_name = "y_"

    def create_symbol(self):
        if not self._added_symbols:
            new_symbol = sp.Symbol(self._base_name + "0")
        else:
            last = str(self.get_last_added())
            new_index = int(last.split("_")[1]) + 1
            new_symbol = sp.Symbol(last[:-1] + str(new_index))

        self._added_symbols.append(new_symbol)
        return new_symbol

    def create_symbol_with_derivative(self):
        new_symbol = self.create_symbol()
        new_symbol_der = sp.Symbol(rf"\dot {new_symbol}")
        return new_symbol, new_symbol_der

    def get_last_added(self):
        return self._added_symbols[-1]

    def get_symbols(self):
        return self._original_symbols + self._added_symbols
