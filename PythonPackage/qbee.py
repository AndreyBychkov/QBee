import sympy as sp
from typing import List, Iterable
from EquationSystem import EquationSystem


def polynomialize(system: List[sp.Eq],
                  parameter_variables: Iterable[sp.Symbol] = None,
                  input_variables: Iterable[sp.Symbol] = None,
                  mode='algebraic') -> List[sp.Eq]:
    eq_system = EquationSystem(system, parameter_variables, input_variables)
    eq_system.polynomialize(mode)
    return eq_system.equations
