import sympy as sp
from typing import List, Iterable
from EquationSystem import EquationSystem


def polynomialize(system: List[sp.Eq],
                  parameter_variables: Iterable[sp.Symbol] = None,
                  input_variables: Iterable[sp.Symbol] = None,
                  mode='algebraic') -> List[sp.Eq]:
    """
        Transforms the system into polynomial form using variable replacement techniques.

        :param system: system of differential equations in form x' = f(x)
        :param parameter_variables: constant parameter variables
        :param input_variables: variables representing input functions
        :param mode: auxiliary equation form.

        Mode
        -----------------
        **algebraic**
            adds auxiliary equations in form y = f(x, y)
        **differential**
             adds auxiliary equations in form y' = f(x, y)

        """
    eq_system = EquationSystem(system, parameter_variables, input_variables)
    eq_system.polynomialize(mode)
    return eq_system.equations
