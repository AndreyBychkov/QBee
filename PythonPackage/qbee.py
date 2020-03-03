import sympy as sp
from typing import List, Iterable, Optional
from EquationSystem import EquationSystem


def polynomialize(system: List[sp.Eq],
                  parameter_variables: Iterable[sp.Symbol] = None,
                  input_variables: Iterable[sp.Symbol] = None,
                  mode='differential') -> List[sp.Eq]:
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


def quadratic_linearize(system: List[sp.Eq],
                        parameter_variables: Iterable[sp.Symbol] = None,
                        input_variables: Iterable[sp.Symbol] = None,
                        mode='differential',
                        method='sqrt-count-first',
                        debug: Optional[str] = None,
                        log_file: Optional[str] = None) -> List[sp.Eq]:
    """
    Transforms the system into quadratic-linear form using variable replacement technique.

    :param system: system of differential equations in form x' = f(x)
    :param parameter_variables: constant parameter variables
    :param input_variables: variables representing input functions
    :param mode: auxiliary equation form.
    :param method: next replacement choice method.
    :param debug: printing mode while quadratic linearization is performed.
    :param log_file: output file for evaluation logging. Must be in 'csv' format.


    Mode
    -----------------
    **algebraic**
        adds auxiliary equations in form y = f(x, y)
    **differential**
         adds auxiliary equations in form y' = f(x, y)

    Method
    -----------------
    **random**
        choose next possible replacement in random way;
    **count-first**
        choose most frequent possible replacement as the next one;
    **sqrt-first**
        choose next possible replacement within variables' squares in random way;
    **sqrt-count-first**
         choose most frequent square replacement as the next one;

    Debug
    ---------------
    **None** or **silent**
        prints nothing;
    **info**
        prints replacement for each iteration;
    **debug**
        prints equations in system with replacement for each iteration;

    """
    eq_system = EquationSystem(system, parameter_variables, input_variables)
    if not eq_system.is_polynomial('full'):
        eq_system.polynomialize(mode)
    eq_system.quadratic_linearized(auxiliary_eq_type=mode, heuristics=method, debug=debug, log_file=log_file)
    return eq_system.equations
