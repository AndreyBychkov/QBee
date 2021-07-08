from qbee import *
from qbee.examples import *


def my_main():
    c = sp.symbols([f"c_{i}" for i in range(1, 5)])
    A, T, Ru, Ea = sp.symbols("A, T, R_u, E_a")
    dc = derivatives(c)
    temp = -A * sp.exp(-Ea / (Ru * T)) * c[0] ** 0.2 * c[1] ** 1.3
    eq = EquationSystem([
        sp.Eq(dc[0], temp),
        sp.Eq(dc[1], 2 * temp),
        sp.Eq(dc[2], -temp),
        sp.Eq(dc[3], -2 * temp)],
        parameter_variables=[A, Ru, Ea],
        input_variables=[T])
    peq = polynomialize(eq)
    peq.print()
    print("-" * 50)
    print(peq.substitution_equations)


if __name__ == '__main__':
    my_main()
