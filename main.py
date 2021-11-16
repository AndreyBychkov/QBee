from qbee import *
from examples import *
from copy import deepcopy
from functools import partial

if __name__ == '__main__':
    x, c = functions("x, c")
    u = multivariable_functions("u", [x])
    p = parameters("p")
    system = [
        (u, p * u.diff(x) ** 3 + c),
    ]
    system.append((u.diff(x), system[0][1].diff(x)))
    system.append((system[1][0].diff(x), system[1][1].diff(x)))
    res = polynomialize_and_quadratize(system, {u.diff(x, 3): 0, c: 0})
    if res:
        print(res)
