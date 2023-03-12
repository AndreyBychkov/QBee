import numpy as np
import sympy as sp
from qbee import *
from qbee import INDEPENDENT_VARIABLE
from qbee.experimental import to_odeint

if __name__ == '__main__':
    x, u = functions("x, u")
    p = parameters("p")
    system = [
        (x, p * sp.sin(x) + u),
    ]
    res = polynomialize_and_quadratize(system, input_der_orders={u: 0})
    print(res)

    ts = np.linspace(0, 100, 10000)
    t = INDEPENDENT_VARIABLE
    my_odeint = to_odeint(system, {"x": 0.1}, {"p": 0.1}, {"u": sp.sin(t)})
    my_odeint_quad = to_odeint(res, {"x": 0.1}, {"p": 0.1}, {"u": sp.sin(t)})

    num_res = my_odeint(ts, rtol=1e-10)
    num_res_quad = my_odeint_quad(ts, rtol=1e-10)
    print(f"L-infinity distance between numerical solutions of the original ODE and the quadratized one: "
          f"{np.linalg.norm(num_res[:, 0] - num_res_quad[:, 0], np.inf)}")
