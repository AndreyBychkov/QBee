Numerical Integration
======================

You can quickly check the quadratization numerically by using a function ``qbee.experimental.to_odeint``
that turns the result into a numerical function that scipy understands. Let's look at the in the example:

We define and quadratize a scalar system
$$
x' = p \sin(x) + u(t)
$$

.. code-block:: python

    >>> import sympy as sp
    >>> from qbee import *
    >>> x, u = functions("x, u")
    >>> p = parameters("p")
    >>> system = [(x, p * sp.sin(x) + u)]
    >>> res = polynomialize_and_quadratize(system, input_free=True)
    >>> res.print(str_function=sp.latex)
    Introduced variables:
    w_{0} = \sin{\left(x \right)}
    w_{1} = \cos{\left(x \right)}

    x' = p w_{0} + u
    w_{0}' = p w_{0} w_{1} + u w_{1}
    w_{1}' = - p w_{0}^{2} - u w_{0}

Here we took advantage of the ability to both get a quadratic system in LaTeX and render it in the documentation.

$$
w_{0} = \sin{\left(x \right)}, w_{1} = \cos{\left(x \right)}
$$
$$
\begin{array}{l}
x' = p w_{0} + u \\
w_{0}' = p w_{0} w_{1} + u w_{1} \\
w_{1}' = - p w_{0}^{2} - u w_{0}
\end{array}
$$

Now we solve the system numerically for $t \in [0; 100]$, $p = 0.1$, $u(t) = \sin(t)$:

.. code-block:: python

    >>> from qbee.experimental import to_odeint
    >>> ts = np.linspace(0, 100, 10000)
    >>> t = INDEPENDENT_VARIABLE
    >>> my_odeint = to_odeint(system, {"x": 0.1}, {"p": 0.1}, {"u": sp.sin(t)})
    >>> my_odeint_quad = to_odeint(res, {"x": 0.1}, {"p": 0.1}, {"u": sp.sin(t)})

    >>> num_res = my_odeint(ts, rtol=1e-10)
    >>> num_res_quad = my_odeint_quad(ts, rtol=1e-10)
    >>> print("L-infinity distance between numerical solutions of the original ODE and the quadratized one: ",
    >>>       np.linalg.norm(num_res[:, 0] - num_res_quad[:, 0], np.inf))
    L-infinity distance between numerical solutions of the original ODE and the quadratized one: 1.5569859268538266e-07


