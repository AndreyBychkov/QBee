Polynomialization
===================

Polynomialization by example
------------------------------

Polynomialization is extremely close in spirit to :ref:`quadratization <Quadratization>`.
The difference is that this time we are converting systems with the right-hand side as elementary functions
to systems with polynomial right-hand sides.
Let's look at an example:
$$
x' = x \sin(a x) + u(t)
$$

.. code-block:: python

    >>> import sympy as sp
    >>> from qbee import *
    >>> x, u = functions("x, u")
    >>> a = parameters("a")
    >>> system = [(x, x * sp.sin(a*x) + u)]  # u is an input variable
    >>> polynomialize(system).print()
    Introduced variables:
    w_0 = sin(a*x)
    w_1 = cos(a*x)

    x' = u + w_0*x
    w_0' = a*w_1*(u + w_0*x)
    w_1' = -a*w_0*(u + w_0*x)

Note that the right-hand side of the resulting system is polynomial.
We will call the set of introduced variables $\{\cos(a \cdot x), \sin(a \cdot x)\}$ a **polynomialization** of the system.

Result of polynomialization can be natively combined with quadratization:

.. code-block:: python

    >>> poly_res = polynomialize(system)
    >>> quadratize(poly_res, new_vars_name="p_").print()
    Introduced variables:
    p_0 = w_1*x
    p_1 = w_0*x

    x' = p_0 + u
    w_0' = -a*p_0*w_1 - a*u*w_1
    w_1' = a*p_1*w_1 + a*u*w_0
    p_0' = a*p_0*p_1 + a*p_1*u + p_0*w_1 + u*w_1
    p_1' = -a*p_0**2 - a*p_0*u + p_1*w_1 + u*w_0

Nevertheless, there is a much simpler way to do it altogether:

.. code-block:: python

    >>> polynomialize_and_quadratize(system).print()
    Introduced variables:
    w_0 = sin(a*x)
    w_1 = cos(a*x)
    w_2 = w_0*x
    w_3 = w_1*x

    x' = u + w_2
    w_0' = a*u*w_1 + a*w_1*w_2
    w_1' = -a*u*w_0 - a*w_0*w_2
    w_2' = a*u*w_3 + a*w_2*w_3 + u*w_0 + w_0*w_2
    w_3' = -a*u*w_2 - a*w_2**2 + u*w_1 + w_1*w_2


Polynomialization and Laurent monomials
---------------------------------------------

Consider the following example:

.. code-block:: python

    >>> import sympy as sp
    >>> from qbee import *
    >>> x = functions("x")
    >>> system = [(x, sp.sin(1 / x))]
    >>> polynomialize(system).print()
    Introduced variables:
    w_0 = sin(1/x)
    w_1 = cos(1/x)

    x' = w_0
    w_0' = -w_0*w_1/x**2
    w_1' = w_0**2/x**2

As you can see, the resulting system is polynomial in :ref:`Laurent sense <Laurent monomials>`, i.e.
have variables with negative integer powers.
Although this form comes in handy for more optimal quadratization in the end
(see :ref:`this section <Optional: Polynomialization of rational powers>`),
you may want to get a system with a polynomial right-hand side in the usual sense.
How to do this is shown in the :ref:`corresponding section <Laurent Monomials>`.

Optional: Polynomialization of rational powers
----------------------------------------------

TODO, use my note

