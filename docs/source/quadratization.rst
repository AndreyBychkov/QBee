Quadratization
================

Quadratization by example
-------------------------

To better understand quadratization, let's consider the following example. We have a scalar ordinary differential equation (ODE) given by:
$$
x' = x^5
$$
Then, introduce a new variable $y = x^4$, leading to an ODE with quadratic right-hand side:
$$
x' = xy
$$
However, this new equation does not provide any additional information.
To add up, we can differentiate the equation of y.
This results in a system of ODEs that is equivalent to the original equation:
$$
\begin{array}{ll}
x' = xy \\
y' = 4x^3 x' = 4x^8 = 4 y^2
\end{array}
$$

Moreover, both equations have a quadratic right-hand side.
We refer to the set $\{y = x^4\}$ as a **quadratization**.
A formal definition of quadratization will be provided later.

But now, let's use QBee to replicate the example:

.. code-block:: python

    >>> from qbee import *
    >>> x = functions("x")
    >>> quadratize([(x, x**5)]).print()  # input follows rule (x, p) => x' = p
    Introduced variables:
    w0 = x**4

    x' = w0*x
    w0' = 4*w0**2

Formal definition
--------------------

Consider a polynomial system of ODEs
$$
\mathbf{\dot x} = \mathbf{p}(\mathbf{x}),
$$
where $\mathbf{x} = \mathbf{x}(t) = [x_1(t), \dots, x_N(t)]^T$ and $\mathbf{p}(\mathbf{x}) = [p_1(\mathbf{x}),\dots, p_N(\mathbf{x})]^T$
with $p_1,\dots, p_N \in \mathbb{C} [\mathbf{x}]$ (are polynomials in $\mathbb{C}$).

Then an $l$-dimensional vector of new variables
$$
\mathbf{w} = \mathbf{w}(\mathbf{x}) \quad \in \mathbb{C}[\mathbf{x}]^l
$$
is said to be a **quadratization** of the given system if there exist vectors $\mathbf{q}_1(\mathbf{x}, \mathbf{w})$
and $\mathbf{q}_2(\mathbf{x}, \mathbf{w})$ of dimension $N$ and $l$ respectively, with entries being polynomials
of total degree at most two (quadratic) such that
$$
\dot{\mathbf{x}} = \mathbf{q}_1(\mathbf{x}, \mathbf{w}) \quad and \quad \dot{\mathbf{w}} = \mathbf{q}_2(\mathbf{x}, \mathbf{w}).
$$

The dimension $l$ of vector $\mathbf{w}$ is called **order of quadratization**.
A quadratization of the smallest possible order is called an **optimal quadratization**.

In the example above $\mathbf{w} = [x^4],\ l = 1,\ \text{and}\ \mathbf{q}_1 = [w_0 x],\ \mathbf{q}_2 = [4w_0^2]$.
Since we introduced only one variable, the resulting quadratization is optimal.


Monomial quadratization
------------------------

If all the polynomials $w_1(\mathbf{x}),\dots,w_l(\mathbf{x})$ are
monomials, the quadratization is called a **monomial quadratization**.
Accordingly, if a monomial quadratization of a system has the smallest possible order among all
the monomial quadratizations of the system, it is called an **optimal monomial quadratization**.

The quadratization in the first example is indeed an optimal monomial quadratization.

The importance of monomial quadratizations is that optimal monomial quadratizations
always exist when there is no :ref:`input functions <Inputs and Parameters>` involved
and $w_1(\mathbf{x}),\dots,w_l(\mathbf{x})$ are non-negative monomials
(:ref:`Otherwise <Laurent monomials>`, quadratization may have a smaller order
but we are not aware of algorithm to find an optimal quadratization in this case).


Laurent monomials
----------------------





