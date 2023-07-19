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

The polynomial quadratization, on the other hand, is generally of a lower order, as seen in the example below.
However, the general algorithm of finding polynomial quadratizations is unknown to us.
Nevertheless, shifting variables can help you to produce better order when using QBee.

**Example**: The system
$$
x' = (x + 1)^{10}
$$
has a polynomial quadratization of order $1$:
$$
\mathbf{w} = [(x+1)^9],
$$
while its optimal monomial quadratization is
$$
\mathbf{v} = [x^2, x^4, x^6, x^8, x^9].
$$
Replacing $x$ with $y = x+1$ will have a result similar to the polynomial case, yet reachable by using QBee:
$$
y' = y^{10}, \quad \mathbf{v} = [x^9].
$$

As a further read, we recommend an article by Foyez Allaudin "`Quadratization of ODEs: Monomial vs. Non-Monomial <https://arxiv.org/abs/2011.03959>`_".


Laurent monomials
----------------------

Laurent monomials are monomials that may be built with both positive and negative powers, for instance $m(x, y) = x^2 y^{-3}$.
By default, QBee introduces new variables as Laurent monomials, since it results in more optimal quadratizations than by using ordinary monomials.

**Example**: Let's quadratize the system below by using QBee:
$$
x' = \exp(1/x).
$$

.. code-block:: python

    >>> import sympy as sp
    >>> from qbee import *
    >>> x = functions("x")
    >>> polynomialize_and_quadratize([(x, sp.exp(1/x))]).print()
    Introduced variables:
    w_0 = exp(1/x)
    w_1 = w_0/x**2  <- It's a Laurent monomial;
    w_2 = w_0/x**3  <- This one as well;

    x' = w_0
    w_0' = -w_0*w_1
    w_1' = -2*w_0*w_2 - w_1**2
    w_2' = -3*w_1**2 - w_1*w_2

**Note**: We first polynomialize the system to bring it to a polynomial form. Read about this technique in :ref:`Polynomialization section <Polynomialization>`.

However, in some contexts, using Laurent monomials is undesirable or even forbidden.
For such cases, we can disable them:

.. code-block:: python

    >>> import sympy as sp
    >>> from qbee import *
    >>> x = functions("x", laurent=False)  # Forbid "x" to be of negative powers in monomials
    >>> polynomialize_and_quadratize([(x, sp.exp(1/x))]).print()
    Introduced variables:  # Note that the quadratization order has grown;
    w_0 = exp(1/x)
    w_1 = 1/x  <- We introduce 1/x to tackle with the lack of Laurent monomials;
    w_2 = w_0*w_1**2
    w_3 = w_0*w_1**3

    x' = w_0
    w_0' = -w_0*w_2
    w_1' = -w_2
    w_2' = -2*w_0*w_3 - w_2**2
    w_3' = -3*w_2**2 - w_2*w_3

For more details and context, check Remark 5.1 in
`Exact and optimal quadratization of nonlinear finite-dimensional non-autonomous dynamical systems <https://doi.org/10.48550/arXiv.2303.10285>`_.