[![Documentation Status](https://readthedocs.org/projects/qbee/badge/?version=latest)](https://qbee.readthedocs.io/en/dev/?badge=latest)
# QBee

### [Online playground](https://huggingface.co/spaces/Armaliltril/qbee)
### [Tutorial in Jupyter](QBee_tutorial.ipynb) ([Colab](https://colab.research.google.com/github/AndreyBychkov/QBee/blob/dev/QBee_tutorial.ipynb))

QBee is a Python library for transforming systems of differential equations into a systems with quadratic right-rand side.

# Installation

## PyPI
`pip install qbee`

## Manual
1. Clone repository: `https://github.com/AndreyBychkov/QBee.git`
   * Or, if you want our bleeding edge version, clone `https://github.com/AndreyBychkov/QBee/tree/dev`
2. Change directory: `cd QBee`
3. Install package: `pip install .`

If you use `poetry` you can alternately install if with
`poetry install`

# What is quadratization?

The problem of *quadratization* is, given a system of ODEs with polynomial right-hand side, reduce the system to a
system with quadratic right-hand side by introducing as few new variables as possible. We will explain it using toy
example. Consider the system

<img alt="\begin{cases} x_1&#39; = x_1 x_2 \\ x_2&#39; = -x_1 x_2^3 \end{cases}" height="135" src="https://latex.codecogs.com/png.latex?\dpi{200}&amp;space;\huge&amp;space;{\color{DarkOrange}&amp;space;\begin{cases}&amp;space;x_1&#39;&amp;space;=&amp;space;x_1&amp;space;x_2&amp;space;\\&amp;space;x_2&#39;&amp;space;=&amp;space;-x_1&amp;space;x_2^3&amp;space;\end{cases}}" width="200"/>

An example of quadratization of this system will be a new variable

<img alt="y = x_1 x_2^2" height="60" src="https://latex.codecogs.com/png.latex?\dpi{200}&amp;amp;amp;space;\huge&amp;amp;amp;space;{\color{DarkOrange}y&amp;amp;amp;space;=&amp;amp;amp;space;x_1&amp;amp;amp;space;x_2^2}" width="150"/>

leading to the following ODE

<img alt="y&#39; = x_2 y - 2y^2" height="50" src="https://latex.codecogs.com/png.latex?\dpi{200}&amp;space;\huge&amp;space;{\color{DarkOrange}y&#39;&amp;space;=&amp;space;x_2&amp;space;y&amp;space;-&amp;space;2y^2}" width="250"/>

Thus, we attained the system with quadratic right-hand side

<img alt="\begin{cases} x_1&#39; = x_1 x_2 \\ x_2&#39; = -x_2 y \\ y&#39; = x_2 y - 2y^2 \end{cases}" height="202" src="https://latex.codecogs.com/png.latex?\dpi{200}&amp;space;\huge&amp;space;{\color{DarkOrange}\begin{cases}&amp;space;x_1&#39;&amp;space;=&amp;space;x_1&amp;space;x_2&amp;space;\\&amp;space;x_2&#39;&amp;space;=&amp;space;-x_2&amp;space;y&amp;space;\\&amp;space;y&#39;&amp;space;=&amp;space;x_2&amp;space;y&amp;space;-&amp;space;2y^2&amp;space;\end{cases}}" width="300"/>

We used only one new variable, so we achieved an *optimal* quadratization.

# Qbee usage

QBee implements algorithms that **take** system of ODEs with elementary functions right-hand side and
**return** *optimal monomial quadratization* - optimal quadratization constructed from monomial substitutions.

We will demonstrate usage of QBee on the example below. Other interactive examples you can find
in [examples section](examples).

### 1. Importing QBee

QBee relies on Sympy for a high-level API.

```python
import sympy as sp
from qbee import *
```

### 2. System definition

For example, we will take the **A1** system from [Swischuk et al.'2020](https://arxiv.org/abs/1908.03620)

<img alt="{\color{DarkOrange} \begin{cases} c_1&#39; = -A \exp(-E_a / (R_u T)) c_1 ^{0.2} c_2^{1.3}\\ c_2&#39; = 2c_1&#39; \\ c_3&#39; = -c_1&#39; \\ c_4&#39; = -2 c_1&#39; \end{cases}}" height="250" src="https://latex.codecogs.com/png.latex?\dpi{200}&amp;space;\huge&amp;space;{\color{DarkOrange}&amp;space;\begin{cases}&amp;space;c_1&#39;&amp;space;=&amp;space;-A&amp;space;\exp(-E_a&amp;space;/&amp;space;(R_u&amp;space;T))&amp;space;c_1&amp;space;^{0.2}&amp;space;c_2^{1.3}\\&amp;space;c_2&#39;&amp;space;=&amp;space;2c_1&#39;&amp;space;\\&amp;space;c_3&#39;&amp;space;=&amp;space;-c_1&#39;&amp;space;\\&amp;space;c_4&#39;&amp;space;=&amp;space;-2&amp;space;c_1&#39;&amp;space;\end{cases}}" width="550"/>

The parameters in the system are `A, Ea and Ru`, and the others are either state variables or inputs.
So, let's define them with the system in code:
```python
A, Ea, Ru = parameters("A, Ea, Ru")
c1, c2, c3, c4, T = functions("c1, c2, c3, c4, T")  

eq1 = -A * sp.exp(-Ea / (Ru * T)) * c1 ** 0.2 * c2 ** 1.3
system = [
    (c1, eq1),
    (c2, 2 * eq1),
    (c3, -eq1),
    (c4, -2 * eq1)
]
```

### 3. Polynomialization and Quadratization

When we work with ODEs with the right-hand side being a general continuous function, 
we utilize the following pipeline: 
```
Input system -> Polynomial System -> Quadratic System
```
and the transformations are called *polynomialization* and *quadratization* accordingly. 

The example system is not polynomial, so we use the most general method for achieving optimal monomial quadratization.

```python
# {T: 2} means than T can have a derivative of order at most two => T''
quadr_system = polynomialize_and_quadratize(system, input_der_orders={T: 2})
if quadr_system:
    quadr_system.print()
```

Sample output:

```
Introduced variables:
w_0 = exp(-Ea/(Ru*T))
w_1 = c1**0.2
w_2 = c2**1.3
w_3 = w_0*w_1
w_4 = T'/T**2
w_5 = T**(-2)
w_6 = T'/T
w_7 = 1/T
w_8 = w_0*w_1*w_2/c1
w_9 = w_0*w_1*w_2/c2

c1' = -A*w_2*w_3
c2' = -2*A*w_2*w_3
c3' = A*w_2*w_3
c4' = 2*A*w_2*w_3
w_0' = Ea*w_0*w_4/Ru
w_1' = -A*w_1*w_8/5
w_2' = -13*A*w_2*w_9/5
w_3' = -A*w_3*w_8/5 + Ea*w_3*w_4/Ru
w_4' = T''*w_5 - 2*w_4*w_6
w_5' = -2*w_5*w_6
w_6' = T''*w_7 - w_6**2
w_7' = -w_6*w_7
w_8' = 4*A*w_8**2/5 - 13*A*w_8*w_9/5 + Ea*w_4*w_8/Ru
w_9' = -A*w_8*w_9/5 - 3*A*w_9**2/5 + Ea*w_4*w_9/Ru
```

Introduced variables are the optimal monomial quadratization.

## Papers

* Optimal Monomial Quadratization for ODE systems: [arxiv](https://arxiv.org/abs/2103.08013), [Springer](https://link.springer.com/chapter/10.1007/978-3-030-79987-8_9)

## Citation

If you find this code useful in your research, please consider citing the above paper that works best for you. 





