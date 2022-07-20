[![Build Status](https://travis-ci.com/AndreyBychkov/QBee.svg?branch=dev)](https://app.travis-ci.com/github/AndreyBychkov/QBee)
[![Documentation Status](https://readthedocs.org/projects/qbee/badge/?version=latest)](https://qbee.readthedocs.io/en/latest/?badge=latest)
# QBee

Python library for transforming systems of ODE equations into a systems with quadratic right-rand side.

# Installation

1. Clone repository: `https://github.com/AndreyBychkov/QBee.git`
   * Or, if you want our bleeding edge version, clone `https://github.com/AndreyBychkov/QBee/tree/dev`
2. Install requirements: `pip install -r requirements.txt`
3. Install the package: `pip install .`

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
import sympy
from qbee import *

sympy.init_printing()  # If you work in Jupyter notebook 
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
    print("Quadratized system:")
    print(quadr_system)
```

Sample output:

```
Variables introduced in polynomialization:
w_{0} = c1**(-0.8)
w_{1} = c2**(-0.7)
w_{2} = 1/T
w_{3} = exp(-Ea*w_{2}/Ru)

Elapsed time: 0.139s.
==================================================
Quadratization result
==================================================
Number of introduced variables: 5
Nodes traversed: 117
Introduced variables:
w_{4} = T'*w_{2}
w_{5} = T'*w_{2}**2
w_{6} = c2**2*w_{0}*w_{1}*w_{3}
w_{7} = w_{2}**2
w_{8} = c1*c2*w_{0}*w_{1}*w_{3}

Quadratized system:
c1' = -A*c2*w_{8}
c2' = -2*A*c2*w_{8}
c3' = A*c2*w_{8}
c4' = 2*A*c2*w_{8}
w_{0}' = 4*A*w_{0}*w_{6}/5
w_{1}' = 7*A*w_{1}*w_{8}/5
w_{2}' = -T'*w_{7}
w_{3}' = Ea*w_{3}*w_{5}/Ru
T' = T'
T'' = T''
T''' = 0
w_{4}' = -T'*w_{5} + T''*w_{2}
w_{5}' = T''*w_{7} - 2*w_{4}*w_{5}
w_{6}' = 4*A*w_{6}**2/5 - 13*A*w_{6}*w_{8}/5 + Ea*w_{5}*w_{6}/Ru
w_{7}' = -2*w_{4}*w_{7}
w_{8}' = -A*w_{6}*w_{8}/5 - 3*A*w_{8}**2/5 + Ea*w_{5}*w_{8}/Ru

Process finished with exit code 0

```

Introduced variables are the optimal monomial quadratization.

### 4. Work inside of package

#### 1. Configuration

Inside of `config.ini` you can change the following arguments:

* `logging_enable = [True | False]`. If enabled, work of algorithm is logged into `logging_file` and `quad_systems_file`
  . Requires memory to work. Is not recommended for long quadratizations.
* `logging_file`: must be in Apache Arrow `.feather` format.
* `quad_systems_file`: dump quadratic systems by using pickle. `.pkl` file format is recommended.
* `progress_bar_enable`: enables progress bar during quadratization.

#### 2. Visualization

In order to visualize work of an algorithm you can pass logging data to `qbee.visualize.visualize_pyvis`:

```python
visualize_pyvis('log.feather', 'quad_systems.pkl')
```

## Papers

* Optimal Monomial Quadratization for ODE systems: [arxiv](https://arxiv.org/abs/2103.08013), [Springer](https://link.springer.com/chapter/10.1007/978-3-030-79987-8_9)

## Citation

If you find this code useful in your research, please consider citing the above paper that works best for you. 





