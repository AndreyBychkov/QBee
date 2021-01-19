[![Build Status](https://travis-ci.com/AndreyBychkov/QBee.svg?branch=master)](https://travis-ci.com/AndreyBychkov/QBee)

# QBee
Python library for converting an ODE system into a quadratic form.

# Installation

1. Clone repository: `https://github.com/AndreyBychkov/QBee.git` 
2. Install requirements: `pip install -r requirements.txt`
3. Install the package: `pip install .`

# What is quadratization?

The problem of *quadratization* is, given a system of ODEs with polynomial right-hand side, reduce the system to a system with quadratic right-hand side by introducing as few new variables as possible.
We will explain it using toy example. Consider the system

![\begin{cases} \dot{x}_1 = x_1 x_2 \\ \dot{x}_2 = -x_1 x_2^3 \end{cases}](https://render.githubusercontent.com/render/math?math=%5Cbegin%7Bcases%7D%20%5Cdot%7Bx%7D_1%20%3D%20x_1%20x_2%20%5C%5C%20%5Cdot%7Bx%7D_2%20%3D%20-x_1%20x_2%5E3%20%5Cend%7Bcases%7D)

An example of quadratization of this system will be a new variable

![y = x_1 x_2^2](https://render.githubusercontent.com/render/math?math=y%20%3D%20x_1%20x_2%5E2)

leading to the following ODE

![\dot{y} = x_2 y - 2y^2](https://render.githubusercontent.com/render/math?math=%5Cdot%7By%7D%20%3D%20x_2%20y%20-%202y%5E2)

Thus, we attained the system with quadratic right-hand side

![\begin{cases} \dot{x}_1 = x_1 x_2 \\ \dot{x}_2 = -x_2 y \\ \dot{y} = x_2 y - 2y^2 \end{cases}](https://render.githubusercontent.com/render/math?math=%5Cbegin%7Bcases%7D%20%5Cdot%7Bx%7D_1%20%3D%20x_1%20x_2%20%5C%5C%20%5Cdot%7Bx%7D_2%20%3D%20-x_2%20y%20%5C%5C%20%5Cdot%7By%7D%20%3D%20x_2%20y%20-%202y%5E2%20%5Cend%7Bcases%7D)

We used only one new variable, so we achieved an *optimal* quadratization. 

# Qbee usage

QBee implements algorithms that **take** system of ODEs with elementary functions right-hand side and
**return** *optimal monomial quadratization* - optimal quadratization constructed from monomial substitutions.

We will demonstrate usage of QBee on the example below. Other interactive examples you can find in [examples section](old_examples). 

### 1. Importing QBee

```python
import sympy
from qbee import *

sympy.init_printing()  # If you work in Jupyter notebook 
```

### 2. General ODEs and polynomialization

First, we introduce variables, using SymPy framework
```python
x = sympy.symbols('x')
dx = derivatives(x)
```

Then we build a ODEs system with elementary function right-hand side
```python
system = EquationSystem([
    sympy.Eq(dx, 1 / (1 + sympy.exp(x)))
])
```

This system is not polynomial, so we need to polynomialize it.

---

We can convert equations right-hand side to polynomials
```python
poly_system = polynomialize(system)
poly_system.print()
```
Output:
```
x' = y_{1}
y_{0}' = y_{0}*y_{1}
y_{1}' = -y_{0}*y_{1}**3
```

**Note:** our implementation of polynomialization is **not** optimal yet. 

### 3. Polynomial ODEs and quadratization

Now, convert polynomial system to quadratic form

First, we define a system of ODEs.
In the list inside `PolynomialSystem` defined right-hand side of input ODE in linear order.
For example, consider the system:

![\begin{cases} \dot{x} = y^3 \\ \dot{y} = x^3 \end{cases}](https://render.githubusercontent.com/render/math?math=%5Cbegin%7Bcases%7D%20%5Cdot%7Bx%7D_1%20%3D%20x_1%20x_2%20%5C%5C%20%5Cdot%7Bx%7D_2%20%3D%20-x_2%20y%20%5C%5C%20%5Cdot%7By%7D%20%3D%20x_2%20y%20-%202y%5E2%20%5Cend%7Bcases%7D)


```python
x, y, _ = sympy.ring(['x', 'y'], sympy.QQ)
system = PolynomialSystem([y ** 3, x **3])
```

Then we need to choose an algorithm. We suggest to use Branch & Bound.

```python
algo = BranchAndBound(system, aeqd_score, [termination_by_best_nvars])
```

Here, we explicitly define heuristics `aeqd_score` and single-value list of early termination criteria `termination_by_best_nvars`.
You can vary these parameters to fit to your system. 

Now, we need to start the algorithm:
```python
result = algo.quadratize()
print(result)
```

Sample output:
```
Elapsed time: 0.209s.
Number of introduced variables: 3
Introduced variables: ['x y', 'x^{2}', 'y^{2}']
Nodes traversed: 17
```

Introduced variables are optimal monomial quadratization.

### 4. Work inside of package

#### 1. Configuration

Inside of `config.ini` you can change the following arguments:
* `logging_enable = [True | False]`. If enabled, work of algorithm is logged into `logging_file` and `quad_systems_file`. Requires memory to work. Is not recommended for long quadratizations.
* `logging_file`: must be in Apache Arrow `.feather` format.
* `quad_systems_file`: dump quadratic systems by using pickle. `.pkl` file format is recommended.
* `progress_bar_enable`: enables progress bar during quadratization.

#### 2. Visualization

In order to visualize work of an algorithm you can pass logging data to `qbee.visualize.visualize_pyvis`:
```python
visualize_pyvis('log/log.feather', 'quad_systems.pkl')
```






