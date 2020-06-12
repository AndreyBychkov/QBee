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







