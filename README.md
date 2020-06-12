# QBee
Python library for converting an ODE system into a quadratic form.

# Installation

1. Clone repository: `https://github.com/AndreyBychkov/QBee.git` 
2. Install requirements: `pip install -r requirements.txt`
3. Install the package: `pip install .`

# What is quadratization?

Consider a system of ODEs

![\begin{cases}     \dot{x}_1 = f_1(\bar{x}),\\     \ldots\\     \dot{x}_n = f_n(\bar{x}),   \end{cases}](https://render.githubusercontent.com/render/math?math=%5Cbegin%7Bcases%7D%20%20%20%20%20%5Cdot%7Bx%7D_1%20%3D%20f_1(%5Cbar%7Bx%7D)%2C%5C%5C%20%20%20%20%20%5Cldots%5C%5C%20%20%20%20%20%5Cdot%7Bx%7D_n%20%3D%20f_n(%5Cbar%7Bx%7D)%2C%20%20%20%5Cend%7Bcases%7D)

where ![\bar{x} = (x_1, \ldots, x_n)](https://render.githubusercontent.com/render/math?math=%5Cbar%7Bx%7D%20%3D%20(x_1%2C%20%5Cldots%2C%20x_n))
and ![f_1, \ldots, f_n \in \mathbb{C}\[\mathbf{x}\]](https://render.githubusercontent.com/render/math?math=f_1%2C%20%5Cldots%2C%20f_n%20%5Cin%20%5Cmathbb%7BC%7D%5B%5Cmathbf%7Bx%7D%5D).
Then a list of new variables 

![y_1 = g_1(\bar{x}), \ldots, y_m = g_m(\bar{x})](https://render.githubusercontent.com/render/math?math=y_1%20%3D%20g_1(%5Cbar%7Bx%7D)%2C%20%5Cldots%2C%20y_m%20%3D%20g_m(%5Cbar%7Bx%7D))

is said to be a *quadratization* of ODEs system if there exists polynomials ![h_1, \ldots, h_{m + n} \in \mathbb{C}\[\bar{x}, \bar{y}\]](https://render.githubusercontent.com/render/math?math=h_1%2C%20%5Cldots%2C%20h_%7Bm%20%2B%20n%7D%20%5Cin%20%5Cmathbb%7BC%7D%5B%5Cbar%7Bx%7D%2C%20%5Cbar%7By%7D%5D)
of degree at most two such that

* ![\dot{x}_i = h_i(\bar{x}, \bar{y})$ for every $1 \leqslant i \leqslant n](https://render.githubusercontent.com/render/math?math=%5Cdot%7Bx%7D_i%20%3D%20h_i(%5Cbar%7Bx%7D%2C%20%5Cbar%7By%7D)%24%20for%20every%20%241%20%5Cleqslant%20i%20%5Cleqslant%20n)
* ![\dot{y}_j = h_{j + n}(\bar{x}, \bar{y})$ for every $1 \leqslant j \leqslant m](https://render.githubusercontent.com/render/math?math=%5Cdot%7By%7D_j%20%3D%20h_%7Bj%20%2B%20n%7D(%5Cbar%7Bx%7D%2C%20%5Cbar%7By%7D)%24%20for%20every%20%241%20%5Cleqslant%20j%20%5Cleqslant%20m)

The number *m* will be called *the order of quadratization*.

A quadratization of the smallest possible order will be called  *an optimal quadratization*.






