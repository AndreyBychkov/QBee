from qbee import *
from qbee.examples import *

if __name__ == '__main__':
    R, x, c = sp.ring("x, c", sp.QQ)
    quad_system = quadratize([c  * x**3], [c])
    print(quad_system)