from qbee import *
from qbee.examples import *

if __name__ == '__main__':
    R, x, y = sp.ring('x, y', sp.QQ)
    system = [y ** 3, x ** 3]
    res = quadratize(system)
    print(res)
