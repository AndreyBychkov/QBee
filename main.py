from qbee.examples import *

if __name__ == '__main__':
    R, x = sp.ring('x', sp.QQ)
    system = [x ** 5, ]
    res = quadratize(system)
    print(res)
