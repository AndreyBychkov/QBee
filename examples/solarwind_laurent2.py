from sympy import *
from qbee import *

if __name__ == '__main__':
    
    R, v1, v2, v3, v4, v5, v6 = ring("v_1, v_2, v_3, v_4, v_5, v_6", QQ)
    #vt = vx * vinv - c * vx
    v1t = (v1 + v2 + v3) * v1**(-1) - (v1 + v2 + v3)
    v2t = (v2) * v2**(-1) + (v2)
    v3t = (v3) * v3**(-1) + (v3)
    v4t = (v4 + v5) * v4**(-1) - (v4 + v5)
    v5t = (v5 + v6) * v5**(-1) - (v5 + v6)
    v6t = v6 * v6**(-1) - v6


    quadr_system = quadratize([v1t, v2t, v3t, v4t, v5t, v6t], generation_strategy=generation_semidiscretized)
    if quadr_system:
        print(quadr_system)
