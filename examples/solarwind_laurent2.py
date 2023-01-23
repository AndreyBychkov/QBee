from sympy import *
from qbee import *

if __name__ == '__main__':
    
    R, v1, v2, v3 = ring("v_1, v_2, v_3", QQ)
    #vt = vx * vinv - c * vx
    v1t = (v1 + v2) * v1**(-1) - (v1 + v2)
    v2t = (v1 + v2 + v3) * v2**(-1) + (v1 + v2 + v3)
    v3t = (v2 + v3) * v3**(-1) + (v2 + v3)

    quadr_system = quadratize([v1t, v2t, v3t], generation_strategy=generation_semidiscretized)
    if quadr_system:
        print(quadr_system)
