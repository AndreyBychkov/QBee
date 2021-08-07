from qbee import *
from qbee.examples import *
from copy import deepcopy

if __name__ == '__main__':
    system = generate_exponential()
    quad_res = quadratize(system, new_vars_name="z")
    print(quad_res)
