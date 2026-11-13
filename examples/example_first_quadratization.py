"""
Stopping the search at the first quadratization found (pruning_by_quadratizations_found), on the 39-state
cardiovascular-respiratory model with 10 inputs of systems.generate_cardiovascular_respiratory.

Inputs are piecewise constant, so new variables must not contain them (input_free=True). The default search for an
optimal quadratization finds nothing within 20 minutes on this system. Stopped at its first quadratization it returns
one with 26 new variables after 139 nodes (about 3 minutes); with the preliminary upper-bound search as well
(calc_upper_bound=True), one with 24 after about 31000 nodes (about 8 minutes).
"""
import time

from qbee import *
from systems import generate_cardiovascular_respiratory

if __name__ == '__main__':
    system = generate_cardiovascular_respiratory()
    print(f"{len(system)} states")
    t0 = time.time()
    res = quadratize(system, input_free=True, calc_upper_bound=False,
                     pruning_functions=[pruning_by_quadratizations_found, *default_pruning_rules])
    print(f"First quadratization: {len(res.introduced_variables)} new variables, {res.nodes_traversed} nodes, "
          f"{time.time() - t0:.0f} s")
    for eq in res.introduced_variables:
        print(str_qbee(eq))
