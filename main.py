from qbee.examples import *
def top_priority():
    """ Set the priority of the process to above-normal."""
    import os
    if os.name == 'posix':
        os.nice(19)
    else:
        import win32api, win32process, win32con

        pid = win32api.GetCurrentProcessId()
        handle = win32api.OpenProcess(win32con.PROCESS_ALL_ACCESS, True, pid)
        win32process.SetPriorityClass(handle, win32process.REALTIME_PRIORITY_CLASS)

def x6x4x3():
    R, x = sp.ring('x', sp.QQ)
    system = PolynomialSystem([x ** 6 + x ** 4 + x ** 3])
    algo = BranchAndBound(system, aeqd_score, [termination_by_best_nvars])
    results = algo.get_optimal_quadratizations()
    for res in results:
        print(res)


def circ4():
    system = generate_circular(4)
    algo = BranchAndBound(system, aeqd_score, [termination_by_best_nvars])
    results = algo.get_optimal_quadratizations()
    for res in results:
        print(res)


def hard3():
    system = generate_hard(3)
    algo = BranchAndBound(system, aeqd_score, [termination_by_best_nvars])
    results = algo.get_optimal_quadratizations()
    for res in results:
        print(res)


if __name__ == '__main__':
    top_priority()
    R, x, y = ring(['x', 'y'], QQ)
    system = generate_lifeware_conjecture(6)
    algo = BranchAndBound(system, aeqd_score, [termination_by_best_nvars])
    algo.newton_polygon_vertices_upper_bound()
    quad_res = algo.quadratize()
    print(quad_res)
