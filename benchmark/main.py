import pandas as pd
import multiprocessing as mp
import subprocess
from qbee import *
from tqdm.autonotebook import tqdm
from typing import List


def benchmark_system(system: PolynomialSystem,
                     system_name,
                     cycles=10,
                     cpu_count=mp.cpu_count() - 2,
                     search_algorithms=('BFS', 'ID-DLS', 'MMDR'),
                     heuristics=('none', 'FF', 'FVC', 'AED', 'AEQD', 'SMD'),
                     initial_max_depth=1,
                     limit_depth=None):
    res = list()
    for algo in search_algorithms:
        for heur in heuristics:
            steps = run_benchmark(system, cycles, system_name, cpu_count=cpu_count, search_algorithm=algo, heuristics=heur, initial_max_depth=initial_max_depth,
                                  limit_depth=limit_depth)
            res_row = {"algorithm": algo, "heuristics": heur, 'steps': steps}
            res.append(res_row)
    res_df = pd.DataFrame(res)
    res_df.to_csv(f"{system_name}_benchmark.csv", index=False)


def run_benchmark(system: PolynomialSystem, cycles: int, system_name, cpu_count=mp.cpu_count() - 2, **params) -> List[int]:
    total_steps = list()
    for _ in tqdm(range(cycles), desc=f"{system_name}: {params['search_algorithm']} | {params['heuristics']}", unit="cycle"):
        with mp.Pool(cpu_count) as pool:
            curr_steps = pool.map(quad_steps, [(system, params)] * cpu_count)
            total_steps += curr_steps
    return total_steps


def quad_steps(system_params):
    system, params = system_params
    return quadratize(system, **params).statistics.steps


if __name__ == '__main__':
    subprocess.call(["python", "xSigmoid.py"], cwd="xSigmoid")
    subprocess.call(["python", "RabinovichFabrikant.py"], cwd="RabinovichFabrikant")
    subprocess.call(["python", "TwoParameterBlueSkyCatastrophe.py"], cwd="TwoParameterBlueSkyCatastrophe")
    # subprocess.call(["python", "4DQiDuChenChenYuan.py"], cwd="4DQiDuChenChenYuan")
