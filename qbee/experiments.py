import sympy as sp
import networkx as nx
import multiprocessing as mp
import platform
from textwrap import dedent
from typing import Optional, Collection, Dict as DictT
from termcolor import colored
from functools import partial
from scipy.spatial import Delaunay, ConvexHull
from sympy.polys.orderings import monomial_key
from collections import OrderedDict
from .heuristics import aeqd_score
from .quadratization import BranchAndBound, PolynomialSystem, pruning_by_best_nvars
from .util import monom2str
from .examples import *


def get_examples():
    examples = dict()
    for i in [3, 4, 5, 6, 7]:
        examples[f'Circular {i}'] = generate_circular(i)
    for i in [2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20]:
        examples[f'Hill {i}'] = generate_hill(i)
    for i in [2, 3, 4]:
        examples[f'Hard {i}'] = generate_hard(i)
    for i in [2, 3]:
        examples[f'Lifeware Conjecture {i}'] = generate_lifeware_conjecture(i)
    for i in [2, 3, 4, 5, 6, 7]:
        examples[f'Cubic Cycle {i}'] = generate_cubic_cycle(i)
    # for i in [6, 7, 8]:
    #     examples[f'Cubic Bicycle {i}'] = generate_cubic_bicycle(i)
    return examples


Seconds = float


def make_benchmark_report(quad_funcs: DictT[str, Callable[[PolynomialSystem], Seconds]],
                          examples=None,
                          n_samples=mp.cpu_count(),
                          use_multiprocessing=True,
                          n_jobs=mp.cpu_count() - 2,
                          title="Benchmark report"):
    if examples is None:
        examples = get_examples()

    def eval_no_mp(func, system, results):
        if n_samples > 1:
            times = [func(system) for i in range(n_samples)]
            min_t, avg_t, std_t = np.min(times), np.average(times), np.std(times)
            results.append((min_t, avg_t, std_t))
        else:
            time_s = func(system)
            results.append(time_s)

    def eval_mp(func, system, results):
        with mp.Pool(n_jobs) as pool:
            times = pool.map(func, [system, ] * n_samples)
            min_t, avg_t, std_t = np.min(times), np.average(times), np.std(times)
            results.append((min_t, avg_t, std_t))

    top_priority()
    report = BenchmarkReport(title, quad_funcs.keys())
    for sys_name, system in examples.items():
        print(f"Processing {sys_name}...")
        results = list()
        for func_name, func in quad_funcs.items():
            if not use_multiprocessing:
                eval_no_mp(func, system, results)
            else:
                eval_mp(func, system, results)
        report.add_row(sys_name, results)
    report.save_to_file()


class BenchmarkReport:
    def __init__(self, title, algo_names):
        self.title = title
        self.content = dedent(f"""
        # {title}
        
        ---
        
        |    System    |    {'    |    '.join(algo_names)}    |
        |{':----------:|' * (len(algo_names) + 1)}""")

    def add_row(self, system, algo_evals):
        if self._is_single_record(algo_evals):
            self.content += dedent(f"""
            |   {system}   |   {'   |   '.join(list(map(self._format_time_no_mp, algo_evals)))}   |""")
        else:
            self.content += dedent(f"""
            |   {system}   |   {'   |   '.join(list(map(self._format_time_mp, algo_evals)))}   |""")

    def save_to_file(self):
        self._add_system_info()
        print(self.content)
        with open(f"{self.title}.md", 'w') as file:
            file.write(self.content)

    def _add_system_info(self):
        self.content += dedent(f"""
        
        ---
        
        * Full info: {platform.platform(True)}
        * CPU: {platform.processor()}
        * System: {platform.system()}
        * Python version: {platform.python_implementation()} {platform.python_version()}
        """)

    def _is_single_record(self, rest_list):
        return isinstance(rest_list[0], float)

    def _format_time_no_mp(self, time):
        if time < 1:
            return f'{np.round(time * 1000, 1)}ms'
        else:
            return f'{np.round(time, 1)}s'

    def _format_time_mp(self, time_stat):
        min_t, avg_t, std_t = time_stat
        min_format = self._format_time_no_mp(min_t)
        avg_format = self._format_time_no_mp(avg_t)
        if avg_t < 1:
            std_format = f'{np.round(std_t * 1000, 1)}'
        else:
            std_format = f'{np.round(std_t, 1)}'
        return f"{min_format}<br>{avg_format}+-{std_format}"


# ======================================================================

def make_graphs(n=2, i=None, j=None):
    """Make a graph recursively, by either including, or skipping each edge.

    Edges are given in lexicographical order by construction."""
    out = []
    if i is None:  # First call
        out = [[(0, 1)] + r for r in make_graphs(n=n, i=0, j=1)]
    elif j < n - 1:
        out += [[(i, j + 1)] + r for r in make_graphs(n=n, i=i, j=j + 1)]
        out += [r for r in make_graphs(n=n, i=i, j=j + 1)]
    elif i < n - 1:
        out = make_graphs(n=n, i=i + 1, j=i + 1)
    else:
        out = [[]]
    return out
