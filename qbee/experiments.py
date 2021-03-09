import sympy as sp
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
from .quadratization import BranchAndBound, PolynomialSystem, termination_by_best_nvars
from .util import monom2str, get_hull_vertices
from .examples import *


def show_convex_hull_report(system: PolynomialSystem, max_depth: Optional[int] = None):
    if max_depth is None:
        algo = BranchAndBound(system, aeqd_score, [termination_by_best_nvars])
        max_depth = algo.quadratize().introduced_vars + __algo_hull_size(algo)
        print(f"Upper bound of search is {max_depth}")

    algo = BranchAndBound(system, aeqd_score, [termination_by_best_nvars])
    quads = __filter_duplicates(algo.get_quadratizations(max_depth))

    print(
        f"Order of quadratization: {colored('green', 'green')} if monom is in convex hull and {colored('red', 'red')} if it's not.")
    for quad in sorted(quads, key=PolynomialSystem.new_vars_count):
        is_in_hull = __points_in_hull(quad, algo.hull)
        colored_quads = __color_quad_in_hull(is_in_hull)

        print(f"{quad.new_vars_count()}: ", end="")
        for q in colored_quads:
            print(q, end=', ')
        print()

    print("=" * 50)

    print(
        f"Order of quadratization: {colored('green', 'green')} if monom is vertex of convex hull and {colored('red', 'red')} if it's not.")
    for quad in sorted(quads, key=PolynomialSystem.new_vars_count):
        belong_to_hull = __is_hull_vertices(quad, algo.hull)
        colored_quads = __color_quad_in_hull(belong_to_hull)

        print(f"{quad.new_vars_count()}: ", end="")
        for q in colored_quads:
            print(q, end=', ')
        print(f"[{len(list(filter(lambda x: x, belong_to_hull.values())))}/{__algo_hull_size(algo)}]")


def __filter_duplicates(quads: Collection[PolynomialSystem]):
    res = list()
    intr_vars = list()
    for quad in quads:
        if quad.introduced_vars not in intr_vars:
            res.append(quad)
            intr_vars.append(quad.introduced_vars)
    return res


def __color_quad_in_hull(is_in_hull: dict):
    res = list()
    for s, pred in is_in_hull.items():
        if pred:
            res.append(colored(s, 'green'))
        else:
            res.append(colored(s, 'red'))
    return res


def __points_in_hull(quad: PolynomialSystem, hull: ConvexHull):
    variables = __sorted_vars(quad)
    return dict(zip(
        list(map(partial(monom2str, gens=quad.gen_syms), variables)),
        Delaunay(hull.points).find_simplex(variables) >= 0))


def __is_hull_vertices(quad: PolynomialSystem, hull: ConvexHull):
    vertices = get_hull_vertices(hull)
    variables = __sorted_vars(quad)
    return dict(zip(
        list(map(partial(monom2str, gens=quad.gen_syms), variables)),
        list(map(lambda m: list(m) in vertices, variables))))


def __algo_hull_size(algo: BranchAndBound):
    return len(get_hull_vertices(algo.hull))


def __sorted_vars(quad: PolynomialSystem):
    return sorted(quad.introduced_vars,
                  key=lambda m: monomial_key('grlex', reversed(quad.gen_syms))(sp.Monomial(m, quad.gen_syms).as_expr()))


# ======================================================================


def get_examples():
    examples = dict()
    for i in [3, 4, 5, 6, 7, 8]:
        examples[f'Circular {i}'] = generate_circular(i)
    for i in [2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20]:
        examples[f'Hill {i}'] = generate_hill(i)
    for i in [2, 3, 4]:
        examples[f'Hard {i}'] = generate_hard(i)
    for i in [2, 3]:
        examples[f'Lifeware Conjecture {i}'] = generate_lifeware_conjecture(i)
    for i in [2, 3, 4, 5, 6, 7]:
        examples[f'Cubic Cycle {i}'] = generate_cubic_cycle(i)
    for i in [6, 7, 8]:
        examples[f'Cubic Bicycle {i}'] = generate_cubic_bicycle(i)
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
