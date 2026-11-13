import pytest
from qbee import *
from examples import *
from functools import partial


@pytest.mark.benchmark
@pytest.mark.parametrize('ord', [3, 4, 5, 6])
def test_circular(benchmark, ord):
    benchmark(partial(quadratize, generate_circular(ord)))


@pytest.mark.benchmark
@pytest.mark.parametrize('ord', [3, 4])
def test_hard(benchmark, ord):
    benchmark(partial(quadratize, generate_hard(ord)))


@pytest.mark.benchmark
@pytest.mark.parametrize('ord', [3, 6, 8, 10, 15])
def test_hill(benchmark, ord):
    benchmark(partial(quadratize, generate_hill(ord)))


@pytest.mark.benchmark
@pytest.mark.parametrize('ord', [2, 3])
def test_lifeware_conjecture(benchmark, ord):
    benchmark(partial(quadratize, generate_lifeware_conjecture(ord)))


@pytest.mark.benchmark
@pytest.mark.parametrize('ord', [6, 8, 10, 12])
def test_x_plus_1(benchmark, ord):
    benchmark(partial(quadratize, generate_x_plus_1(ord)))


def _first_quadratization(system):
    return quadratize(system, input_free=True, calc_upper_bound=False,
                      pruning_functions=[pruning_by_quadratizations_found, *default_pruning_rules])


@pytest.mark.benchmark
@pytest.mark.parametrize('kwargs', [dict(gas_exchange=False), dict(n_sites=1, drug_conductance=False)])
def test_cardiovascular_respiratory_first(benchmark, kwargs):
    # 11 states (about 4 s) and 27 states (about 40 s)
    system = generate_cardiovascular_respiratory(**kwargs)
    benchmark.pedantic(_first_quadratization, args=(system,), rounds=1)


@pytest.mark.benchmark
@pytest.mark.expensive
def test_cardiovascular_respiratory_first_full(benchmark):
    # 39 states, 10 inputs: about 3 minutes; the search for an optimal quadratization does not finish in 20 minutes
    benchmark.pedantic(_first_quadratization, args=(generate_cardiovascular_respiratory(),), rounds=1)
