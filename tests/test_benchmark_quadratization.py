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
