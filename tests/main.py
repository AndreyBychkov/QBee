import os


def run_without_benchmark_systems():
    os.system('python3 -m pytest -m "not benchmark and not experimental"')


def run_all():
    os.system('python3 -m pytest')


if __name__ == '__main__':
    run_without_benchmark_systems()
