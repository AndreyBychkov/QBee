import os


def run_default():
    os.system('python3 -m pytest -m "not benchmark and not experimental and not expensive"')


def run_all():
    os.system('python3 -m pytest')


if __name__ == '__main__':
    run_default()
