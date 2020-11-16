import os

def platform_python():
    if os.name == 'posix':
        return "python3"
    else:
        return "python"

def call_x_plus_1_series():
    os.system(rf"{platform_python()} -m pytest x_plus_1_series.py --benchmark-only --benchmark-group-by=param:ord")

if __name__ == '__main__':
    call_x_plus_1_series()