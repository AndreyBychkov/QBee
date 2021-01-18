import os


def top_priority():
    """ Set the priority of the process to above-normal."""

    if os.name == 'posix':
        os.nice(19)
    else:
        import win32api, win32process, win32con

        pid = win32api.GetCurrentProcessId()
        handle = win32api.OpenProcess(win32con.PROCESS_ALL_ACCESS, True, pid)
        win32process.SetPriorityClass(handle, win32process.REALTIME_PRIORITY_CLASS)


def platform_python():
    if os.name == 'posix':
        return "python3"
    else:
        return "python"


def call_x_plus_1_series():
    os.system(
        rf"{platform_python()} -m pytest x_plus_1_series.py --benchmark-only --benchmark-group-by=param:ord --benchmark-autosave")


def call_circular():
    os.system(
        rf"{platform_python()} -m pytest circular.py --benchmark-only --benchmark-group-by=param:ord --benchmark-autosave")


def call_hard():
    os.system(
        rf"{platform_python()} -m pytest hard.py --benchmark-only --benchmark-group-by=param:ord --benchmark-autosave")


def call_hill():
    os.system(
        rf"{platform_python()} -m pytest hill.py --benchmark-only --benchmark-group-by=param:ord --benchmark-autosave")


def call_long_monomial():
    os.system(
        rf"{platform_python()} -m pytest long_monomial.py --benchmark-only --benchmark-group-by=param:ord --benchmark-autosave")


if __name__ == '__main__':
    top_priority()

    call_x_plus_1_series()
    call_circular()
    call_hill()

    # Takes too long by now
    # call_hard()
    # call_long_monomial()
