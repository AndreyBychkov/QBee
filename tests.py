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
    return "python3" if os.name == "posix" else "python"


def benchmark_everything(save_results=True, make_histogram=True, sort_by="name", min_rounds=10, workers=1):
    query = [f"{platform_python()} -m pytest -v",
             "--benchmark-only",
             "--benchmark-group-by=param:ord",
             "--benchmark-warmup=on",
             f"--benchmark-min-rounds={min_rounds}",
             f"--benchmark-sort={sort_by}"]
    if save_results:
        query.append("--benchmark-autosave")
    if make_histogram:
        query.append("--benchmark-histogram")

    os.system(" ".join(query))


def run_tests():
    top_priority()
    os.system(f'{platform_python()} -m pytest -v -m "not benchmark and not experimental and not expensive"')


def run_benchmarks():
    top_priority()
    benchmark_everything()


if __name__ == '__main__':
    run_tests()
