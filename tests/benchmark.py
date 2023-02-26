import os
import multiprocessing as mp
from psutil import cpu_count
from glob import glob


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


def benchmark_file(filename, save_results=True, make_histogram=True, sort_by="name", min_rounds=5):
    query = [f"{platform_python()} -m pytest {filename}",
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


def benchmark_everything(save_results=True, make_histogram=True, sort_by="name", min_rounds=10, workers=1):
    # TODO: spawn processes to accelerate
    # TODO: Rewrite with pytest -m "benchmark" or smth like this
    [benchmark_file(file, save_results, make_histogram, sort_by, min_rounds)
     for file in glob(r"*.py") if file not in ["main.py", "__init__.py"]]


if __name__ == '__main__':
    top_priority()
    benchmark_everything()
