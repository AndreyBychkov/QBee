import os

if __name__ == '__main__':
    os.system(rf"pytest-benchmark compare 0005 --group-by=param:ord")
