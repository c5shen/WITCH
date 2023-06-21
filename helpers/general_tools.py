import psutil, os

# return memory usage of python process by MB
def memoryUsage():
    process = psutil.Process(os.getpid())
    mem = process.memory_info().rss / float(2 ** 20)
    return mem
