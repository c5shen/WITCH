import psutil, os
import argparse

# return memory usage of python process by MB
def memoryUsage():
    process = psutil.Process(os.getpid())
    mem = process.memory_info().rss / float(2 ** 20)
    return mem

# reformat argparse help text formatting
class SmartHelpFormatter(argparse.HelpFormatter):
    def _split_lines(self, text, width):
        if '\n' in text:
            temp = text.split('\n')
            ret = []
            for _splice in [argparse.HelpFormatter._split_lines(self, x, width)
                    for x in temp]:
                ret.extend(_splice)
            return ret
        return argparse.HelpFormatter._split_lines(self, text, width)
