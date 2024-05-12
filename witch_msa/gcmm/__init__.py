from witch_msa.configs import Configs
from concurrent.futures import ProcessPoolExecutor
import os, sys, inspect

'''
Customized ProcessPoolExecutor class to handle callbacks and monitor current
progress in query alignments
'''
class WITCHProcessPoolExecutor(ProcessPoolExecutor):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._running_jobs = 0
        self._submitted_jobs = 0
        self._finished_jobs = 0
    
    def submit(self, *args, **kwargs):
        future = super().submit(*args, **kwargs)
        self._running_jobs += 1
        self._submitted_jobs += 1
        #future.add_done_callback(self._worker_is_done)
        return future

    def _worker_is_done(self, future):
        self._running_jobs -= 1
        self._finished_jobs += 1
        print('Finished jobs: {}/{}'.format(
            self._finished_jobs, self._submitted_jobs), end='\r', flush=True)
    
    def get_pool_usage(self):
        return self._running_jobs
    
    def get_finished_jobs(self):
        return self._finished_jobs

'''
Simple function for notifying user of errors encountered
'''
def notifyError(location):
    print('Encountered an error at {}\n\tcheck {}'.format(
        location, Configs.error_path))
    exit(1)

'''
Simple function to obtain current line number of the caller
'''
def getLineInfo():
    items = inspect.stack()[1][1:4]
    return '{}:{} - Line {}'.format(items[0], items[2], items[1])

'''
Simple function for sanity-checking all output files of a given list that:
    (1) they exists
    (2) they have size > 0
'''
def sanityCheckFileCreation(files):
    ret = []    # list of problematic files
    for f in files:
        if os.path.exists(f) and os.stat(f).st_size > 0:
            pass
        else:
            ret.append(f)
    return ret
