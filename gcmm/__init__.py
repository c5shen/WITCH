from configs import Configs
from concurrent.futures import ProcessPoolExecutor

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

def notifyError(location):
    print('Encountered an error at {}, check {}'.format(
        location, Configs.error_path))
    exit(1)
