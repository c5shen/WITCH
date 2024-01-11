import time
from witch_msa.configs import Configs, tqdm_styles
from tqdm import tqdm
import itertools

import concurrent.futures

'''
A class defining a generic task object that will be used for job submission 
'''
class MyTask(object):
    # (required) the list of parameters used in the task
    # (optional) parent(s) of the task
    #               (i.e., the other tasks that depends on this task)
    # (optional) children of this task
    #               (i.e., the other tasks that this task depends on)
    def __init__(self, *args, **kwargs):
        self.args = tuple(*args)

        _valid = ['parent', 'children']
        for k, v in kwargs.items():
            if k in _valid:
                if isinstance(v, MyTask):
                    setattr(self, k, [v])
                elif isinstance(v, list):
                    setattr(self, k, v)
                else:
                    raise TypeError(type(v))
    
    def get_args(self):
        return self.args
    
    # currently unused
    def get_parent(self):
        if 'parent' in self.__dict__:
            return self.parent
        else:
            return None
    
    # currently unused
    def get_children(self):
        if 'children' in self.__dict__:
            return self.children
        else:
            return None

'''
Helper function to convert a list of lists of arguments to a list of MyTask objects
Assumption: all elements of args are of the same length
Return:     a generator of MyTask
'''
def getTasks(*args):
    total_length = len(args[0])
    for i in range(total_length):
        _args = [x[i] for x in args]
        yield MyTask(_args)

'''
Variant of getTasks to use given index positions to select
Also append the index to each yielded element at back
'''
def getTasksWithIndexes(indexes, *args):
    for i in indexes:
        _args = [x[i] for x in args] + [i]
        yield MyTask(_args)

'''
Helper function to handle a single future object with any return values.
Run additional callbacks with the return values and the additional callback
arguments supplemented.
Return:     the runtime to run the additional callback
'''
def handleFuture(future, success, ignored, retry, i_retry,
        callback_func, callback_args):
    s1 = time.time()
    ret = future.result()
    
    # first four fields of any callbacks will be: success=<list>,
    # ignored=<list>, retry=<list>, i_retry=<int>
    if callback_func:
        callback_func(success, ignored, retry, i_retry, 
                *ret, *callback_args)
        return time.time() - s1
    else:
        # default behavior: attach ret to success
        success.append(ret)
        return 0.

'''
Helper function to run tasks defined by a list of MyTask objects
Required:   the function that defines the task
            the process pool to submit to
            a generator of MyTask objects
            the number of MyTask objects
Optional:   max_concurrent_jobs=<int>   # default will submit all tasks at once
            i_retry=<int>               # number of retry for failed tasks
            callback_func=<function>    # callback to run after a future is handled
                                          should take future return values as the
                                          first set of arguments
            callback_args=<*args>       # additional arguments for the callback
Return:     success, ignored, retry
            total runtime (seconds) for handling/running the tasks
'''
def runTasks(func, pool, mytasks, num_tasks, **kwargs):
    handle_runtime = 0.
    max_concurrent_jobs = kwargs.get('max_concurrent_jobs', num_tasks)
    i_retry = kwargs.get('i_retry', 0) 
    callback_func = kwargs.get('callback_func', None)
    callback_args = kwargs.get('callback_args', [])

    success, ignored, retry = [], [], []
    with tqdm(total=num_tasks, **tqdm_styles) as pbar:
        futures = {
                #pool.submit(func, *task.get_args()): task.get_id()
                pool.submit(func, *task.get_args()): task.get_parent()
                for task in itertools.islice(mytasks, max_concurrent_jobs)
        }
        while futures:
            # wait for the next future to complete
            done, _ = concurrent.futures.wait(
                    futures, return_when=concurrent.futures.FIRST_COMPLETED)
            
            for future in done:
                # depending on kwargs, allow re-adding some failed tasks back
                # to queue
                handle_runtime += handleFuture(future, success, ignored, retry,
                        i_retry, callback_func, callback_args)
                _ = futures.pop(future)
            pbar.update(len(done))

            # schedule the next batch of tasks, no more than the number of tasks
            # that just finished
            for task in itertools.islice(mytasks, len(done)):
                future = pool.submit(func, *task.get_args())
                futures[future] = task.get_parent()
    return success, ignored, retry, handle_runtime
