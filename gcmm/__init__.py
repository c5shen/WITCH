from configs import Configs

def notifyError(location):
    print('Encountered an error at {}, check {}'.format(
        location, Configs.error_path))
    exit(1)
