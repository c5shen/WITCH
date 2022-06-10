import os, sys
try:
    import configparser
except ImportError:
    import ConfigParser as configparser
from argparse import ArgumentParser, Namespace
from platform import platform

_root_dir = os.path.dirname(os.path.realpath(__file__))
_config_path = os.path.join(_root_dir, 'default.config')

def setup(platform=platform()): 
    config_defaults = []
    cparser = configparser.ConfigParser()
    cparser.optionxform = str
    assert os.path.exists('{}'.format(_config_path)), \
            "default config file {} missing! Please redownload from Github\n".format(
                    _config_path)

    main_config_path = os.path.join(_root_dir, 'main.config')
    if os.path.exists(main_config_path):
        print('Main configuration file {} exists, skipping...\n'.format(
            main_config_path))
        exit(0)

    # if platform is linux then we just copy the default config file
    # as main.config
    if 'macOS' not in platform:
        print('System is {}, copying default config as main.config...'.format(
            platform))
        os.system('cp {} {}'.format(_config_path, main_config_path))
    else:
        print("System is {}, reconfiguring main.config...".format(platform))
        # set Backbone to the same as default
        default_config = configparser.ConfigParser()
        with open(_config_path, 'r') as f:
            default_config.read_file(f)
        
        for section in default_config.sections():
            cparser.add_section(section)
            for k, v in default_config[section].items():
                cparser.set(section, k, v)

        # configure MAGUS to use macOS compatible executables
        _macOS_dir = os.path.join(_root_dir, 'tools', 'macOS')

        binaries = os.popen('ls {}'.format(_macOS_dir)).read().split('\n')[:-1]
        for binary in binaries:
            path = os.path.join(_macOS_dir, binary)
            if 'FastTreeMP' in path:
                cparser.set('MAGUS', 'fasttreepath', path)
            else:
                cparser.set('MAGUS', '{}path'.format(binary), path)

        with open(main_config_path, 'w') as f:
            cparser.write(f)
    print('\n(Done) main.config written to {}'.format(main_config_path))

def main():
    parser = ArgumentParser()
    parser.add_argument('-c', '--configure', default=None,
            choices=['macOS', 'linux'], required=False,
            help='Specify which system the configuration file should resolve ' \
                 'to. This is an optional field and by default, the software ' \
                 'will resolve this automatically. Regardless, users can ' \
                 'specify either macOS or linux.')
    args = parser.parse_args()
    
    if args.configure:
        setup(args.configure)
    else:
        setup()

if __name__ == "__main__":
    main()
