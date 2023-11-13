# setup.py for PyPI
import setuptools
from setuptools import find_packages
import os, sys, shutil
try:
    import configparser
except ImportError:
    import ConfigParser as configparser
from argparse import ArgumentParser, Namespace

_root_dir = os.path.dirname(os.path.realpath(__file__))
_config_path = os.path.join(_root_dir, 'default.config')

def get_version(path):
    with open(_root_dir + '/' + path, 'r') as f:
        lines = f.read().splitlines()
        for line in lines:
            if line.startswith('__version__'):
                delim = '"' if '"' in line else "'"
                return line.split(delim)[1]

classifiers = [
        'Development Status :: 4 - Beta',
        'Intended Audience :: Scientists',
        'Topic :: Bioinformatics',
        'License :: GPL-3.0',
        'Programming Language :: Python :: 3',
        ]
requirements = [
        'cython>=0.29', 'configparser>=5.0.0',
        'Dendropy>=4.4.0,<4.6.0', 'numpy>=1.15',
        'psutil>=5.0', 'tqdm>=4.0.0',
        ]

#def setup(prioritize_user_software, platform=platform()): 
def main():
    from platform import platform
    prioritize_user_software = True
    platform = platform()

    config_defaults = []
    cparser = configparser.ConfigParser()
    cparser.optionxform = str
    assert os.path.exists('{}'.format(_config_path)), \
            "default config file {} missing! Please redownload from Github\n".format(
                    _config_path)

    main_config_path = os.path.join(_root_dir, 'main.config')
    if os.path.exists(main_config_path):
        print('Main configuration file {} exists...'.format(main_config_path))
        print('Overwriting existing configuration file...')
        #ans = input('Do you wish to regenerate the file? [yes(y)/no(n)]:')
        #if ans != 'yes' and ans != 'y':
        #    print('Exiting WITCH configuration setup...')
        #    exit(0)

    print('\n\n')
    # initialize main config file using default config file
    default_config = configparser.ConfigParser()
    with open(_config_path, 'r') as f:
        default_config.read_file(f)
    for section in default_config.sections():
        cparser.add_section(section)
        for k, v in default_config[section].items():
            cparser.set(section, k, v)

    # if platform is linux then we just copy the default config file
    # as main.config
    set_sections = ['Basic', 'MAGUS']
    if 'macOS' not in platform:
        print('System is {}, using default config as main.config...'.format(
            platform))
        # use existing binaries from MAGUS subfolder (reduce redundancy of
        # duplicated binaries)
        _magus_tools_dir = os.path.join(_root_dir, 'tools', 'magus', 'tools') 
        for _section in set_sections:
            # mafftpath
            cparser.set(_section, 'mafftpath',
                    os.path.join(_magus_tools_dir, 'mafft', 'mafft'))
            # mclpath
            cparser.set(_section, 'mclpath',
                    os.path.join(_magus_tools_dir, 'mcl', 'bin', 'mcl'))
            # fasttreepath
            cparser.set(_section, 'fasttreepath',
                    os.path.join(_magus_tools_dir, 'fasttree', 'FastTreeMP'))
            # hmmer packages
            for hmmer_pkg in ['hmmsearch', 'hmmalign', 'hmmbuild']:
                cparser.set(_section, '{}path'.format(hmmer_pkg),
                        os.path.join(_magus_tools_dir, 'hmmer', hmmer_pkg))
    else:
        if 'x86' not in platform:
            print('Warning: system is not using x86 architecture.',
                    'Some softwares such as FastTreeMP need to be',
                    'self-provided. See {} [Basic] '.format(_config_path),
                    'section for more information.')
        print("System is {}, reconfiguring main.config...".format(platform))

        # configure MAGUS to use macOS compatible executables
        _macOS_dir = os.path.join(_root_dir, 'tools', 'macOS')

        binaries = os.popen('ls {}'.format(_macOS_dir)).read().split('\n')[:-1]
        for binary in binaries:
            path = os.path.join(_macOS_dir, binary)
            for _section in set_sections:
                if 'FastTreeMP' in path:
                    cparser.set(_section, 'fasttreepath', path)
                else:
                    cparser.set(_section, '{}path'.format(binary), path)

    # binaries from the user's environment will be used in priority
    # if they exist
    if prioritize_user_software:
        print('Detecting existing software from the user\'s environment...')
        software = ['mafft', 'mcl', 
                'hmmsearch', 'hmmalign', 'hmmbuild', 'FastTreeMP']
        for soft in software:
            print('\t{}: {}'.format(soft, shutil.which(soft)))
            for _section in set_sections:
                if shutil.which(soft):
                    if soft != 'FastTreeMP':
                        cparser.set(_section, '{}path'.format(soft),
                                shutil.which(soft))
                    else:
                        cparser.set(_section, 'fasttreepath',
                                shutil.which(soft))

    with open(main_config_path, 'w') as f:
        cparser.write(f)
    print('\n(Done) main.config written to {}'.format(main_config_path))
    print('If you would like to make manual changes, please directly edit {}'.format(
        main_config_path))

    setuptools.setup(name="witch-msa",
            packages=find_packages(),
            version=get_version('witch.py'),
            description="WITCH - A Multiple Sequence Alignment Tool",
            url="https://github.com/c5shen/WITCH",
            license="GPL-3.0",
            classifiers=classifiers,
            install_requires=requirements,
            python_requires='>=3.7')
        #long_description='text/x-rst',

#def main():
#    parser = ArgumentParser()
#    parser.add_argument('-c', '--configure', default=None,
#            choices=['macOS', 'linux'], required=False,
#            help='Specify which system the configuration file should resolve ' \
#                 'to. This is an optional field and by default, the software ' \
#                 'will resolve this automatically. Regardless, users can ' \
#                 'specify either macOS (x86) or linux.')
#    parser.add_argument('-p', '--prioritize-user-software', default="true",
#            choices=['true', 'false'], required=False,
#            help='Whether to prioritize using software in the user\'s ' \
#                 'environment. default: True')
#    args = parser.parse_args()
#    
#    if args.configure:
#        setup(args.prioritize_user_software == "true", args.configure) 
#    else:
#        setup(args.prioritize_user_software == "true")
#
if __name__ == "__main__":
    main()
