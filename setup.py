#!/usr/bin/env python

# setup.py for PyPI
import setuptools
from setuptools import Command, find_packages
import os, sys, shutil
try:
    import configparser
except ImportError:
    import ConfigParser as configparser
from argparse import ArgumentParser, Namespace
from platform import platform

_root_dir = os.path.dirname(__file__)
_config_path = os.path.join(_root_dir, 'witch_msa', 'default.config')

def get_version(path):
    with open(_root_dir + '/' + path, 'r') as f:
        lines = f.read().splitlines()
        for line in lines:
            if line.startswith('__version__'):
                delim = '"' if '"' in line else "'"
                return line.split(delim)[1]

def get_tools_dir(where=None):
    platform_name = platform()
    if where is None:
        if 'macos' in platform_name.lower():
            where = os.path.join("macOS")
        else:
            where = os.path.join("magus", "tools")
    path = os.path.join(os.getcwd(), "witch_msa", "tools", where)
    if not os.path.exists(path):
        raise OSError("Tool is not available on {}!".format(platform_name))
    return path

class ConfigWITCH(Command):
    description = "Configure WITCH for the current user and create main.config."
    user_options = [
            ('contained', 'c',
                ('Whether WITCH should be install to ~/.witch_msa '
                    'or in a self-contained manner.'))]

    def initopts(self):
        self.contained = None
        self.configfile = None
        self.basepath = None

    def initpath(self, name='main.config'):
        # install to WITCH directory
        if self.contained:
            self.configfile = os.path.expanduser(
                    os.path.abspath(os.path.join('.witch_msa', name)))
        else:
            self.configfile = os.path.expanduser('~/.witch_msa/{}'.format(name))
        self.basepath = os.path.dirname(self.configfile)
        # write to local for installation to system
        # will read in during runs to find the main.config file
        with open('witch_msa/home.path', 'w') as f:
            f.write(self.basepath)
            f.close()
        #return basepath, configfile

    def initialize_options(self):
        self.initopts()

    def finalize_options(self):
        self.initpath('main.config')
        print('\nCreating WITCH main.config at {} and tools at {}'.format(
            self.configfile, self.basepath))

    def get_tools_dest(self):
        return os.path.join(self.basepath, "tools")

    def copy_tool_to_lib(self, tool, _from, bits=True):
        _to = os.path.join(self.get_tools_dest(), tool)
        if not os.path.isdir(_from): 
            if os.path.exists(_to):
                os.rmove(_to)
            shutil.copy2(_from, _to)
        else:
            if os.path.exists(_to):
                shutil.rmtree(_to)
            shutil.copytree(_from, _to, False, None)
        return _to

    def run(self):
        '''
            prioritize software that exists in user's environment
        '''
        prioritize_user_software = True

        if not os.path.isdir(self.basepath):
            # avoid making nested directories just in case
            os.mkdir(self.basepath)
        #if not os.path.isdir(self.get_tools_dest()):
        #    os.mkdir(self.get_tools_dest())

        # create main.config file at self.configfile
        config_defaults = []
        cparser = configparser.ConfigParser()
        cparser.optionxform = str
        assert os.path.exists('{}'.format(_config_path)), \
                "default config file {} missing! Please redownload from Github\n".format(
                        _config_path)

        main_config_path = self.configfile
        if os.path.exists(main_config_path):
            print('Main configuration file {} exists...'.format(main_config_path))
            print('Overwriting existing configuration file...')

        print('\n')
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
        platform_name = platform()
        tools_dir = get_tools_dir()
        set_sections = ['Basic', 'MAGUS']

        # copy magus directory to tools/
        magus_dir = os.path.join(os.getcwd(), "witch_msa", "tools", "magus")
        #copied_magus_dir = self.copy_tool_to_lib(
        #        'magus', os.path.join(os.getcwd(), "tools", "magus"))
        cparser.set('Basic', 'magus_path', magus_dir + '/magus.py')

        if 'macos' not in platform_name.lower():
            print('System is {}, using default config as main.config...'.format(
                platform_name))
            # use existing binaries from MAGUS subfolder (reduce redundancy of
            # duplicated binaries)
            #_magus_tools_dir = os.path.join(_root_dir, 'tools', 'magus', 'tools') 
            for _section in set_sections:
                # mafftpath
                cparser.set(_section, 'mafftpath',
                        os.path.join(magus_dir, 'tools', 'mafft', 'mafft'))
                # mclpath
                cparser.set(_section, 'mclpath',
                        os.path.join(magus_dir, 'tools', 'mcl', 'bin', 'mcl'))
                # fasttreepath
                cparser.set(_section, 'fasttreepath',
                        os.path.join(magus_dir, 'tools', 'fasttree', 'FastTreeMP'))
                # hmmer packages
                for hmmer_pkg in ['hmmsearch', 'hmmalign', 'hmmbuild']:
                    cparser.set(_section, '{}path'.format(hmmer_pkg),
                            os.path.join(magus_dir, 'tools', 'hmmer', hmmer_pkg))
        else:
            if 'x86' not in platform_name:
                print('Warning: system is not using x86 architecture.',
                        'Some softwares such as FastTreeMP need to be',
                        'self-provided. See {} [Basic] '.format(_config_path),
                        'section for more information.')
            print("System is {}, reconfiguring main.config...".format(platform_name))

            # configure MAGUS to use macOS compatible executables
            #_macOS_dir = os.path.join(_root_dir, 'tools', 'macOS')
            #binaries = os.popen('ls {}'.format(_macOS_dir)).read().split('\n')[:-1]
            binaries = os.listdir(tools_dir)
            for binary in binaries:
                path = os.path.join(tools_dir, binary)
                #path = os.path.join(_macOS_dir, binary)
                for _section in set_sections:
                    if 'FastTreeMP' in path:
                        cparser.set(_section, 'fasttreepath', path)
                                #self.copy_tool_to_lib('FastTreeMP', path))
                    else:
                        cparser.set(_section, '{}path'.format(binary), path)
                                #self.copy_tool_to_lib(binary, path))

        # binaries from the user's environment will be used in priority
        # if they exist
        if prioritize_user_software:
            print('Detecting existing software from the user\'s environment...')
            software = ['mafft', 'mcl', 
                    'hmmsearch', 'hmmalign', 'hmmbuild', 'FastTreeMP']
            print('\tDetected:\n')
            for soft in software:
                if shutil.which(soft):
                    print('\t{}: {}'.format(soft, shutil.which(soft)))
                    for _section in set_sections:
                        if soft == 'FastTreeMP':
                            cparser.set(_section, 'fasttreepath',
                                    shutil.which(soft))
                        elif soft == 'magus':
                            cparser.set('Basic', 'magus_path',
                                    shutil.which(soft))
                        else:
                            cparser.set(_section, '{}path'.format(soft),
                                    shutil.which(soft))

        with open(main_config_path, 'w') as f:
            cparser.write(f)
        print('\n(Done) main.config written to {}'.format(main_config_path))
        print('If you would like to make manual changes, please directly edit {}'.format(
            main_config_path))

#def create_main_config_path(prioritize_user_software=True):
#    platform_name = platform()
#    #prioritize_user_software = True
#
#    config_defaults = []
#    cparser = configparser.ConfigParser()
#    cparser.optionxform = str
#    assert os.path.exists('{}'.format(_config_path)), \
#            "default config file {} missing! Please redownload from Github\n".format(
#                    _config_path)
#
#    # initialize main.config location and record in home.path
#    basepath, configfile = initpath()
#    if not os.path.isdir(basepath):
#        os.makedirs(basepath)
#
#    #main_config_path = os.path.join(_root_dir, 'main.config')
#    main_config_path = configfile
#    if os.path.exists(main_config_path):
#        print('Main configuration file {} exists...'.format(main_config_path))
#        print('Overwriting existing configuration file...')
#        #ans = input('Do you wish to regenerate the file? [yes(y)/no(n)]:')
#        #if ans != 'yes' and ans != 'y':
#        #    print('Exiting WITCH configuration setup...')
#        #    exit(0)
#
#    print('\n\n')
#    # initialize main config file using default config file
#    default_config = configparser.ConfigParser()
#    with open(_config_path, 'r') as f:
#        default_config.read_file(f)
#    for section in default_config.sections():
#        cparser.add_section(section)
#        for k, v in default_config[section].items():
#            cparser.set(section, k, v)
#
#    # if platform is linux then we just copy the default config file
#    # as main.config
#    set_sections = ['Basic', 'MAGUS']
#    if 'macOS' not in platform_name:
#        print('System is {}, using default config as main.config...'.format(
#            platform_name))
#        # use existing binaries from MAGUS subfolder (reduce redundancy of
#        # duplicated binaries)
#        _magus_tools_dir = os.path.join(_root_dir, 'tools', 'magus', 'tools') 
#        for _section in set_sections:
#            # mafftpath
#            cparser.set(_section, 'mafftpath',
#                    os.path.join(_magus_tools_dir, 'mafft', 'mafft'))
#            # mclpath
#            cparser.set(_section, 'mclpath',
#                    os.path.join(_magus_tools_dir, 'mcl', 'bin', 'mcl'))
#            # fasttreepath
#            cparser.set(_section, 'fasttreepath',
#                    os.path.join(_magus_tools_dir, 'fasttree', 'FastTreeMP'))
#            # hmmer packages
#            for hmmer_pkg in ['hmmsearch', 'hmmalign', 'hmmbuild']:
#                cparser.set(_section, '{}path'.format(hmmer_pkg),
#                        os.path.join(_magus_tools_dir, 'hmmer', hmmer_pkg))
#    else:
#        if 'x86' not in platform_name:
#            print('Warning: system is not using x86 architecture.',
#                    'Some softwares such as FastTreeMP need to be',
#                    'self-provided. See {} [Basic] '.format(_config_path),
#                    'section for more information.')
#        print("System is {}, reconfiguring main.config...".format(platform_name))
#
#        # configure MAGUS to use macOS compatible executables
#        _macOS_dir = os.path.join(_root_dir, 'tools', 'macOS')
#
#        binaries = os.popen('ls {}'.format(_macOS_dir)).read().split('\n')[:-1]
#        for binary in binaries:
#            path = os.path.join(_macOS_dir, binary)
#            for _section in set_sections:
#                if 'FastTreeMP' in path:
#                    cparser.set(_section, 'fasttreepath', path)
#                else:
#                    cparser.set(_section, '{}path'.format(binary), path)
#
#    # binaries from the user's environment will be used in priority
#    # if they exist
#    if prioritize_user_software:
#        print('Detecting existing software from the user\'s environment...')
#        software = ['mafft', 'mcl', 
#                'hmmsearch', 'hmmalign', 'hmmbuild', 'FastTreeMP']
#        for soft in software:
#            print('\t{}: {}'.format(soft, shutil.which(soft)))
#            for _section in set_sections:
#                if shutil.which(soft):
#                    if soft != 'FastTreeMP':
#                        cparser.set(_section, '{}path'.format(soft),
#                                shutil.which(soft))
#                    else:
#                        cparser.set(_section, 'fasttreepath',
#                                shutil.which(soft))
#
#    with open(main_config_path, 'w') as f:
#        cparser.write(f)
#    print('\n(Done) main.config written to {}'.format(main_config_path))
#    print('If you would like to make manual changes, please directly edit {}'.format(
#        main_config_path))

#def setup(prioritize_user_software, platform=platform()): 
#def main():
    #create_main_config_path()
    #b_recreate_main = False
    #if len(sys.argv) > 1:
    #    b_recreate_main = sys.argv[1] == 'main'

setuptools.setup(
        #name="witch-msa",
        packages=["witch_msa", "witch_msa/gcmm", "witch_msa/helpers"], #find_packages(),
        package_data={"witch_msa": ["witch_msa/tools"]},
        include_package_data=True,
        cmdclass={"config": ConfigWITCH},
        #version=get_version('witch.py'),
        #description="WITCH - A Multiple Sequence Alignment Tool",
        long_description=open('README.rst', 'r').read(),
        long_description_content_type='text/x-rst',
        #author='Chengze Shen',
        #author_email='chengze5@illinois.edu',
        #url="https://github.com/c5shen/WITCH",
        #license="GPL-3.0",
        #classifiers=classifiers,
        #install_requires=requirements,
        #python_requires='>=3.7',
        scripts=['bin/witch-msa', 'witch.py'],
        #data_files=[("", "home.path")]
        )

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
#if __name__ == "__main__":
#    main()
