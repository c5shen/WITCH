import os, sys, shutil
try:
    import configparser
except ImportError:
    import ConfigParser as configparser
from argparse import ArgumentParser, Namespace
from platform import platform

def find_main_config(homepath):
    with open(homepath, 'r') as f:
        _root_dir = f.read().strip()
        main_config_path = os.path.join(_root_dir, 'main.config')
        return _root_dir, main_config_path

'''
first time run, need user to initialize the main.config
if it is not installed through github (i.e., python setup.py config)
will be needed if installed through pip/pypi
'''
def init_config_file(homepath, prioritize_user_software=True):
    # read from sys.argv to find if "-y" or "--bypass-setup" exists
    args = sys.argv[1:]
    bypass_setup = False
    if '-y' in args or '--bypass-setup' in args:
        bypass_setup = True

    # initialize a home.path that points to local user main.config
    # if it exists then pass on
    if os.path.exists(homepath):
        # detecting old home.path file based on creation time
        if os.stat(homepath).st_mtime >= os.stat(__file__).st_mtime:
            return find_main_config(homepath)
        else:
            print('Found old home.path and regenerating...')
            os.remove(homepath)
    else:
        print('Cannot find home.path: {}'.format(homepath))

    # install to user local directory
    # bypassing the setup step to directly use the default path
    _root_dir = ''
    if not bypass_setup:
        _root_dir = input('Create main.config file at [default: ~/.witch_msa/]: ')

    if _root_dir == '':
        _root_dir = os.path.expanduser('~/.witch_msa')
    else:
        _root_dir = os.path.abspath(_root_dir)
    main_config_path = os.path.join(_root_dir, 'main.config')
    print('Initializing main configuration file: {}...'.format(main_config_path))

    # write to local for installation to system
    # will read in during runs to find the main.config file
    if not os.path.isdir(_root_dir):
        os.mkdir(_root_dir)
    with open(homepath, 'w') as f:
        f.write(_root_dir)

    # create main.config file at configfile using default.config
    _config_path = os.path.join(os.path.dirname(__file__), 'default.config')
    config_defaults = []
    cparser = configparser.ConfigParser()
    cparser.optionxform = str
    assert os.path.exists('{}'.format(_config_path)), \
            "default config file {} missing! Please redownload from Github\n".format(
                    _config_path)

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
    tools_dir = os.path.join(os.path.dirname(__file__), 'tools')
    set_sections = ['Basic', 'MAGUS']

    # copy magus directory to tools/
    magus_dir = os.path.join(tools_dir, 'magus')
    cparser.set('Basic', 'magus_path', magus_dir + '/magus.py')

    if 'macos' not in platform_name.lower():
        print('System is {}, using default config as main.config...'.format(
            platform_name))
        # use existing binaries from MAGUS subfolder (reduce redundancy of
        # duplicated binaries)
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
        binaries = os.listdir(tools_dir + '/macOS')
        for binary in binaries:
            path = os.path.join(tools_dir, 'macOS', binary)
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
    # DO NOT EXIT; can start running WITCH with any given commands now
    #exit(0)
    return _root_dir, main_config_path
