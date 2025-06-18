import configparser
import os

#config = configparser.ConfigParser(inline_comment_prefixes="#")
#config.read('init_simulation.cfg')

path_current_directory = os.path.dirname(__file__)
path_config_file = os.path.join(path_current_directory, 'init_simulation.cfg')
config = configparser.ConfigParser()
config.read(path_config_file)


PATH = path_current_directory