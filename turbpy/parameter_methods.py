import yaml
import os


def yaml_dict(yaml_path):
    '''
    Reads .yml parameterization files.
    INPUT:
        yaml_path - path to the .yml file to read
    OUTPUT:
        param_dict - python dictionary of parameter values and methods.
    '''
    if not os.path.exists(yaml_path):
        raise IOError('Could not find yml file at ' + yaml_path)

    with open(yaml_path, 'r') as ymlfile:
        param_dict = yaml.load(ymlfile)
    return param_dict
