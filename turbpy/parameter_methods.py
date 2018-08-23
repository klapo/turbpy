import yaml
import os


def yaml_dict(yaml_path):
    '''
    Reads .yml parameterization files (specifically for the default values).
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


def get_params(param_dict):
    '''
    Reads the user provided parameter dictionary variable and compares it to
    the available methods and determines which parameter values are necessary.
    Missing values are filled in with the defaults specified in
    default_params.yml
    '''

    # Get the default parameter/parameterizations.
    defaults = yaml_dict('default_params.yml')
    # Initialize the output parameter dictionary
    out_param_dict = defaults

    # Determine if a stability method was provided:
    try:
        stab_method = param_dict['stability_method']
    except KeyError:
        stab_method = defaults['stability_method']

    # Determine if the stability method name matches the names available.
    if stab_method not in defaults['available_stab_methods']:
        print('Unrecognized stability choice: '
              + stab_method
              + '\nValid stability options: '
              + ', '.join(defaults['available_stab_methods']))

    # Determine any tunable parameters.
    methods_need_params = ['standard', 'louis', 'mahrt']
    if stab_method in methods_need_params:
        try:
            stab_param = param_dict['stability_params'][stab_method]
        except KeyError:
            stab_param = defaults['stability_params'][stab_method]

    # Deal with Monin-Obukhov method implementation seperately.
    if stab_method == 'monin_obukhov':
        # Determine which gradient function to use.
        try:
            gradient_function = param_dict['monin_obukhov']['gradient_function']
        except KeyError:
            gradient_function = defaults['monin_obukhov']['gradient_function']

        # Determine if the gradient function name is correct
        if gradient_function not in defaults['monin_obukhov']['available_gradient_funcs']:
            print('Unrecognized gradient function choice for Monin-Obukhov: '
                  + gradient_function
                  + '\nValid stability options: '
                  + ', '.join(defaults['monin_obukhov']['available_gradient_funcs']))

        # Catch the one case with a tunable parameter.
        if gradient_function == 'Webb':
            try:
                stab_param = param_dict['stability_params']['Webb']
            except KeyError:
                stab_param = defaults['stability_params']['Webb']
        else:
            stab_param = None

        # Find out what capping should be used.
        try:
            capping = param_dict['monin_obukhov']['capping']
        except KeyError:
            capping = defaults['monin_obukhov']['capping']

        # Determine if the capping option is correct
        if capping not in defaults['monin_obukhov']['capping']:
            print('Unrecognized capping choice for Monin-Obukhov: '
                  + capping
                  + '\nValid stability options: '
                  + ', '.join(defaults['monin_obukhov']['available_capping']))

        # Assign MO values to the parameter dictionary to return to turbpy
        out_param_dict['monin_obukhov']['gradient_function'] = gradient_function
        out_param_dict['monin_obukhov']['capping'] = capping

    # Assign last values for returning.
    out_param_dict['stability_methods'] = stab_method
    out_param_dict['stability_params'][stab_method] = stab_param

    return(out_param_dict)
