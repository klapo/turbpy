import yaml
import os
import pkg_resources


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
    defaults_yml = pkg_resources.resource_filename('turbpy',
                                                   'default_params.yml')
    defaults = yaml_dict(defaults_yml)
    # Initialize the output parameter dictionary
    out_param_dict = defaults

    # Determine if a stability method was provided:
    try:
        stab_method = param_dict['stability_method']
    except KeyError:
        stab_method = defaults['stability_method']

    # Determine if the stability method name matches the names available.
    if stab_method not in defaults['available_stab_methods']:
        raise ValueError('Unrecognized stability choice: ' + stab_method
                         + '\nValid stability options: '
                         + ', '.join(defaults['available_stab_methods']))

    # Determine any tunable parameters.
    methods_need_params = ['standard', 'louis', 'mahrt']
    if stab_method in methods_need_params:
        try:
            stab_param = param_dict['stability_params'][stab_method]
        except KeyError:
            stab_param = defaults['stability_params'][stab_method]

    # Find out what capping should be used.
    try:
        capping = param_dict['capping']
    except KeyError:
        capping = defaults['capping']
    # Determine if the capping option is correct
    if capping not in defaults['available_capping']:
        raise ValueError('Unrecognized capping choice: ' + capping
                         + '\nValid stability options: '
                         + ', '.join(defaults['available_capping']))

    if capping == 'louis_Ri_capping' and not stab_method == 'louis':
        raise ValueError('Louis capping (MJ98) can only be implemented \
                          with the Louis stability scheme.')

    # Deal with Monin-Obukhov method implementation seperately.
    if stab_method == 'monin_obukhov':
        # Determine which gradient function to use.
        try:
            gradient_function = param_dict['monin_obukhov']['gradient_function']
        except KeyError:
            gradient_function = defaults['monin_obukhov']['gradient_function']

        # Determine if the gradient function name is correct
        if gradient_function not in defaults['monin_obukhov']['available_gradient_funcs']:
            raise ValueError('Unrecognized gradient function choice\
                             for Monin-Obukhov: ' + gradient_function
                             + '\nValid stability options: '
                             + ', '.join(defaults['monin_obukhov']['available_gradient_funcs']))

        # Catch the one case with a tunable parameter.
        if gradient_function == 'webb':
            try:
                stab_param = param_dict['stability_params']['webb']
            except KeyError:
                stab_param = defaults['stability_params']['webb']
        else:
            stab_param = None

        # Determine if the capping option is correct
        if capping not in defaults['monin_obukhov']['available_capping']:
            raise ValueError('Unrecognized capping choice for Monin-Obukhov: '
                             + capping
                             + '\nValid stability options: '
                             + ', '.join(defaults['monin_obukhov']['available_capping']))

        # Determine the roughness length parameterizations
        try:
            roughness_function = param_dict['monin_obukhov']['roughness_function']
        except KeyError:
            roughness_function = defaults['monin_obukhov']['roughness_function']

        # Determine if the roughness function is correct
        if roughness_function not in defaults['monin_obukhov']['available_roughness_funcs']:
            raise ValueError('Unrecognized roughness length option for Monin-Obukhov: '
                             + roughness_function
                             + '\nValid roughness length options: '
                             + ', '.join(defaults['monin_obukhov']['available_roughness_funcs']))

        # Determine how the conductance is calculated
        try:
            conductance_approx = param_dict['monin_obukhov']['conductance_approx']
        except KeyError:
            conductance_approx = defaults['monin_obukhov']['conductance_approx']

        # Determine if the roughness function is correct
        if conductance_approx not in defaults['monin_obukhov']['available_cond_funcs']:
            raise ValueError('Unrecognized conductance option for Monin-Obukhov: '
                             + conductance_approx
                             + '\nValid conductance options: '
                             + ', '.join(defaults['monin_obukhov']['available_cond_funcs']))

        # Assign MO values to the parameter dictionary to return to turbpy
        out_param_dict['monin_obukhov']['gradient_function'] = gradient_function
        out_param_dict['monin_obukhov']['capping'] = capping
        out_param_dict['monin_obukhov']['roughness_function'] = roughness_function
        out_param_dict['monin_obukhov']['conductance_approx'] = conductance_approx

    # Assign last values for returning.
    out_param_dict['stability_method'] = stab_method
    out_param_dict['stability_params'][stab_method] = stab_param
    out_param_dict['capping'] = capping

    return(out_param_dict)
