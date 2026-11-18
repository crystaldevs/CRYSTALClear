#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""


def _settings_path():
    """
    The file the executable paths are kept in, in the user's home directory,
    so that they outlast the session and are shared by every installation.

    Returns:
        path (str): ``~/.config/CRYSTALClear/settings.json``, expanded.
    """
    import os

    return os.path.join(os.path.expanduser('~'), '.config', 'CRYSTALClear',
                        'settings.json')


def _read_settings():
    """
    Read the saved settings.

    Returns:
        settings (dict): The saved settings, empty if there are none.
    """
    import json

    try:
        with open(_settings_path(), 'r') as file:
            settings = json.load(file)
    except (OSError, ValueError):
        return {}
    return settings if isinstance(settings, dict) else {}


def _save_setting(key, value):
    """
    Save one setting, keeping the others.

    Args:
        key (str): Name of the setting.
        value (str): Its value.
    """
    import json
    import os

    settings = _read_settings()
    settings[key] = value
    os.makedirs(os.path.dirname(_settings_path()), exist_ok=True)
    with open(_settings_path(), 'w') as file:
        json.dump(settings, file, indent=4)


def set_runcry_path(path):
    """
    Set the path for the Runcry executable.

    The path is saved in ``~/.config/CRYSTALClear/settings.json`` and used by
    every later call of :func:`runcry`. An empty path leaves it unchanged.

    Args:
        path (str): The path to the Runcry executable.
    """
    if len(path) > 0:
        _save_setting('runcry_path', path)


def set_runprop_path(path):
    """
    Set the path for the Runprop executable.

    The path is saved in ``~/.config/CRYSTALClear/settings.json`` and used by
    every later call of :func:`runprop`. An empty path leaves it unchanged.

    Args:
        path (str): The path to the Runprop executable.
    """
    if len(path) > 0:
        _save_setting('runprop_path', path)


def runcry(file_name, guessp=None):
    """
    Run Runcry calculation.

    Args:
        file_name (str): The name of the file to run the calculation.
        guessp (str, optional): The guessp parameter. Default is None.

    Returns:
        message (str | None): The result of the calculation or an error message.
    """
    runcry_path = _read_settings().get('runcry_path')
    if runcry_path is None:
        return 'Please set the runcry path before calling this function'

    import re
    import subprocess
    import sys

    converged = False
    index = 0

    # file_name = file_name.split('.')[0]

    if converged is False:
        if index > 3:
            return None
        else:
            if guessp is not None:
                run_calc = runcry_path + ' ' + file_name + ' ' + guessp
            else:
                run_calc = runcry_path + ' ' + file_name
            process = subprocess.Popen(
                run_calc.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            output, error = process.communicate()
            index += 1
            try:
                file = open('%s.out' % file_name, 'r')
                data = file.readlines()
                file.close()
            except:
                print('EXITING: a .out file needs to be specified')
                sys.exit(1)
            for i, line in enumerate(data[::-1]):
                if re.match(r'^ EEEEEEEEEE TERMINATION', line) is not None:
                    converged = True
                    return '%s.out calculation successfully completed' % file_name


def runprop(prop_name, wf_file):
    """
    Run Runprop calculation.

    Args:
        prop_name (str): The name of the property to calculate.
        wf_file (str): The name of the wavefunction file.

    Returns:
        message (str | None): The result of the calculation or an error message.
    """
    runprop_path = _read_settings().get('runprop_path')
    if runprop_path is None:
        return ('Please set the runprop path before calling it')

    import re
    import subprocess
    import sys

    wf_file = wf_file.split('.')[0]
    prop_name = prop_name.split('.')[0]
    converged = False

    index = 0
    while index < 3:
        run_calc = runprop_path + ' ' + prop_name + ' ' + wf_file
        process = subprocess.Popen(
            run_calc.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output, error = process.communicate()
        index = index + 1
        try:
            file = open(prop_name+'.outp', 'r')
            data = file.readlines()
            file.close()
        except:
            print('EXITING: a .out file needs to be specified')
            sys.exit(1)
        for i, line in enumerate(data[::-1]):
            if re.match(r'^ EEEEEEEEEE TERMINATION', line) != None:
                converged = True
                return '%s.outp calculation successfully completed' % prop_name

    return None
