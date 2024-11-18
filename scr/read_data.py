import os

import pandas as pd
import yaml


def read_data(cfg):
    """Read input data.

    Parameters
    ----------
    cfg : dict
        Configuration information

    Return
    ------
    data : pandas.dataFrame

    """
    path = cfg["input_path"]
    filename = os.path.join(path, cfg["filename"])
    data = pd.read_excel(filename)
    if "Remarks" in data.columns:
        del data["Remarks"]
    return data


def read_config(filename):
    """Read configuration file.

    Parameters
    ----------
    filename : string
        Path to configuration file.

    Return
    -----
    conf_file : dict
        Configuration information.
    """
    with open(filename, "r") as file:
        conf_file = yaml.safe_load(file)
    return conf_file


def read_constant(filename):
    """Read constant file.

    Parameters
    ----------
    filename : string
        Path to constant file.

    Return
    -----
    conf_file : dict
        Constant data.
    """
    with open(filename, "r") as file:
        conf_file = yaml.safe_load(file)
    return conf_file
