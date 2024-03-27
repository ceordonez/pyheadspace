import os

import yaml
import pandas as pd


def read_data(cfg):
    path = cfg['input_path']
    filename = os.path.join(path, cfg['filename'])
    data = pd.read_excel(filename)
    if 'Remarks' in data.columns:
        del data['Remarks']
    return data

def read_config(filename):
    with open(filename, 'r') as file:
        conf_file = yaml.safe_load(file)
    return conf_file
