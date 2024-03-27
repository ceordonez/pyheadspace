import os
import pandas as pd

def write_data(cfg, data):
    """Write data to excel file.

    Parameters
    ----------
    cfg : TODO
    data : TODO

    Returns
    -------
    TODO

    """
    filename = os.path.join(cfg['output_path'], cfg['filename_out'])
    writer = pd.ExcelWriter(filename, datetime_format='mmm d yyyy hh:mm:ss',
                            engine='xlsxwriter')
    data.to_excel(writer, sheet_name = 'Results', index=False)#, float_format='%.3f')
    workbook = writer.book
    worksheet = writer.sheets['Results']
    dateformat = workbook.add_format({'num_format': 'yyyy/mm/dd hh:mm'})
    floatformat = workbook.add_format({'num_format': '0.000'})
    worksheet.set_column('B:B', None, dateformat)
    worksheet.set_column('F:X', None, floatformat)
    writer._save()
