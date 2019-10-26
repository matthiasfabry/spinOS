"""
Module that handles the loading of the relevant data for the solver.
"""

import sys
import numpy as np


def spinOSparser(pointerfile: str):
    """

    :param pointerfile:
    :return:
    """
    wd = pointerfile
    token = wd[-1]
    while token is not '/':
        wd = wd[:-1]
        token = wd[-1]

    readin = np.genfromtxt(pointerfile, dtype=None, encoding='utf-8')
    filetypes = readin[:, 0]
    datafiles = readin[:, 1]
    data_dict = dict()
    if 'guessfile' not in readin[:, 0]:
        print('no guesses supplied, cannot do a Levenberg–Marquardt algorithm without an initial guess; stopping')
        exit()
    for i in range(len(filetypes)):
        if filetypes[i] == 'RVfile1':
            data_dict['RV1'] = np.loadtxt(wd + datafiles[i])
        elif filetypes[i] == 'RVfile2':
            data_dict['RV2'] = np.loadtxt(wd + datafiles[i])
        elif filetypes[i] == 'ASfile':
            data_dict['AS'] = np.loadtxt(wd + datafiles[i])
        elif filetypes[i] == 'guessfile':
            data_dict['guesses'] = np.genfromtxt(wd + datafiles[i], dtype=None, filling_values=np.nan, encoding='utf-8')
        else:
            print('did not understood file pointer on line {},'
                  ' only \'RVfile1\' \'RVfile2\' \'ASfile\' and \'guessfile\' '
                  'are supported'.format(i+1))
    return data_dict


print(spinOSparser(sys.argv[1]))
