"""
Module that handles the loading of the relevant data for the solver.
"""

import sys
import numpy as np
from enum import Enum, auto


def spinOSparser(pointerfile: str):
    """
    parses the parameter files which points to the different data files and guessfiles

    :param pointerfile: a text file containing the relative paths to the relevant datafiles, as well as a guess file
    :return: a dictionary containing the data, as well as a tag for the program to decide what actions to perform
    """
    print('Reading in your data...')
    wd = pointerfile
    token = wd[-1]
    remaining_tokens = len(pointerfile)
    while (token is not '/') and remaining_tokens > 0:
        remaining_tokens -= 1
        wd = wd[:-1]
        token = wd[-1]
    filetypes, datafiles = np.genfromtxt(pointerfile, dtype=None, encoding='utf-8', unpack=True)
    if filetypes.size == 1:
        filetypes, datafiles = [filetypes], [datafiles]
    if 'RVfile1' and 'RVfile2' and 'ASfile' in filetypes:
        tag = SpinOStag.SB2AS
    elif 'RVfile1' and 'RVfile2' in filetypes:
        tag = SpinOStag.SB2
    elif 'RVfile1' and 'ASfile' in filetypes:
        tag = SpinOStag.SB1AS
    elif 'RVfile1' in filetypes:
        tag = SpinOStag.SB1
    elif 'ASfile' in filetypes:
        tag = SpinOStag.AS
    else:
        tag = SpinOStag.PLOTONLY
    data_dict = dict()
    if 'guessfile' not in filetypes:
        print('no guesses supplied, cannot do a Levenberg–Marquardt minimization without an initial guess nor plot'
              ' something; stopping')
        exit()
    for i in range(len(filetypes)):
        if filetypes[i] == 'RVfile1':
            data_dict['RV1'] = np.loadtxt(wd + datafiles[i])
        elif filetypes[i] == 'RVfile2':
            data_dict['RV2'] = np.loadtxt(wd + datafiles[i])
        elif filetypes[i] == 'ASfile':
            data_dict['AS'] = np.loadtxt(wd + datafiles[i])
        elif filetypes[i] == 'guessfile':
            guesses = np.genfromtxt(wd + datafiles[i], dtype=None, filling_values=np.nan, encoding='utf-8')
            data_dict['guesses'] = dict()
            data_dict['fix_flags'] = dict()
            for guess in guesses:
                if guess[0] == 'i' or guess[0] == 'omega' or guess[0] == 'Omega':
                    data_dict['guesses'][guess[0]] = guess[1] * np.pi / 180
                else:
                    data_dict['guesses'][guess[0]] = guess[1]
                data_dict['fix_flags'][guess[0]] = guess[2]
        else:
            print('did not understood file pointer on line {}, only \'RVfile1\', \'RVfile2\', \'ASfile\' and '
                  '\'guessfile\' are supported'.format(i + 1))
    print('Data read complete!')
    return data_dict, tag


class SpinOStag(Enum):
    """
    Enumerative class containing tags for the program to decide what actions to undertake in the analysis
    """
    SB1 = auto()
    SB2 = auto()
    SB1AS = auto()
    SB2AS = auto()
    AS = auto()
    PLOTONLY = auto()
