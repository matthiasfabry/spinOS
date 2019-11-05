"""
Module that handles the loading of the relevant data for the solver.
"""
import numpy as np


def spinOSparser(pointerfile: str):
    """
    parses the parameter file which points to the different data files and guessfiles

    :param pointerfile: a text file containing the relative paths to the relevant datafiles
    :return: wd: the working directory
            plotonly: a bool letting the program know that no observational data is supplied
            filetypes: list of datafile types
            filenames: the relative path to the files
    """
    # Determine working directory
    wd = pointerfile
    token = wd[-1]
    remaining_tokens = len(pointerfile)
    while (token is not '/') and remaining_tokens > 0:
        remaining_tokens -= 1
        wd = wd[:-1]
        token = wd[-1]
    # parse the pointer file
    filetypes, filenames = np.genfromtxt(pointerfile, dtype=None, encoding='utf-8', unpack=True)
    # determine whether guesses were supplied
    if 'guessfile' not in filetypes:
        exit('no guesses supplied, cannot do a Levenberg–Marquardt minimization without an initial guess nor plot'
             ' something; stopping')
        return
    guessfile = filenames[filetypes == 'guessfile'][0]
    # determine whether any data is supplied
    if 'RVfile1' or 'RVfile2' or 'ASfile' in filetypes:
        plotonly = False
        # load data
        data_dict = data_loader(wd, filetypes, filenames)
        return wd, plotonly, guess_loader(wd, guessfile), data_dict
    else:
        plotonly = True
        return wd, plotonly, guess_loader(wd, guessfile)


def guess_loader(wd: str, guessfile: str) -> dict:
    # parse guesses
    guesses = np.genfromtxt(wd + guessfile, dtype=None, filling_values=np.nan, encoding='utf-8')
    guessdict = dict()
    guessdict['guesses'] = dict()
    guessdict['varying'] = dict()
    print('Reading guesses...')
    for guess in guesses:
        # convert degrees to radians
        if guess[0] == 'i' or guess[0] == 'omega' or guess[0] == 'Omega':
            guessdict['guesses'][guess[0]] = guess[1] * np.pi / 180
        else:
            guessdict['guesses'][guess[0]] = guess[1]
        # set flags whether to vary a parameters
        guessdict['varying'][guess[0]] = guess[2]
    print('Guess reading complete!\n')
    return guessdict


def data_loader(wd: str, filetypes: list, filenames: list) -> dict:
    data_dict = dict()
    print('Reading data...')
    for i in range(len(filetypes)):
        if filetypes[i] == 'RV1file':
            data_dict['RV1'] = np.loadtxt(wd + filenames[i])
        elif filetypes[i] == 'RV2file':
            data_dict['RV2'] = np.loadtxt(wd + filenames[i])
        elif filetypes[i] == 'ASfile':
            data_dict['AS'] = np.loadtxt(wd + filenames[i])
    print('Data reading complete!\n')
    return data_dict
