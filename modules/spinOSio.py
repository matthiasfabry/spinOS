"""
Module that handles the loading of the relevant data for the solver.
This module is developed using numpy 1.18.1

Author:
Matthias Fabry, Instituut voor Sterrekunde, KU Leuven, Belgium

Date:
21 Jan 2020
"""
import os
import numpy as np
from modules import constants as c


def spinOSparser(pointerfile: str, doseppaconversion: bool = True):
    """
    parses the parameter file which points to the different data files and guessfiles

    :param doseppaconversion: boolean indicating whether the data is in sep/pa format (True) so that it converts it to
                                east/north format
    :param pointerfile: a text file containing the relative paths to the relevant datafiles
    :return: a dictionary containing the guessed parameters of the system
             a dictionary containing the observational data (possibly empty)
    """

    # parse the pointer file
    wd = os.path.dirname(os.path.abspath(pointerfile))+'/'
    filetypes, filenames = np.genfromtxt(pointerfile, dtype=None, encoding='utf-8', unpack=True)

    # determine whether guesses were supplied
    if 'guessfile' not in filetypes:
        exit('no guesses supplied, cannot do a Levenberg–Marquardt minimization without an initial guess; stopping')
        return
    if filetypes.size > 1:
        guessfile = filenames[filetypes == 'guessfile'][0]
    else:
        guessfile = filenames
    # determine whether any data is supplied
    if 'RVfile1' or 'RVfile2' or 'ASfile' in filetypes:
        # load data
        data_dict = data_loader(wd, filetypes, filenames, doseppaconversion)
        return wd, guess_loader(wd, guessfile), data_dict
    else:
        return wd, guess_loader(wd, guessfile), dict()


def guess_loader(wd: str, guessfile: str) -> dict:
    """
    parses the guess file and determines values and flags for each guess
    :param wd: the working directory
    :param guessfile: pathname (relative to wd) pointing to the file containing guesses
    :return: dictionary containing the guesses and flags for each parameter
    """
    guesses = np.genfromtxt(wd + guessfile, dtype=None, filling_values=np.nan, encoding='utf-8')
    guessdict = dict()
    for guess in guesses:
        guessdict[guess[0]] = (guess[1], guess[2])
    return guessdict


def guess_saver(wd: str, guess_dict: dict) -> None:
    with open(wd+'guesses.txt', 'w') as guessfile:
        for param, guess in guess_dict.items():
            guessfile.write(param + ' {} {}\n'.format(guess[0], str(guess[1])))


def data_loader(wd: str, filetypes: list, filenames: list, doseppaconversion: bool = True) -> dict:
    data_dict = dict()
    for i in range(len(filetypes)):
        if filetypes[i] == 'RV1file':
            data = np.loadtxt(wd + filenames[i])
            data_dict['RV1'] = dict()
            data_dict['RV1']['hjds'] = data[:, 0]
            data_dict['RV1']['RVs'] = data[:, 1]
            data_dict['RV1']['errors'] = data[:, 2]
        elif filetypes[i] == 'RV2file':
            data = np.loadtxt(wd + filenames[i])
            data_dict['RV2'] = dict()
            data_dict['RV2']['hjds'] = data[:, 0]
            data_dict['RV2']['RVs'] = data[:, 1]
            data_dict['RV2']['errors'] = data[:, 2]
        elif filetypes[i] == 'ASfile':
            data = np.loadtxt(wd + filenames[i])
            data_dict['AS'] = dict()
            data_dict['AS']['hjds'] = data[:, 0]
            data_dict['AS']['easterrors'], data_dict['AS']['northerrors'] = \
                convert_error_ellipse(data[:, 3], data[:, 4], data[:, 5] * c.deg2rad)
            data_dict['AS']['majors'] = data[:, 3]
            data_dict['AS']['minors'] = data[:, 4]
            data_dict['AS']['pas'] = data[:, 5]
            if doseppaconversion:
                data_dict['AS']['easts'] = data[:, 1] * np.sin(data[:, 2] * c.deg2rad)
                data_dict['AS']['norths'] = data[:, 1] * np.cos(data[:, 2] * c.deg2rad)
            else:
                data_dict['AS']['easts'] = data[:, 1]
                data_dict['AS']['norths'] = data[:, 2]
    return data_dict


def convert_error_ellipse(major, minor, angle):
    """
    Converts error ellipses to actual east and north errors by a sampling the error ellipse in a monte-carlo way and
    then taking the variance in the east and north directions.
    :param major: length of the major axis of the error ellipse
    :param minor: length of the minor axis of the error ellipse
    :param angle: position angle east of north of the major axis
    :return: east and north error
    """
    num = 1000
    east_error = np.zeros(len(major))
    north_error = np.zeros(len(major))
    for i in range(len(east_error)):
        cosa = np.cos(angle[i])
        sina = np.sin(angle[i])
        temp_major = np.random.randn(num) * major[i]
        temp_minor = np.random.randn(num) * minor[i]
        rotated_temp = np.matmul(np.array([[cosa, sina], [-sina, cosa]]), [temp_major, temp_minor])
        east_error[i] = np.std(rotated_temp[0])
        north_error[i] = np.std(rotated_temp[1])
    return east_error, north_error
