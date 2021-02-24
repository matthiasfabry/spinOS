"""
Copyright 2020, 2021 Matthias Fabry
This file is part of spinOS.

spinOS is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

spinOS is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with spinOS.  If not, see <https://www.gnu.org/licenses/>.


Module that handles the loading of the relevant data for the solver.
"""
import numpy as np


def guess_loader(wd: str, guessfile: str) -> dict:
    """
    parses the guess file and determines values and flags for each guess
    :param wd: the working directory
    :param guessfile: pathname (relative to wd) pointing to the file containing guesses
    :return: dictionary containing the guesses and flags for each parameter
    """
    wd = check_slash(wd)
    guesses = np.genfromtxt(wd + guessfile, dtype=None, filling_values=np.nan, usecols=(0, 1, 2),
                            encoding='utf-8')
    guessdict = dict()
    for i in range(12):
        guessdict[guesses[i][0]] = (guesses[i][1], guesses[i][2])
    return guessdict


def guess_saver(wd: str, name: str, guess_dict: dict) -> None:
    """
    saves guesses to a file
    :param wd: working directory
    :param name: file name
    :param guess_dict: guesses to save
    """
    wd = check_slash(wd)
    with open(wd + name + '.txt', 'w') as guessfile:
        for param, guess in guess_dict.items():
            guessfile.write(param + ' {} {}\n'.format(guess[0], str(guess[1])))


def data_loader(wd: str, filetypes: list, filenames: list) -> dict:
    """
    loads data from files into a dictionary
    :param wd: working directory where the files are
    :param filetypes: data types to load, must be 'RV1file', 'RV2file', or 'ASfile'
    :param filenames: names of the files in question
    :return: data in a dictionary
    """
    wd = check_slash(wd)
    data_dict = dict()
    for i in range(len(filetypes)):
        if filetypes[i] == 'RV1file':
            data = np.loadtxt(wd + filenames[i])
            data_dict['RV1'] = dict()
            data_dict['RV1']['hjds'] = data[:, 0]
            data_dict['RV1']['RVs'] = data[:, 1]
            try:
                data_dict['RV1']['errors'] = data[:, 2]
            except IndexError:
                # put dummy error if none found in data
                data_dict['RV1']['errors'] = data[:, 1] * 0.05
        elif filetypes[i] == 'RV2file':
            data = np.loadtxt(wd + filenames[i])
            data_dict['RV2'] = dict()
            data_dict['RV2']['hjds'] = data[:, 0]
            data_dict['RV2']['RVs'] = data[:, 1]
            try:
                data_dict['RV2']['errors'] = data[:, 2]
            except IndexError:
                # put dummy error if none found in data
                data_dict['RV2']['errors'] = data[:, 1] * 0.05
        elif filetypes[i] == 'ASfile':
            data = np.loadtxt(wd + filenames[i])
            data_dict['AS'] = dict()
            data_dict['AS']['hjds'] = data[:, 0]
            data_dict['AS']['majors'] = data[:, 3]
            data_dict['AS']['minors'] = data[:, 4]
            data_dict['AS']['pas'] = data[:, 5]
            data_dict['AS']['eastsorsep'] = data[:, 1]
            data_dict['AS']['northsorpa'] = data[:, 2]
    return data_dict


def convert_error_ellipse(major, minor, angle):
    """
    Converts error ellipses to actual east and north errors by a sampling the error ellipse monte-carlo style and
    then taking the variance in the east and north directions.
    :param major: length of the major axis of the error ellipse
    :param minor: length of the minor axis of the error ellipse
    :param angle: position angle east of north of the major axis
    :return: east and north error
    """
    num = 1000
    cosa = np.cos(angle)
    sina = np.sin(angle)
    temp_major = np.random.randn(num) * major
    temp_minor = np.random.randn(num) * minor
    rotated_temp = np.matmul(np.array([[cosa, sina], [-sina, cosa]]), [temp_major, temp_minor])
    east_error = np.std(rotated_temp[0])
    north_error = np.std(rotated_temp[1])
    return east_error, north_error


def check_slash(wd):
    if len(wd) == 0:
        return wd
    if wd[-1] != '/':
        wd += '/'
    return wd
