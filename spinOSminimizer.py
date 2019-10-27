"""
Module that performs a non-linear least squares minimization of the spectrescopic and the astrometric data
using the Levenberg=Marquardt algorithm
"""


import scipy.optimize as spopt
import numpy as np
from spinOSloader import spinOStag


def LMminimizer(datadict: dict, tag: spinOStag):
    if tag == spinOStag.SB2AS:
        pass
    if tag == spinOStag.SB2:
        pass
    if tag == spinOStag.SB1AS:
        pass
    if tag == spinOStag.SB1:
        pass
    if tag == spinOStag.AS:
        pass

