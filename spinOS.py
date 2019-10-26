"""
presenting spinOS: the SPectroscopic and INterferometric Orbital Solution finder.

author:
    Matthias Fabry
    Instituut voor Sterrekunde, KU Leuven, Belgium

date:
    23 Oct 2019

version:
    alpha

goal:
    spinOS computes the best fit orbital solution given:
     1) radial and systemic velocity data for either or both components of a spectroscopic binary and/or
     2) astrometric data containing the separtions and position angles, and
     3) an initial guess of the parameters:
        a) a        the semi major axis,
        b) e        the eccentricity
        c) i        the inclination
        d) omega    the longitude of the periastron with respect to the ascending node of the secondary
        e) Omega    the longitude of the ascending node east of north
        f) T        the time of periastron passage
        g) P        the period of the binary
        h) gamma1   the (apparent) systemic velocity of the primary
        i) gamma2   the (apparent) systemic velocity of the secondary
        j) d        the distance

acknowledgements:
    This python3 implementation is heavily based on an earlier IDL implementation by Hugues Sana.

"""
import sys
import numpy as np
import spinOSloader

# read in files
paramfile = sys.argv[1]
data_dict = spinOSloader.spinOSparser(paramfile)


# compute best elements
    # compute lots of models
    # do chisq with data


# compute model of these elements


# plot resulting RV curve


# plot resulting apparent orbit


# print best orbital params

