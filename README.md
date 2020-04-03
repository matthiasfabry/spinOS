<b>
presenting spinOS: the SPectroscopic and INterferometric Orbital Solution finder.
</b>

<i>Goal:</i>

spinOS computes the best fit orbital solution given:
 1) radial and systemic velocity data for either or both components of a spectroscopic binary and/or
 2) astrometric data containing the separtions and error ellipses, and
 3) an initial guess of the parameters:
     
        - e        the eccentricity
        - i        the inclination (deg)
        - omega    the longitude of the periastron with respect to the ascending node of the secondary (deg)
        - Omega    the longitude of the ascending node of the seconday measured east of north (deg)
        - t0       the time of periastron passage (day)
        - k1       the semiamplitude of the radial velocity curve of the primary (km/s)
        - k2       the semiamplitude of the radial velocity curve of the secondary (km/s)
        - p        the period of the binary (day)
        - gamma1   the peculiar velocity of the primary (km/s)
        - gamma2   the peculiar velocity of the secondary (km/s)
        - d        the distance (pc)
       
    For each parameter, a tag True/False should be supplied to decide for the minimizer to vary this parameter. (So
    False means it will keep it fixed at the supplied value)

This package provides a commandline script interface as well as a GUI (recommended) to easily visualize your data

<i>Usage:</i>

To use the GUI, simply run:
    
    python spinOS.py [working directory (optional)]

For the command line version, information on usage can be found in the docstring of spinOScommandline.py

<i>Dependencies:</i>

    python 3.7.7
    numpy 1.18.1
    scipy 1.4.1
    lmfit 1.0.0
    matplotlib 3.1.3
    emcee 3.0.2  (for MCMC calculations)
    corner 2.0.1 (for plotting of an MCMC diagram)

<i>Author:</i>

    Matthias Fabry
    Instituut voor Sterrekunde, KU Leuven, Belgium

<i>Date:</i>

    3 Apr 2020

<i>Version:</i>

    2.3

<i>Acknowledgements:</i>

    This python3 implementation is heavily based on an earlier IDL implementation by Hugues Sana.
    We thank the authors of lmfit for the development of their package.
