import numpy as np

AU2KM = 1.495978707e8  # (km)
PC2KM = 3.085677581e13  # (km)
MSUN = 1.9885e30  # (kg)
RSUN = 6.957e5  # (km)
G = 6.67430e-20  # (km3 kg-1 s-2)
DEG2RAD = np.pi / 180
RAD2DEG = 180 / np.pi
MAS2RAD = 1e-3 / 3600 * DEG2RAD
DAY2SEC = 86400
RAD2MAS = RAD2DEG * 3600 * 1e3
TITLESIZE = 20
HCOLOR = '#3399ff'
TIME_STR = r'time [day]'
PHASE_STR = r'orbital phase'
PARAM_LIST = ['p', 'e', 'i', 'omega', 'Omega', 't0', 'd', 'k1', 'k2', 'gamma1', 'gamma2', 'mt']
