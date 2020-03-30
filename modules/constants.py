import numpy as np

au2km = 1.495978707e8  # (km)
pc2km = 3.085677581e13  # (km)
m_sun = 1.9885e30  # (kg)
r_sun = 6.957e5  # (km)
G = 6.67430e-20  # (km3 kg-1 s-2)
deg2rad = np.pi / 180
rad2deg = 180 / np.pi
mas2rad = 1e-3 / 3600 * deg2rad
day2sec = 86400
rad2mas = rad2deg * 3600 * 1e3
