"""
Copyright 2020, 2021, 2022 Matthias Fabry
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
"""
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
NORMALSIZE = 12
HCOLOR = '#3399ff'
BGCOLOR = '#d9d9d9'
FONTCOLOR = "#000000"
TIME_STR = r'time [day]'
PHASE_STR = r'orbital phase'
PARAM_LIST = ['p', 'e', 'i', 'omega', 'Omega', 't0', 'd', 'k1', 'k2', 'gamma1', 'gamma2', 'mt']
START_LIMS = [-50, 50, -0.15, 1.15, -10, 10, -10, 10]
LIM_STRINGS = ['RV y lower limit',
               'RV y upper limit',
               'RV x lower limit',
               'RV x upper limit',
               'Orbit y lower limit',
               'Orbit y upper limit',
               'Orbit x lower limit',
               'Orbit x upper limit']
RV1COLORS = ['mediumblue', 'indigo', 'fuchsia', 'midnightblue']
RV2COLORS = ['indianred', 'maroon', 'orangered', 'goldenrod']
ASCOLORS = ['darkred', 'green', 'navy', 'mediumvioletred']
ASDISTCOLORS = ['lightcoral', 'limegreen', 'slateblue', 'orchid']
VERSION = "2.7.3"
