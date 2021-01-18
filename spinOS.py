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

Main script for launching spinOS
"""
import sys
import getopt


def hhelp():
    """
    help function
    """
    print('spinOS.py [<dir>]                                             for the GUI')
    print('or')
    print('spinOS.py [<dir>] -i <pointer> [-p] [-s] [-m [-t <steps>]]    for the commandline utility.')


try:
    opts, args = getopt.getopt(sys.argv[1:], "w:hi:psmt:")
except getopt.GetoptError:
    hhelp()
    sys.exit(2)

for opt, arg in opts:
    if opt == '-h':
        hhelp()
        sys.exit()

if __name__ == '__main__':
    if not opts:
        try:
            wd = args[0]
        except IndexError:
            wd = None
        import modules.spinOSGUI as gui

        gui.run(wd)
    else:
        import modules.spinOScommandline as cm

        cm.run(opts)
