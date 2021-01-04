"""
Main script for launching spinOS
"""
import sys
import getopt

modpath = './modules'
if modpath not in sys.path:
    sys.path.append(modpath)


def hhelp():
    """
    help function
    """
    print('spinOS.py [<dir>]                                             for the GUI')
    print('or')
    print('spinOS.py [<dir>] -i <pointer> [-p] [-s] [-m [-t <steps>]]    for the commandline utility.')


if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "w:hi:psmt:")
    except getopt.GetoptError:
        hhelp()
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            hhelp()
            sys.exit()

    if not opts:
        try:
            wd = args[0]
        except IndexError:
            wd = None
        import modules.spinOSGUI as app

        app.run(wd)
    else:
        import modules.spinOScommandline as cmline

        cmline.run(opts)
