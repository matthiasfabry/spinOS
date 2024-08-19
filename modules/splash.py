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


Defines the splash screen class
"""

import tkinter as tk
import time
from sys import platform


class Splash:

    def __init__(self, rroot, file, wait, width, height):
        self.__root = rroot
        self.__file = file
        self.__wait = wait + time.time()
        self.__scrW, self.__scrH = width, height

    def __enter__(self):
        # hide main window
        self.__root.withdraw()
        self.__window = tk.Toplevel(self.__root)
        self.__window.lift()
        self.__splash = tk.PhotoImage(master=self.__window, file=self.__file)
        if platform == 'darwin':
            # noinspection PyProtectedMember
            self.__window.tk.call("::tk::unsupported::MacWindowStyle", "style", self.__window._w,
                                  "plain", "none")
        # geometry
        imgW = self.__splash.width()
        imgH = self.__splash.height()
        Xpos = (self.__scrW - imgW) // 2
        Ypos = (self.__scrH - imgH) // 2
        self.__window.geometry('+{}+{}'.format(Xpos, Ypos))
        # put image
        tk.Label(self.__window, image=self.__splash).grid()
        # display before .mainloop()
        self.__window.update_idletasks()
        self.__window.overrideredirect(True)
        self.__window.update()

    def __exit__(self, *args):
        if time.time() < self.__wait:
            time.sleep(self.__wait - time.time())
        del self.__splash
        self.__window.destroy()
        self.__root.update_idletasks()
        self.__root.deiconify()
