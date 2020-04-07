import tkinter as tk
import time


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
        # noinspection PyProtectedMember
        self.__window.tk.call("::tk::unsupported::MacWindowStyle", "style", self.__window._w, "plain", "none")
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

