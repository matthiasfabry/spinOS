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
"""
import tkinter as tk
from tkinter import ttk
import sys


def getName(file=None):
    msg = 'Name of new data set' if file is None else 'Name of dataset from \'{}\''.format(file)
    name = tk.simpledialog.askstring('name', msg)
    if name == '' and file is not None:
        name = file
    return name


class VerticalScrolledFrame(tk.Frame):

    def __init__(self, parent):
        tk.Frame.__init__(self, parent)
        self.canvas = tk.Canvas(self, borderwidth=0)
        self.frame = tk.Frame(self.canvas)
        self.vsb = ttk.Scrollbar(self, orient="vertical", command=self.canvas.yview)
        self.canvas.configure(yscrollcommand=self.vsb.set)

        self.vsb.pack(side=tk.RIGHT, fill=tk.Y)
        self.canvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        # TODO: try to center frame within canvas
        self.canvas.create_window((0, 0), window=self.frame, anchor="nw",
                                  tags="self.frame")
        self.canvas.config(yscrollcommand=self.vsb.set, scrollregion=(0, 0, 100, 100))
        self.vsb.lift(self.frame)
        self.frame.bind("<Configure>", self.onFrameConfigure)
        self.frame.bind('<Enter>', self._bound_to_mousewheel)
        self.frame.bind('<Leave>', self._unbound_to_mousewheel)

    def onFrameConfigure(self, event):
        self.canvas.configure(scrollregion=self.canvas.bbox("all"))

    def _bound_to_mousewheel(self, event):
        self.canvas.bind_all("<MouseWheel>", self._on_mousewheel)

    def _unbound_to_mousewheel(self, event):
        self.canvas.unbind_all("<MouseWheel>")

    def _on_mousewheel(self, event):
        if sys.platform.startswith('windows'):
            self.canvas.yview_scroll(int(-1*(event.delta / 120)), "units")
        else:
            self.canvas.yview_scroll(-int(event.delta*2), "units")