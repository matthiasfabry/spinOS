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
import tkinter as tk
import tkinter.ttk as ttk
from abc import ABC, abstractmethod
from typing import TYPE_CHECKING

import numpy as np

import modules.constants as cst
import modules.spinOS_io as spl
import modules.utils as util

if TYPE_CHECKING:
    import modules.gui as gui


class DataManager:
    
    def __init__(self, ggui: 'gui.SpinOSGUI'):
        self.gui = ggui
        self.datasets = {'RV1': [], 'RV2': [], 'AS': []}
        self.defWeight = None
    
    def buildSets(self):
        for dataset in self.datasets['RV1']:
            dataset.setData()
        for dataset in self.datasets['RV2']:
            dataset.setData()
        for dataset in self.datasets['AS']:
            dataset.setData()
        self.setDefWeight()
        self.gui.def_as_weight.set(np.round(self.defWeight, 4))
    
    def getBuiltRV1s(self):
        data = None
        for dataset in self.datasets['RV1']:
            if data is None and dataset.getData() is not None:
                data = dataset.getData()
            elif dataset.getData() is not None:
                data = np.vstack((data, dataset.getData()))
        return data
    
    def getBuiltRV2s(self):
        data = None
        for dataset in self.datasets['RV2']:
            if data is None and dataset.getData() is not None:
                data = dataset.getData()
            elif dataset.getData() is not None:
                data = np.vstack((data, dataset.getData()))
        return data
    
    def getBuiltASs(self):
        data = None
        for dataset in self.datasets['AS']:
            if data is None and dataset.getData() is not None:
                data = dataset.getData()
            elif dataset.getData() is not None:
                data = np.vstack((data, dataset.getData()))
        return data
    
    def get_all_data(self):
        datadict = {}
        if self.hasRV1():
            datadict['RV1'] = self.getBuiltRV1s()
        if self.hasRV2():
            datadict['RV2'] = self.getBuiltRV2s()
        if self.hasAS():
            datadict['AS'] = self.getBuiltASs()
        return datadict
    
    def hasRV1(self):
        return len(self.datasets['RV1']) > 0
    
    def hasRV2(self):
        return len(self.datasets['RV2']) > 0
    
    def hasAS(self):
        return len(self.datasets['AS']) > 0
    
    def loaddataintoSets(self):
        """
        loads the data from the currently selected files into datasets
        """
        filetypes = list()
        filenames = list()
        if self.gui.rv1_file.get() != '' and self.gui.load_rv1.get():
            filetypes.append('RV1file')
            filenames.append(self.gui.rv1_file.get())
        if self.gui.rv2_file.get() != '' and self.gui.load_rv2.get():
            filetypes.append('RV2file')
            filenames.append(self.gui.rv2_file.get())
        if self.gui.as_file.get() != '' and self.gui.load_as.get():
            filetypes.append('ASfile')
            filenames.append(self.gui.as_file.get())
        try:
            datahere = spl.data_loader(self.gui.wd.get(), filetypes, filenames)
        except (OSError, KeyError) as e:
            print(e)
            return None
        if self.gui.rv1_file.get() != '' and self.gui.load_rv1.get():
            name = util.getString(
                'Name of dataset from {}'.format(self.gui.rv1_file.get()),
                default=self.gui.rv1_file.get())
            if name is not None:
                newset = RVDataSet(self, self.gui.rv1_tab, self.gui.rv1book,
                                   self.gui.rv1_file.get(),
                                   name, tpe='RV1')
                newset.setentriesfromfile(datahere['RV1'])
                self.datasets['RV1'].append(newset)
                self.gui.plotter.rv1data_lines.append(None)
                self.gui.load_rv1.set(False)
                self.gui.toggle_rv1()
        if self.gui.rv2_file.get() != '' and self.gui.load_rv2.get():
            name = util.getString(
                'Name of dataset from {}'.format(self.gui.rv2_file.get()),
                default=self.gui.rv2_file.get())
            if name is not None:
                newset = RVDataSet(self, self.gui.rv2_tab, self.gui.rv2book,
                                   self.gui.rv2_file.get(),
                                   name, tpe='RV2')
                newset.setentriesfromfile(datahere['RV2'])
                self.datasets['RV2'].append(newset)
                self.gui.load_rv2.set(False)
                self.gui.toggle_rv2()
        if self.gui.as_file.get() != '' and self.gui.load_as.get():
            name = util.getString(
                'Name of dataset from {}'.format(self.gui.as_file.get()),
                default=self.gui.as_file.get())
            if name is not None:
                newset = ASDataSet(self, self.gui.as_tab, self.gui.asbook,
                                   self.gui.as_file.get(),
                                   name, seppa=self.gui.seppa.get())
                newset.setentriesfromfile(datahere['AS'])
                self.datasets['AS'].append(newset)
                self.gui.load_as.set(False)
                self.gui.toggle_as()
        self.gui.set_RV_or_AS_mode()
        self.buildSets()
    
    def emptyrv1DataSet(self):
        newset = RVDataSet(self, self.gui.rv1_tab, self.gui.rv1book,
                           util.getString('Name of new data set',
                                          default='rv1'), tpe='RV1')
        if newset is not None:
            self.datasets['RV1'].append(newset)
            self.gui.toggle_rv1()
    
    def emptyrv2DataSet(self):
        newset = RVDataSet(self, self.gui.rv2_tab, self.gui.rv2book,
                           util.getString('Name of new data set',
                                          default='rv2'), tpe='RV2')
        if newset is not None:
            self.datasets['RV2'].append(newset)
            self.gui.toggle_rv2()
    
    def emptyasDataSet(self, seppa=False):
        newset = ASDataSet(self, self.gui.as_tab, self.gui.asbook,
                           util.getString('Name of new data set',
                                          default='astrometry'), seppa=seppa)
        if newset is not None:
            self.datasets['AS'].append(newset)
            self.gui.toggle_as()
    
    def setDefWeight(self):
        if not self.hasAS():
            self.defWeight = 0
        else:
            numas = sum(
                (len(dataset.getData()) if dataset.getData() is not None else 0
                 for dataset in
                 self.datasets['AS']))
            if not self.hasRV1():
                numrv1 = 0
            else:
                numrv1 = sum(
                    (
                        len(dataset.getData()) if dataset.getData() is not None
                        else 0
                        for dataset in
                        self.datasets['RV1']))
            if not self.hasRV2():
                numrv2 = 0
            else:
                numrv2 = sum(
                    (
                        len(dataset.getData()) if dataset.getData() is not None
                        else 0
                        for dataset in
                        self.datasets['RV2']))
            self.defWeight = numas / (numrv1 + numrv2 + numas)


class DataSet(ABC):
    # TODO: add support for deletion of datasets
    
    def __init__(self, dataman, tab, book, filename, name=None, tpe='AS'):
        if name is None:
            name = filename
        self.dataman: DataManager = dataman
        self.gui = self.dataman.gui
        self.fname = filename
        self.tpe = tpe
        newset = ttk.Frame(tab)
        self.name_var = tk.StringVar()
        self.name_var.set(name)
        self.name_label = ttk.Label(newset, textvariable=self.name_var,
                                    font=('', 13, 'bold'))
        self.name_label.pack()
        self.datagrid = util.VerticalScrolledFrame(newset)
        self.datagrid.pack(fill=tk.BOTH, expand=1)
        self.selall = tk.BooleanVar(value=True)
        tk.Checkbutton(self.datagrid, var=self.selall, bg=cst.BGCOLOR,
                       command=self.selAll).grid()
        but = ttk.Frame(newset)
        ttk.Button(but, text='+', command=self.addentry).grid()
        ttk.Button(but, text='Rename', command=self.setName) \
            .grid(row=0, column=1)
        but.pack()
        self.data = None
        self.book = book
        self.book.add(newset, text=self.fname)
        self.id = len(book.tabs()) - 1
        self.entries = []
        # I add line objects to the plotter but don't bind them to these
        # dataset objects, makes deletion hard
        if self.tpe == 'RV1':
            self.gui.plotter.rv1data_lines.append(None)
        elif self.tpe == 'RV2':
            self.gui.plotter.rv2data_lines.append(None)
        else:
            self.gui.plotter.asdata_lines.append(None)
            self.gui.plotter.as_ellipses.append(None)
            self.gui.plotter.as_dist_lines.append(None)
    
    def setName(self) -> None:
        self.name_var.set(util.getString('new name for this dataset'))
    
    def setData(self) -> None:
        self.data = None
        for entry in self.entries:
            if self.data is None and entry.toInclude():
                self.data = entry.getData()
            elif entry.toInclude():
                self.data = np.vstack((self.data, entry.getData()))
        if self.data is not None and self.data.ndim == 1:
            self.data = self.data.reshape((1, self.data.size))
    
    @abstractmethod
    def setentriesfromfile(self, data) -> None:
        raise NotImplementedError
    
    @abstractmethod
    def addentry(self) -> None:
        raise NotImplementedError
    
    def getData(self) -> np.ndarray:
        return self.data
    
    def selAll(self):
        for entry in self.entries:
            entry.include.set(self.selall.get())
        self.dataman.buildSets()


class RVDataSet(DataSet):
    
    def __init__(self, dataman, tab, book, filename, name=None, **kwargs):
        super().__init__(dataman, tab, book, filename, name, **kwargs)
        ttk.Label(self.datagrid, text='date').grid(row=0, column=1)
        ttk.Label(self.datagrid, text='RV').grid(row=0, column=2)
        ttk.Label(self.datagrid, text='error').grid(row=0, column=3)
    
    def addentry(self):
        self.entries.append(RVEntry(self.datagrid, len(self.entries) + 1))
    
    def setentriesfromfile(self, datadict):
        self.entries = []  # delete all present entries here
        for i in range(len(datadict['hjds'])):
            self.entries.append(
                RVEntry(self.datagrid, i + 1, hjdin=datadict['hjds'][i],
                        rvin=datadict['RVs'][i],
                        errorin=datadict['errors'][i]))


class ASDataSet(DataSet):
    
    def __init__(self, dataman, tab, book, filename, name=None, seppa=False,
                 **kwargs):
        super().__init__(dataman, tab, book, filename, name, **kwargs)
        self.seppa = seppa
        
        ttk.Label(self.datagrid, text='date').grid(row=0, column=1)
        ttk.Label(self.datagrid, text='Sep' if self.seppa else 'East').grid(
            row=0, column=2)
        ttk.Label(self.datagrid, text='PA' if self.seppa else 'North').grid(
            row=0, column=3)
        ttk.Label(self.datagrid, text='major').grid(row=0, column=4)
        ttk.Label(self.datagrid, text='minor').grid(row=0, column=5)
        ttk.Label(self.datagrid, text='error PA').grid(row=0, column=6)
    
    def addentry(self):
        self.entries.append(
            ASEntry(self.datagrid.frame, len(self.entries) + 1))
    
    def setentriesfromfile(self, datadict):
        self.entries = []  # delete all present entries here
        for i in range(len(datadict['hjds'])):
            self.entries.append(
                ASEntry(self.datagrid, i + 1, hjdin=datadict['hjds'][i],
                        eastorsepin=datadict['eastsorsep'][i],
                        northorpain=datadict['northsorpa'][i],
                        majorin=datadict['majors'][i],
                        minorin=datadict['minors'][i],
                        pa_in=datadict['pas'][i], seppa=self.seppa))


class Entry(ABC):
    
    def __init__(self, datagrid, i, hjdin=None):
        self.include = tk.BooleanVar()
        check = tk.Checkbutton(datagrid, bg=cst.BGCOLOR, var=self.include)
        check.grid(row=i, sticky=tk.E)
        self.hjdvar = tk.DoubleVar()
        if hjdin is not None:
            self.hjdvar.set(hjdin)
            self.include.set(True)
        hjd = ttk.Entry(datagrid, textvariable=self.hjdvar, width=10)
        hjd.grid(row=i, column=1)
    
    def toInclude(self):
        return self.include.get()
    
    def getHjd(self):
        return self.hjdvar.get()


class RVEntry(Entry):
    
    def __init__(self, datagrid, i, rvin=None, errorin=None, hjdin=None):
        super().__init__(datagrid, i, hjdin)
        self.rvvar = tk.DoubleVar()
        if rvin is not None:
            self.rvvar.set(rvin)
        rv = ttk.Entry(datagrid, textvariable=self.rvvar, width=10)
        rv.grid(row=i, column=2)
        self.errorvar = tk.DoubleVar()
        if errorin is not None:
            self.errorvar.set(errorin)
        error = ttk.Entry(datagrid, textvariable=self.errorvar, width=10)
        error.grid(row=i, column=3)
    
    def getData(self):
        return np.array([self.getHjd(), self.getRV(), self.getError()])
    
    def getRV(self):
        return self.rvvar.get()
    
    def getError(self):
        return self.errorvar.get()


class ASEntry(Entry):
    
    def __init__(self, datagrid, i, seppa=False,
                 hjdin=None, eastorsepin=None, northorpain=None, majorin=None,
                 minorin=None, pa_in=None):
        super().__init__(datagrid, i, hjdin)
        self.seppa = seppa
        self.eastorsepvar = tk.DoubleVar()
        if eastorsepin is not None:
            self.eastorsepvar.set(eastorsepin)
        eastorsep = ttk.Entry(datagrid, textvariable=self.eastorsepvar,
                              width=5)
        eastorsep.grid(row=i, column=2)
        self.northorpavar = tk.DoubleVar()
        if northorpain is not None:
            self.northorpavar.set(northorpain)
        northorpa = ttk.Entry(datagrid, textvariable=self.northorpavar,
                              width=5)
        northorpa.grid(row=i, column=3)
        self.majorvar = tk.DoubleVar()
        if majorin is not None:
            self.majorvar.set(majorin)
        major = ttk.Entry(datagrid, textvariable=self.majorvar, width=5)
        major.grid(row=i, column=4)
        self.minorvar = tk.DoubleVar()
        if minorin is not None:
            self.minorvar.set(minorin)
        minor = ttk.Entry(datagrid, textvariable=self.minorvar, width=5)
        minor.grid(row=i, column=5)
        self.pavar = tk.DoubleVar()
        if pa_in is not None:
            self.pavar.set(pa_in)
        pa = ttk.Entry(datagrid, textvariable=self.pavar, width=5)
        pa.grid(row=i, column=6)
    
    def getData(self):
        if self.seppa:
            east = self.getEastorsep() * np.sin(
                self.getNorthorpa() * cst.DEG2RAD)
            north = self.getEastorsep() * np.cos(
                self.getNorthorpa() * cst.DEG2RAD)
        else:
            east = self.getEastorsep()
            north = self.getNorthorpa()
        easterror, northerror = spl.convert_error_ellipse(self.getMajor(),
                                                          self.getMinor(),
                                                          self.getPA())
        return np.array([self.getHjd(), east, north, easterror, northerror,
                         self.getMajor(), self.getMinor(), self.getPA()])
    
    def getEastorsep(self):
        return self.eastorsepvar.get()
    
    def getNorthorpa(self):
        return self.northorpavar.get()
    
    def getMajor(self):
        return self.majorvar.get()
    
    def getMinor(self):
        return self.minorvar.get()
    
    def getPA(self):
        return self.pavar.get()
