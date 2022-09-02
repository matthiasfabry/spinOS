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
import pathlib
import tkinter as tk
from tkinter import ttk

import lmfit as lm
import numpy as np

import modules.binary_system as bsys
import modules.constants as cst
import modules.minimizer as spm
import modules.spinOS_io as spl
import modules.splash as splash
import modules.utils as util
from modules.data_manager import DataManager
from modules.plotting import Plotting


class SpinOSGUI:
    """
    class specifying the main spinOS tk implementation
    """
    
    # IMPROVEMENT: fix fonts to some standard one?
    def __init__(self, master, wwd, width, height):
        
        # set homogenous style
        s = ttk.Style()
        s.theme_use('default')
        s.configure("TFrame", background=cst.BGCOLOR)
        s.configure("TEntry", padding=4)
        s.configure("TButton", padding=5)
        s.configure("TLabel", padding=3)
        
        # FRAME STRUCTURE #
        # set the root frame
        tabs = ttk.Notebook(master)
        
        # set the data frame
        data_frame_tab = ttk.Frame(tabs)
        data_frame = util.VerticalScrolledFrame(data_frame_tab)
        data_frame.pack(expand=1, fill=tk.BOTH, anchor=tk.N)
        
        # set the guess frame
        guess_infer_tab = ttk.Frame(tabs)
        guess_infer_top = util.VerticalScrolledFrame(guess_infer_tab)
        infer_frame = ttk.Frame(guess_infer_top)
        guess_frame = ttk.Frame(guess_infer_top)
        guess_frame.pack()
        refreshframe1 = ttk.Frame(guess_infer_top)
        refreshframe1.pack(pady=10)
        infer_frame.pack()
        guess_infer_top.pack(expand=1, fill=tk.BOTH, anchor=tk.N)
        
        # set the minimization frame
        min_frame_tab = ttk.Frame(tabs)
        min_frame = util.VerticalScrolledFrame(min_frame_tab)
        min_frame.pack(expand=1, fill=tk.BOTH, anchor=tk.N)
        
        # set the plot window controls frame
        plt_frame_tab = ttk.Frame(tabs)
        plt_frame_top = util.VerticalScrolledFrame(plt_frame_tab)
        plt_frame = ttk.Frame(plt_frame_top)
        plt_frame_top.pack(expand=1, fill=tk.BOTH, anchor=tk.N)
        
        # GLOBAL GUI VARS #
        self.w, self.h = width, height
        self.wd = wwd
        
        # the plotter instance
        self.plotter = Plotting(self)
        
        # DATA FRAME #
        self.datamanager = DataManager(self)
        filesframe = ttk.Frame(data_frame)
        firstlabel = ttk.Label(filesframe, text='DATA',
                               font=('', cst.TITLESIZE, 'underline'))
        firstlabel.grid(columnspan=5, sticky=tk.N)
        
        # define inlcusion variables
        self.load_rv1 = tk.BooleanVar()
        self.load_rv2 = tk.BooleanVar()
        self.load_as = tk.BooleanVar()
        
        # assign to checkbuttons
        rv1check = ttk.Checkbutton(filesframe, var=self.load_rv1,
                                   command=self.toggle_rv1)
        rv2check = ttk.Checkbutton(filesframe, var=self.load_rv2,
                                   command=self.toggle_rv2)
        ascheck = ttk.Checkbutton(filesframe, var=self.load_as,
                                  command=self.toggle_as)
        
        # put them in a nice grid
        rv1check.grid(row=2, sticky=tk.E)
        rv2check.grid(row=3, sticky=tk.E)
        ascheck.grid(row=4, sticky=tk.E)
        
        # define labels
        ttk.Label(filesframe, text='load?').grid(row=1, column=0, sticky=tk.E)
        ttk.Label(filesframe, text='Working directory').grid(row=1, column=1,
                                                             sticky=tk.E)
        self.rv1_label = ttk.Label(filesframe, text='Primary RV file',
                                   state=tk.DISABLED)
        self.rv1_label.grid(row=2, column=1, sticky=tk.E)
        self.rv2_label = ttk.Label(filesframe, text='Secondary RV file',
                                   state=tk.DISABLED)
        self.rv2_label.grid(row=3, column=1, sticky=tk.E)
        self.as_label = ttk.Label(filesframe, text='Astrometric data file',
                                  state=tk.DISABLED)
        self.as_label.grid(row=4, column=1, sticky=tk.E)
        
        # define entries
        self.wd = ttk.Entry(filesframe)
        self.rv1_file = ttk.Entry(filesframe)
        self.rv2_file = ttk.Entry(filesframe)
        self.as_file = ttk.Entry(filesframe)
        
        # put some mock values
        if wwd:
            self.wd.insert(0, wwd)
        self.rv1_file.insert(0, 'primary_vels.txt')
        self.rv2_file.insert(0, 'secondary_vels.txt')
        self.as_file.insert(0, 'relative_astrometry.txt')
        
        # disable them after inserting stuff
        self.rv1_file.config(state=tk.DISABLED)
        self.rv2_file.config(state=tk.DISABLED)
        self.as_file.config(state=tk.DISABLED)
        
        # put in a nice grid
        self.wd.grid(row=1, column=2)
        self.rv1_file.grid(row=2, column=2)
        self.rv2_file.grid(row=3, column=2)
        self.as_file.grid(row=4, column=2)
        
        self.seppa = tk.BooleanVar(value=True)
        self.seppa_but = tk.Radiobutton(filesframe, text='Sep/PA',
                                        variable=self.seppa, value=True,
                                        fg=cst.FONTCOLOR, bg=cst.BGCOLOR,
                                        state=tk.DISABLED)
        self.seppa_but.grid(row=4, column=3)
        self.en_but = tk.Radiobutton(filesframe, text='E/N',
                                     variable=self.seppa, value=False,
                                     fg=cst.FONTCOLOR, bg=cst.BGCOLOR,
                                     state=tk.DISABLED)
        self.en_but.grid(row=4, column=4, sticky=tk.W)
        
        ttk.Button(filesframe, text='Load Data',
                   command=self.datamanager.loaddataintoSets,
                   ).grid(row=5, columnspan=5, pady=10)
        filesframe.grid_columnconfigure(0, weight=1)
        filesframe.grid_columnconfigure(4, weight=1)
        filesframe.pack(expand=1, fill=tk.BOTH, anchor=tk.N)
        
        data_manager_frame = ttk.Frame(data_frame)
        ttk.Label(data_manager_frame, text='Data Manager',
                  font=('', cst.TITLESIZE, 'underline')).pack()
        managernotebook = ttk.Notebook(data_manager_frame)
        self.rv1_tab = ttk.Frame(managernotebook)
        self.rv1book = ttk.Notebook(self.rv1_tab)
        self.rv1book.pack(expand=1, fill=tk.BOTH)
        ttk.Button(self.rv1_tab, text='Add Dataset',
                   command=self.datamanager.emptyrv1DataSet).pack()
        managernotebook.add(self.rv1_tab, text='Primary RVs')
        
        self.rv2_tab = ttk.Frame(managernotebook)
        self.rv2book = ttk.Notebook(self.rv2_tab)
        self.rv2book.pack(expand=1, fill=tk.BOTH)
        ttk.Button(self.rv2_tab, text='Add Dataset',
                   command=self.datamanager.emptyrv2DataSet).pack()
        managernotebook.add(self.rv2_tab, text='Secondary RVs')
        
        self.as_tab = ttk.Frame(managernotebook)
        self.asbook = ttk.Notebook(self.as_tab)
        self.asbook.pack(expand=1, fill=tk.BOTH)
        self.seppa_dtmgr = tk.BooleanVar(value=True)
        self.seppa_but_dtmgr = tk.Radiobutton(self.as_tab, text='Sep/PA',
                                              variable=self.seppa,
                                              bg=cst.BGCOLOR, fg=cst.FONTCOLOR,
                                              value=True)
        self.seppa_but_dtmgr.pack(side=tk.LEFT, expand=1)
        self.en_but_dtmgr = tk.Radiobutton(self.as_tab, text='E/N',
                                           variable=self.seppa,
                                           bg=cst.BGCOLOR, fg=cst.FONTCOLOR,
                                           value=False)
        self.en_but_dtmgr.pack(side=tk.LEFT, expand=1)
        ttk.Button(self.as_tab, text='Add Dataset',
                   command=lambda: self.datamanager.emptyasDataSet(
                       seppa=self.seppa_dtmgr.get())) \
            .pack(side=tk.LEFT, expand=1)
        managernotebook.add(self.as_tab, text='Astrometry')
        
        managernotebook.pack(expand=1, fill=tk.BOTH, side=tk.TOP)
        
        data_manager_frame.pack(expand=1, fill=tk.BOTH, side=tk.TOP)
        
        # GUESS FRAME #
        self.system = None
        self.param_dict = None
        self.guess_dict = None
        columns = 7
        labelcolumn = 1
        entrycolumn = 2
        varycheckcolumn = 3
        transfercolumn = 4
        minresultcolumn = 5
        errorcolumn = 6
        
        numofparams = 12
        rparams = range(numofparams)
        
        titlesrow = 2
        paramgridrow = titlesrow + 1
        buttonrow = paramgridrow + numofparams
        
        # print the labels in the guess frame
        ttk.Label(guess_frame, text='SYSTEM PARAMETERS',
                  font=('', cst.TITLESIZE, 'underline')).grid(
            columnspan=columns)
        
        ttk.Label(guess_frame, text='Guess file').grid(row=1,
                                                       column=labelcolumn,
                                                       sticky=tk.E)
        self.guess_file = ttk.Entry(guess_frame, width=15)
        self.guess_file.insert(0, 'guesses.txt')
        self.guess_file.grid(row=1, column=entrycolumn, sticky=tk.S,
                             columnspan=2)
        
        ttk.Label(guess_frame, text='Guesses').grid(row=titlesrow,
                                                    column=entrycolumn)
        ttk.Label(guess_frame, text='Vary?').grid(row=titlesrow,
                                                  column=varycheckcolumn)
        ttk.Label(guess_frame, text='Transfer').grid(row=titlesrow,
                                                     column=transfercolumn)
        ttk.Label(guess_frame, text='Result').grid(row=titlesrow,
                                                   column=minresultcolumn)
        ttk.Label(guess_frame, text='Error').grid(row=titlesrow,
                                                  column=errorcolumn)
        
        self.lock_gs = tk.BooleanVar(False)
        self.locked_image = tk.PhotoImage(
            file=pathlib.Path(__file__).parent.parent.joinpath('rsc/lock.png'))
        self.unlocked_image = tk.PhotoImage(
            file=pathlib.Path(__file__).parent.parent.joinpath(
                'rsc/unlock.png'))
        self.lock_gs_button = ttk.Button(guess_frame, image=self.locked_image,
                                         command=self.toggle_lock)
        self.lock_gs_button.grid(row=paramgridrow + 10)
        
        self.q_mode = tk.BooleanVar(False)
        self.lock_q_button = ttk.Button(guess_frame, width=1, text='q',
                                        command=self.toggle_q)
        self.lock_q_button.grid(row=paramgridrow + 8)
        
        self.param_var_list = [tk.StringVar() for _ in rparams]
        self.param_var_list[0].set('p (days) =')
        self.param_var_list[1].set('e =')
        self.param_var_list[2].set('i (deg) =')
        self.param_var_list[3].set('omega (deg) =')
        self.param_var_list[4].set('Omega (deg) =')
        self.param_var_list[5].set('t0 (JD) =')
        self.param_var_list[6].set('d (pc) =')
        self.param_var_list[7].set('k1 (km/s) =')
        self.param_var_list[8].set('k2 (km/s) =')
        self.param_var_list[9].set('gamma1 (km/s) =')
        self.param_var_list[10].set('gamma2 (km/s) =')
        self.param_var_list[11].set('M_tot (Msun) =')
        
        self.param_label_list = [
            ttk.Label(guess_frame, textvariable=self.param_var_list[i]) for i
            in
            rparams]
        
        for i in rparams:
            self.param_label_list[i].grid(row=paramgridrow + i,
                                          column=labelcolumn, sticky=tk.E)
        
        # initialize the entry variables
        self.guess_var_list = [tk.StringVar(value='0') for _ in rparams]
        
        # define entry boxes
        self.guess_entry_list = [
            ttk.Entry(guess_frame, textvariable=self.guess_var_list[i],
                      width=10) for i in rparams]
        # put in a nice grid
        for i in rparams:
            self.guess_entry_list[i].grid(row=paramgridrow + i,
                                          column=entrycolumn)
        
        # define the vary state variables
        self.vary_var_list = [tk.BooleanVar() for _ in rparams]
        
        # define checkbuttons for vary states
        self.vary_button_list = [
            ttk.Checkbutton(guess_frame, var=self.vary_var_list[i]) for i in
            rparams]
        
        # put the checkbuttons in a nice grid
        for i in rparams:
            self.vary_button_list[i].grid(row=paramgridrow + i,
                                          column=varycheckcolumn)
        
        # define the transfer buttons
        # for this semantic to work, we need to wrap the lambda function
        # into another one, so that
        # each command
        # references to its own number 'y', rather than the outer 'i' of the
        # list comprehension
        self.transfer_button_list = [ttk.Button(guess_frame, text='<-',
                                                command=(lambda y: (
                                                    lambda: self.transfer(y)))(
                                                    i)).grid(
            row=paramgridrow + i,
            column=transfercolumn) for i in
                                     rparams]
        
        # define the minimized parameter variables
        self.mininimzed_var_list = [tk.StringVar() for _ in rparams]
        
        # define the labels the minimized parameters will go in
        self.min_label_list = [
            ttk.Label(guess_frame, textvariable=self.mininimzed_var_list[i],
                      width=8) for i in
            rparams]
        
        for i in rparams:
            self.min_label_list[i].grid(row=paramgridrow + i,
                                        column=minresultcolumn)
        
        # define the error variables
        self.error_var_list = [tk.StringVar() for _ in rparams]
        
        # define the labels the errors will go in
        self.error_label_list = [
            ttk.Label(guess_frame, textvariable=self.error_var_list[i],
                      width=8)
            for i in rparams]
        for i in rparams:
            self.error_label_list[i].grid(row=paramgridrow + i,
                                          column=errorcolumn)
        
        # define the buttons in this frame
        ttk.Button(guess_frame, text='Load guesses',
                   command=self.load_guesses).grid(row=buttonrow,
                                                   column=labelcolumn)
        ttk.Button(guess_frame, text='Save guesses', command=self.save_guesses,
                   ).grid(row=buttonrow, column=entrycolumn)
        ttk.Button(guess_frame, text='Save parameters',
                   command=self.save_params,
                   ).grid(row=buttonrow, column=minresultcolumn,
                          columnspan=2)
        ttk.Button(refreshframe1, text='Refresh Plots & Inferred Parameters',
                   command=self.update).pack()
        
        # INFER FRAME #
        self.mprimary = tk.StringVar()
        self.msecondary = tk.StringVar()
        self.semimajord = tk.StringVar()
        self.semimajork1k2 = tk.StringVar()
        self.totalmass = tk.StringVar()
        
        # define labels
        ttk.Label(infer_frame, text='INFERRED PARAMETERS (from guesses)',
                  font=('', cst.TITLESIZE, 'underline')).grid(columnspan=4,
                                                              sticky=tk.N)
        ttk.Label(infer_frame, text='From k1/k2',
                  font=('', 13, 'underline')).grid(row=1,
                                                   columnspan=2)
        ttk.Label(infer_frame, text='M1 (M_sun) =').grid(row=3, sticky=tk.E)
        ttk.Label(infer_frame, text='M2 (M_sun) =').grid(row=4, sticky=tk.E)
        ttk.Label(infer_frame, text='M (M_sun) =').grid(row=5, sticky=tk.E)
        ttk.Label(infer_frame, text='a (AU) =').grid(row=2, sticky=tk.E)
        ttk.Label(infer_frame, textvariable=self.mprimary).grid(row=3,
                                                                column=1)
        ttk.Label(infer_frame, textvariable=self.msecondary).grid(row=4,
                                                                  column=1)
        ttk.Label(infer_frame, textvariable=self.semimajork1k2).grid(row=2,
                                                                     column=1)
        ttk.Label(infer_frame, textvariable=self.totalmass).grid(row=5,
                                                                 column=1)
        ttk.Separator(infer_frame).grid(column=2, row=2, rowspan=5,
                                        sticky=tk.NS)
        ttk.Label(infer_frame, text='From d/M_tot:',
                  font=('', 13, 'underline')).grid(row=1,
                                                   column=3,
                                                   columnspan=2)
        ttk.Label(infer_frame, text='a (AU) =').grid(row=2, column=3,
                                                     sticky=tk.E)
        ttk.Label(infer_frame, textvariable=self.semimajord).grid(row=2,
                                                                  column=4)
        
        # MINIMIZATION FRAME #
        self.minimization_run_number = 0
        self.minresult = None
        self.didmcmc = False
        self.method = tk.StringVar(value='leastsq')
        self.redchisq = tk.DoubleVar()
        self.dof = tk.IntVar()
        self.steps = tk.IntVar(value=1000)
        self.walkers = tk.IntVar(value=100)
        self.burn = tk.IntVar(value=100)
        self.thin = tk.IntVar(value=1)
        
        self.rms_rv1 = tk.DoubleVar()
        self.rms_rv2 = tk.DoubleVar()
        self.rms_as = tk.DoubleVar()
        self.do_custom_weight = tk.BooleanVar(value=False)
        self.def_as_weight = tk.DoubleVar()
        self.custom_as_weight = tk.DoubleVar()
        self.hops = tk.IntVar()
        
        # define labels and buttons in a grid
        methodframe = ttk.Frame(min_frame)
        ttk.Label(methodframe, text='MINIMIZATION',
                  font=('', cst.TITLESIZE, 'underline')).grid(
            columnspan=4)
        ttk.Label(methodframe, text='Method:').grid(sticky=tk.E)
        tk.Radiobutton(methodframe, text='Levenberg-Marquardt',
                       variable=self.method,
                       fg=cst.FONTCOLOR,
                       value='leastsq', bg=cst.BGCOLOR,
                       command=self.toggle_method). \
            grid(row=1, column=1, sticky=tk.W)
        tk.Radiobutton(methodframe, text='Basinhopping', variable=self.method,
                       value='basinhopping', bg=cst.BGCOLOR, fg=cst.FONTCOLOR,
                       command=self.toggle_method).grid(row=2, column=1,
                                                        sticky=tk.W)
        tk.Radiobutton(methodframe, text='LM+MCMC', variable=self.method,
                       value='emcee',
                       bg=cst.BGCOLOR, fg=cst.FONTCOLOR,
                       command=self.toggle_method).grid(row=3, column=1,
                                                        sticky=tk.W)
        self.hops_label = ttk.Label(methodframe, text='# of hops:',
                                    state=tk.DISABLED)
        self.hops_label.grid(row=2, column=2)
        self.hops_entry = ttk.Entry(methodframe, textvariable=self.hops,
                                    width=5, state=tk.DISABLED)
        self.hops_entry.grid(row=2, column=3)
        methodframe.pack()
        mcframe = ttk.Frame(min_frame)
        ttk.Label(mcframe, text='MCMC params: ').grid(row=0, column=1)
        self.steps_label = ttk.Label(mcframe, text='# of steps:',
                                     state=tk.DISABLED)
        self.steps_label.grid(row=0, column=2, sticky=tk.E)
        self.steps_entry = ttk.Entry(mcframe, textvariable=self.steps, width=5,
                                     state=tk.DISABLED)
        self.steps_entry.grid(row=0, column=3)
        self.walkers_label = ttk.Label(mcframe, text='# of walkers:',
                                       state=tk.DISABLED)
        self.walkers_label.grid(row=0, column=4, sticky=tk.E)
        self.walkers_entry = ttk.Entry(mcframe, textvariable=self.walkers,
                                       width=5,
                                       state=tk.DISABLED)
        self.walkers_entry.grid(row=0, column=5)
        self.burn_label = ttk.Label(mcframe, text='Burn:', state=tk.DISABLED)
        self.burn_label.grid(row=0, column=6, sticky=tk.E)
        self.burn_entry = ttk.Entry(mcframe, textvariable=self.burn, width=5,
                                    state=tk.DISABLED)
        self.burn_entry.grid(row=0, column=7)
        self.thin_label = ttk.Label(mcframe, text='Thin:', state=tk.DISABLED)
        self.thin_label.grid(row=0, column=8, sticky=tk.E)
        self.thin_entry = ttk.Entry(mcframe, textvariable=self.thin, width=5,
                                    state=tk.DISABLED)
        self.thin_entry.grid(row=0, column=9)
        mcframe.pack()
        self.mc_widg = {self.steps_label, self.steps_entry, self.walkers_entry,
                        self.walkers_label,
                        self.burn_entry, self.burn_label, self.thin_label,
                        self.thin_entry}
        otherminframe = ttk.Frame(min_frame)
        ttk.Label(otherminframe, text='astrometric weight from data = ').grid(
            row=2, column=1,
            sticky=tk.E)
        self.def_weight_label = ttk.Label(otherminframe,
                                          textvariable=self.def_as_weight)
        self.def_weight_label.grid(row=2, column=2, sticky=tk.W)
        self.as_weight_button = ttk.Checkbutton(otherminframe,
                                                var=self.do_custom_weight,
                                                command=self.toggle_weights)
        self.as_weight_button.grid(row=3, sticky=tk.E)
        self.weight_label = ttk.Label(otherminframe,
                                      text='Custom astrometric weight =',
                                      state=tk.DISABLED)
        self.weight_label.grid(row=3, column=1, sticky=tk.E)
        self.weight_slider = tk.Scale(otherminframe,
                                      variable=self.custom_as_weight,
                                      from_=0, to=1, orient=tk.HORIZONTAL,
                                      resolution=0.001,
                                      state=tk.DISABLED, length=180,
                                      bg=cst.BGCOLOR,
                                      fg=cst.FONTCOLOR)
        self.weight_slider.grid(row=3, column=2, columnspan=2, sticky=tk.W)
        
        ttk.Button(otherminframe, text='Minimize!',
                   command=self.minimize).grid(row=4, columnspan=4)
        ttk.Label(otherminframe, text='Results',
                  font=('', cst.TITLESIZE, 'underline')) \
            .grid(row=5, columnspan=4)
        ttk.Label(otherminframe, text='Red. Chi Sqrd =').grid(row=6,
                                                              sticky=tk.E)
        ttk.Label(otherminframe, text='Deg. of frdm =').grid(row=7,
                                                             sticky=tk.E)
        ttk.Label(otherminframe, textvariable=self.redchisq).grid(row=6,
                                                                  column=1,
                                                                  sticky=tk.W)
        ttk.Label(otherminframe, textvariable=self.dof).grid(row=7, column=1,
                                                             sticky=tk.W)
        ttk.Label(otherminframe, text='RMS Primary (km/s) =').grid(row=6,
                                                                   column=2,
                                                                   sticky=tk.E)
        ttk.Label(otherminframe, text='RMS Secondary (km/s) =').grid(row=7,
                                                                     column=2,
                                                                     sticky=tk.E)
        ttk.Label(otherminframe, text='RMS Rel. Orbit (mas) =').grid(row=8,
                                                                     column=2,
                                                                     sticky=tk.E)
        ttk.Label(otherminframe, textvariable=self.rms_rv1).grid(row=6,
                                                                 column=3,
                                                                 sticky=tk.W)
        ttk.Label(otherminframe, textvariable=self.rms_rv2).grid(row=7,
                                                                 column=3,
                                                                 sticky=tk.W)
        ttk.Label(otherminframe, textvariable=self.rms_as).grid(row=8,
                                                                column=3,
                                                                sticky=tk.W)
        self.min_save_button = ttk.Button(otherminframe,
                                          text='Save minimization result',
                                          command=self.save_params,
                                          state=tk.DISABLED)
        self.min_save_button.grid(row=9, columnspan=4)
        self.mcplotbutton = ttk.Button(otherminframe,
                                       text='Make MCMC scatterplot matrix',
                                       command=self.plotter.make_corner_diagram,
                                       state=tk.DISABLED)
        self.mcplotbutton.grid(row=10, columnspan=4)
        otherminframe.pack()
        
        # PLOT CONTROLS #
        ttk.Label(plt_frame, text='PLOT CONTROLS',
                  font=('', cst.TITLESIZE, 'underline')).grid(
            columnspan=6)
        # UI elements
        self.phase_label = ttk.Label(plt_frame, text='phase =',
                                     state=tk.DISABLED)
        self.phase_label.grid(row=1, column=1, sticky=tk.E)
        self.phase_slider = tk.Scale(plt_frame, variable=self.plotter.phase,
                                     from_=0, to=1,
                                     orient=tk.HORIZONTAL, bg=cst.BGCOLOR,
                                     fg=cst.FONTCOLOR,
                                     resolution=0.005, length=300,
                                     state=tk.DISABLED)
        self.phase_slider.grid(row=1, column=2, columnspan=4)
        self.phase_button = ttk.Checkbutton(plt_frame,
                                            var=self.plotter.do_phasedot,
                                            command=self.toggle_dot,
                                            state=tk.DISABLED)
        self.phase_button.grid(row=1)
        
        self.plot_rv1data_label = ttk.Label(plt_frame, text='Primary RV data',
                                            state=tk.DISABLED)
        self.plot_rv1data_label.grid(row=2, column=1)
        self.plot_rv1data_button = ttk.Checkbutton(plt_frame,
                                                   var=self.plotter.do_datarv1,
                                                   state=tk.DISABLED)
        self.plot_rv1data_button.grid(row=2)
        
        self.plot_rv2data_label = ttk.Label(plt_frame,
                                            text='Secondary RV data',
                                            state=tk.DISABLED)
        self.plot_rv2data_label.grid(row=3, column=1)
        self.plot_rv2data_button = ttk.Checkbutton(plt_frame,
                                                   var=self.plotter.do_datarv2,
                                                   state=tk.DISABLED)
        self.plot_rv2data_button.grid(row=3)
        
        self.plot_asdata_label = ttk.Label(plt_frame, text='Astrometric data',
                                           state=tk.DISABLED)
        self.plot_asdata_label.grid(row=4, column=1)
        self.plot_asdata_button = ttk.Checkbutton(plt_frame,
                                                  var=self.plotter.do_dataas,
                                                  state=tk.DISABLED)
        self.plot_asdata_button.grid(row=4)
        
        self.plot_rv1model_label = ttk.Label(plt_frame,
                                             text='Primary RV model',
                                             state=tk.DISABLED)
        self.plot_rv1model_label.grid(row=2, column=3)
        self.plot_rv1model_button = ttk.Checkbutton(plt_frame,
                                                    var=self.plotter.do_modelrv1,
                                                    state=tk.DISABLED)
        self.plot_rv1model_button.grid(row=2, column=2)
        
        self.plot_gamma1_label = ttk.Label(plt_frame, text=r'Gamma 1',
                                           state=tk.DISABLED)
        self.plot_gamma1_label.grid(row=3, column=3)
        self.plot_gamma1_button = ttk.Checkbutton(plt_frame,
                                                  var=self.plotter.do_modelgamma1,
                                                  state=tk.DISABLED)
        self.plot_gamma1_button.grid(row=3, column=2)
        
        self.plot_rv2model_label = ttk.Label(plt_frame,
                                             text='Secondary RV model',
                                             state=tk.DISABLED)
        self.plot_rv2model_label.grid(row=4, column=3)
        self.plot_rv2model_button = ttk.Checkbutton(plt_frame,
                                                    var=self.plotter.do_modelrv2,
                                                    state=tk.DISABLED)
        self.plot_rv2model_button.grid(row=4, column=2)
        
        self.plot_gamma2_label = ttk.Label(plt_frame, text='Gamma 2',
                                           state=tk.DISABLED)
        self.plot_gamma2_label.grid(row=5, column=3)
        self.plot_gamma2_button = ttk.Checkbutton(plt_frame,
                                                  var=self.plotter.do_modelgamma2,
                                                  state=tk.DISABLED)
        self.plot_gamma2_button.grid(row=5, column=2)
        
        self.plot_asmodel_label = ttk.Label(plt_frame, text='Model Orbit',
                                            state=tk.DISABLED)
        self.plot_asmodel_label.grid(row=2, column=5)
        self.plot_asmodel_button = ttk.Checkbutton(plt_frame,
                                                   var=self.plotter.do_modelas,
                                                   state=tk.DISABLED)
        self.plot_asmodel_button.grid(row=2, column=4)
        
        self.plot_nodeline_label = ttk.Label(plt_frame, text='Line of nodes',
                                             state=tk.DISABLED)
        self.plot_nodeline_label.grid(row=3, column=5)
        self.plot_nodeline_button = ttk.Checkbutton(plt_frame,
                                                    var=self.plotter.do_nodeline,
                                                    state=tk.DISABLED)
        self.plot_nodeline_button.grid(row=3, column=4)
        
        self.plot_semimajor_label = ttk.Label(plt_frame,
                                              text='Semi-major axis',
                                              state=tk.DISABLED)
        self.plot_semimajor_label.grid(row=4, column=5)
        self.plot_semimajor_button = ttk.Checkbutton(plt_frame,
                                                     var=self.plotter.do_semimajor,
                                                     state=tk.DISABLED)
        self.plot_semimajor_button.grid(row=4, column=4)
        
        self.plot_peri_label = ttk.Label(plt_frame, text='Periastron',
                                         state=tk.DISABLED)
        self.plot_peri_label.grid(row=5, column=5)
        self.plot_peri_button = ttk.Checkbutton(plt_frame,
                                                var=self.plotter.do_peri,
                                                state=tk.DISABLED)
        self.plot_peri_button.grid(row=5, column=4)
        
        self.as_dist_label = ttk.Label(plt_frame, text='Astrometric errors',
                                       state=tk.DISABLED)
        self.as_dist_label.grid(row=6, column=5)
        self.as_dist_button = ttk.Checkbutton(plt_frame,
                                              var=self.plotter.do_as_dist,
                                              state=tk.DISABLED)
        self.as_dist_button.grid(row=6, column=4)
        
        legend_button = ttk.Checkbutton(plt_frame, var=self.plotter.do_legend)
        legend_button.grid(row=6)
        ttk.Label(plt_frame, text='Legend').grid(row=6, column=1)
        
        self.pphase_but = tk.Radiobutton(plt_frame, text='phase',
                                         command=self.toggle_phase_time,
                                         bg=cst.BGCOLOR, fg=cst.FONTCOLOR,
                                         variable=self.plotter.plot_vs_phase,
                                         value=True,
                                         state=tk.DISABLED)
        self.ptime_but = tk.Radiobutton(plt_frame, text='time',
                                        command=self.toggle_phase_time,
                                        bg=cst.BGCOLOR, fg=cst.FONTCOLOR,
                                        variable=self.plotter.plot_vs_phase,
                                        value=False,
                                        state=tk.DISABLED)
        self.pphase_but.grid(row=6, column=2)
        self.ptime_but.grid(row=6, column=3)
        self.modelwidgets = {self.plot_asmodel_label, self.plot_asmodel_button,
                             self.plot_rv1model_button,
                             self.plot_rv2model_button,
                             self.plot_rv1model_label,
                             self.plot_rv2model_label,
                             self.plot_semimajor_button,
                             self.plot_semimajor_label,
                             self.plot_nodeline_button,
                             self.plot_nodeline_label,
                             self.plot_peri_label,
                             self.plot_peri_button, self.as_dist_button,
                             self.as_dist_label,
                             self.pphase_but, self.ptime_but,
                             self.plot_gamma1_button,
                             self.plot_gamma1_label, self.plot_gamma2_button,
                             self.plot_gamma2_label}
        plt_frame.pack()
        
        ttk.Button(plt_frame_top, text='Refresh Plots', command=self.update,
                   ).pack()
        
        settings_frame = ttk.Frame(plt_frame_top)
        entrycol = 1
        ttk.Label(settings_frame, text='PLOT SETTINGS',
                  font=('', cst.TITLESIZE, 'underline')).grid(
            columnspan=3)
        
        ttk.Label(settings_frame, text='Axis label size').grid(row=1,
                                                               sticky=tk.E)
        ttk.Entry(settings_frame, textvariable=self.plotter.axeslabelsize,
                  width=15) \
            .grid(row=1, column=entrycol, columnspan=2)
        
        ttk.Label(settings_frame, text='Tick label size').grid(row=2,
                                                               sticky=tk.E)
        ttk.Entry(settings_frame, textvariable=self.plotter.ticklabelsize,
                  width=15) \
            .grid(row=2, column=entrycol, columnspan=2)
        
        ttk.Label(settings_frame, text='Grid visibility').grid(row=3,
                                                               sticky=tk.E)
        
        self.grid_button = ttk.Checkbutton(settings_frame,
                                           var=self.plotter.do_grids)
        self.grid_button.grid(row=3, column=1)
        
        ttk.Label(settings_frame, text='Axis limits').grid(row=4)
        tk.Radiobutton(settings_frame, text='Auto',
                       var=self.plotter.limcontrol, value=True,
                       bg=cst.BGCOLOR, fg=cst.FONTCOLOR,
                       command=self.togglelimcontrol).grid(row=4, column=1)
        tk.Radiobutton(settings_frame, text='Manual',
                       var=self.plotter.limcontrol, value=False,
                       bg=cst.BGCOLOR, fg=cst.FONTCOLOR,
                       command=self.togglelimcontrol).grid(row=4, column=2)
        
        limitrow = 5
        self.limit_labels = []
        for i in range(8):
            self.limit_labels.append(
                ttk.Label(settings_frame, text=cst.LIM_STRINGS[i],
                          state=tk.DISABLED))
        
        for i in range(8):
            self.limit_labels[i].grid(row=limitrow + i, sticky=tk.E)
        
        self.limit_entries = []
        
        for i in range(8):
            self.limit_entries.append(
                ttk.Entry(settings_frame, textvariable=self.plotter.limits[i],
                          state=tk.DISABLED))
        for i in range(8):
            self.limit_entries[i].grid(row=limitrow + i, column=entrycol,
                                       columnspan=2)
        
        settings_frame.pack(pady=10)
        
        ttk.Button(plt_frame_top, text='Refresh Plots',
                   command=self.plotter.update_plots).pack(pady=10)
        ttk.Button(plt_frame_top, text='New plot windows',
                   command=self.plotter.init_plots).pack(pady=10)
        
        # setup the plotting windows
        self.plotter.init_plots()
        
        # display in the root frame
        tabs.add(data_frame_tab, text='Data Files')
        tabs.add(guess_infer_tab, text='System/Parameters')
        tabs.add(min_frame_tab, text='Minimization')
        tabs.add(plt_frame_tab, text='Plot Controls')
        tabs.pack(expand=1, fill=tk.BOTH)
    
    @staticmethod
    def toggle(widg, boolvalue):
        """
        toggles widget widg to be disabled or not given a boolvalue
        :param widg: widget to toggle
        :param boolvalue: bool
        """
        if boolvalue:
            widg.config(state=tk.NORMAL)
        else:
            widg.config(state=tk.DISABLED)
    
    def toggle_phase_time(self):
        if not self.plotter.plot_vs_phase.get() and \
                self.plotter.do_phasedot.get():
            self.plotter.do_phasedot.set(False)
        self.toggle(self.phase_button, self.plotter.plot_vs_phase.get())
    
    def togglelimcontrol(self):
        for widg in self.limit_labels:
            self.toggle(widg, not self.plotter.limcontrol.get())
        for widg in self.limit_entries:
            self.toggle(widg, not self.plotter.limcontrol.get())
        if not self.plotter.limcontrol.get():
            self.plotter.matchLimits()
    
    def toggle_rv1(self):
        """
        toggles the RV1 widgets
        """
        for widg in self.rv1_file, self.rv1_label:
            self.toggle(widg, self.load_rv1.get())
        for widg in self.plot_rv1data_label, self.plot_rv1data_button:
            self.toggle(widg, self.datamanager.hasRV1())
    
    def toggle_rv2(self):
        """
        toggles the RV2 widgets
        """
        for widg in self.rv2_file, self.rv2_label:
            self.toggle(widg, self.load_rv2.get())
        for widg in self.plot_rv2data_button, self.plot_rv2data_label:
            self.toggle(widg, self.datamanager.hasRV2())
    
    def toggle_as(self):
        """
        toggles the AS widgets
        """
        for widg in self.as_file, self.as_label, self.seppa_but, self.en_but:
            self.toggle(widg, self.load_as.get())
        for widg in self.plot_asdata_label, self.plot_asdata_button:
            self.toggle(widg, self.datamanager.hasAS())
    
    def toggle_dot(self):
        for widg in self.phase_slider, self.phase_label:
            self.toggle(widg, self.plotter.do_phasedot.get())
    
    def toggle_q(self):
        self.q_mode.set(not self.q_mode.get())
        if self.q_mode.get():
            self.param_var_list[8].set('q = k1/k2 =')
            self.toggle(self.vary_button_list[8], False)
            self.lock_q_button.config(text='k')
            if float(self.guess_var_list[8].get()) != 0:
                self.guess_var_list[8].set(np.round(
                    float(self.guess_var_list[7].get()) / float(
                        self.guess_var_list[8].get()), 3))
        else:
            self.param_var_list[8].set('k2 (km/s) =')
            self.toggle(self.vary_button_list[8], True)
            self.lock_q_button.config(text='q')
            if float(self.guess_var_list[8].get()) != 0:
                self.guess_var_list[8].set(np.round(
                    float(self.guess_var_list[7].get()) / float(
                        self.guess_var_list[8].get()), 3))
    
    def toggle_lock(self):
        self.lock_gs.set(not self.lock_gs.get())
        if self.lock_gs.get():
            self.lock_gs_button.config(image=self.unlocked_image)
            for widgset in self.param_label_list, self.vary_button_list, \
                           self.guess_entry_list:
                self.toggle(widgset[10], False)
        
        else:
            self.lock_gs_button.config(image=self.locked_image)
            for widgset in self.param_label_list, self.vary_button_list, \
                           self.guess_entry_list:
                self.toggle(widgset[10], True)
            self.guess_var_list[10].set(str(self.guess_var_list[9].get()))
    
    def toggle_method(self):
        """
        toggles the appropriate method widgets
        """
        for widg in self.mc_widg:
            if self.method.get() == 'emcee':
                self.toggle(widg, True)
            else:
                self.toggle(widg, False)
        for widg in self.hops_label, self.hops_entry:
            if self.method.get() == 'basinhopping':
                self.toggle(widg, True)
            else:
                self.toggle(widg, False)
    
    def toggle_weights(self):
        """
        toggles the weight widgets
        """
        if not (self.datamanager.hasAS() and (
                self.datamanager.hasRV1() or self.datamanager.hasRV2())):
            self.do_custom_weight.set(False)
        
        self.toggle(self.weight_label, self.do_custom_weight.get())
        self.toggle(self.weight_slider, self.do_custom_weight.get())
    
    def set_RV_or_AS_mode(self):
        """
        sets the parameters in the correct inference mode
        """
        for lst in self.param_label_list, self.vary_button_list:
            if self.datamanager.hasRV1() and self.datamanager.hasRV2() and \
                    self.datamanager.hasAS():
                for i in {2, 4, 6, 7, 8, 9, 10, 11}:
                    lst[i].config(state=tk.NORMAL)
            elif self.datamanager.hasRV1() and self.datamanager.hasRV2():
                for i in {7, 8, 9, 10}:
                    lst[i].config(state=tk.NORMAL)
                for i in {2, 4, 6, 11}:
                    lst[i].config(state=tk.DISABLED)
            elif self.datamanager.hasRV1() and self.datamanager.hasAS():
                for i in {2, 4, 6, 7, 9, 11}:
                    lst[i].config(state=tk.NORMAL)
                for i in {8, 10}:
                    lst[i].config(state=tk.DISABLED)
            elif self.datamanager.hasRV1():
                for i in {7, 9, 11}:
                    lst[i].config(state=tk.NORMAL)
                for i in {2, 4, 6, 8, 10, 11}:
                    lst[i].config(state=tk.DISABLED)
            elif self.datamanager.hasAS():
                for i in {2, 4, 6, 11}:
                    lst[i].config(state=tk.NORMAL)
                for i in {7, 8, 9, 10}:
                    lst[i].config(state=tk.DISABLED)
            else:
                for i in {2, 4, 6, 7, 8, 9, 10}:
                    lst[i].config(state=tk.NORMAL)
                lst[11].config(state=tk.DISABLED)
    
    def transfer(self, varno):
        """
        pushes a minimization result to the parameter column
        :param varno: number in the parameter list
        """
        self.guess_var_list[varno].set(self.mininimzed_var_list[varno].get())
    
    def load_guesses(self):
        """
        load guesses from a file to the guess column
        """
        try:
            self.guess_dict = spl.guess_loader(self.wd.get(),
                                               self.guess_file.get())
        except IOError:
            print('cannot find your guess file!')
            self.guess_dict = None
            return
        except ValueError as e:
            print(
                'your guessfile seems to have badly formatted data or '
                'something...')
            print(e)
            self.guess_dict = None
            return
        try:
            for i in range(len(cst.PARAM_LIST)):
                self.guess_var_list[i].set(
                    self.guess_dict[cst.PARAM_LIST[i]][0])
                self.vary_var_list[i].set(
                    str(self.guess_dict[cst.PARAM_LIST[i]][1]))
        except (ValueError, TypeError) as e:
            print('some parameter has not been set properly:', e)
            self.guess_dict = None
            return
        self.set_system()
    
    def set_guess_dict_from_entries(self):
        """
        builds the guess dict from the guess column
        """
        self.guess_dict = {}
        for i in range(len(cst.PARAM_LIST)):
            if i == 8:
                self.guess_dict[cst.PARAM_LIST[i]] = (
                    float(self.guess_var_list[
                              8].get()) if not self.q_mode.get() else float(
                        self.guess_var_list[7].get()) / float(
                        self.guess_var_list[8].get()),
                    self.vary_var_list[8].get())
            else:
                self.guess_dict[cst.PARAM_LIST[i]] = (
                    float(self.guess_var_list[i].get()),
                    self.vary_var_list[i].get())
        
        self.param_dict = dict()
        for param, value in self.guess_dict.items():
            self.param_dict[param] = value[0]
    
    def set_system(self):
        """
        sets the system from the current guess column
        """
        try:
            self.set_guess_dict_from_entries()
            self.system = bsys.BinarySystem(self.param_dict)
        except ValueError:
            print('invalid model!')
            self.guess_dict = None
            self.system = None
            self.plotter.plot_vs_phase.set(False)
            self.toggle_phase_time()
            for widg in self.modelwidgets:
                self.toggle(widg, False)
            return False
        else:
            self.toggle_phase_time()
            for widg in self.modelwidgets:
                self.toggle(widg, True)
            return True
    
    def minimize(self):
        """
        launches a minimization run
        """
        self.set_guess_dict_from_entries()
        self.datamanager.buildSets()
        data_dict = self.datamanager.get_all_data()
        if self.guess_dict is not None and len(data_dict) > 0:
            # calculate best parameters
            try:
                if self.do_custom_weight.get():
                    w = self.custom_as_weight.get()
                else:
                    w = None
                self.minresult, rms_rv1, rms_rv2, rms_as \
                    = spm.LMminimizer(self.guess_dict, data_dict,
                                      self.method.get(),
                                      self.hops.get(), self.steps.get(),
                                      self.walkers.get(),
                                      self.burn.get(), self.thin.get(), w,
                                      self.lock_gs.get(),
                                      self.q_mode.get())
                if self.method.get() == 'emcee':
                    self.didmcmc = True
                    self.toggle(self.mcplotbutton, True)
                else:
                    self.didmcmc = False
                self.minimization_run_number += 1
                self.toggle(self.min_save_button, True)
                pars = self.minresult.params
                # fill in the entries
                self.mininimzed_var_list[0].set(np.round(pars['p'].value, 3))
                if pars['p'].vary:
                    self.error_var_list[0].set(
                        np.round(0 if pars['p'].stderr is None else pars[
                            'p'].stderr, 3))
                self.mininimzed_var_list[1].set(np.round(pars['e'].value, 3))
                if pars['e'].vary:
                    self.error_var_list[1].set(
                        np.round(0 if pars['e'].stderr is None else pars[
                            'e'].stderr, 3))
                self.mininimzed_var_list[2].set(
                    np.round(pars['i'].value % 360, 3))
                if pars['i'].vary:
                    self.error_var_list[2].set(
                        np.round(0 if pars['i'].stderr is None else pars[
                            'i'].stderr, 3))
                self.mininimzed_var_list[3].set(
                    np.round(pars['omega'].value % 360, 3))
                if pars['omega'].vary:
                    self.error_var_list[3].set(
                        np.round(0 if pars['omega'].stderr is None else pars[
                            'omega'].stderr, 3))
                self.mininimzed_var_list[4].set(
                    np.round(pars['Omega'].value % 360, 3))
                if pars['Omega'].vary:
                    self.error_var_list[4].set(
                        np.round(0 if pars['Omega'].stderr is None else pars[
                            'Omega'].stderr, 3))
                self.mininimzed_var_list[5].set(np.round(pars['t0'].value, 3))
                if pars['t0'].vary:
                    self.error_var_list[5].set(
                        np.round(0 if pars['t0'].stderr is None else pars[
                            't0'].stderr, 3))
                self.mininimzed_var_list[6].set(np.round(pars['d'].value, 3))
                if pars['d'].vary:
                    self.error_var_list[6].set(
                        np.round(0 if pars['d'].stderr is None else pars[
                            'd'].stderr, 3))
                self.mininimzed_var_list[7].set(np.round(pars['k1'].value, 3))
                if pars['k1'].vary:
                    self.error_var_list[7].set(
                        np.round(0 if pars['k1'].stderr is None else pars[
                            'k1'].stderr, 3))
                if self.q_mode.get():
                    self.mininimzed_var_list[8].set(
                        np.round(pars['q'].value, 3))
                else:
                    self.mininimzed_var_list[8].set(
                        np.round(pars['k2'].value, 3))
                    if pars['k2'].vary:
                        self.error_var_list[8].set(
                            np.round(0 if pars['k2'].stderr is None else pars[
                                'k2'].stderr, 3))
                self.mininimzed_var_list[9].set(
                    np.round(pars['gamma1'].value, 3))
                if pars['gamma1'].vary:
                    self.error_var_list[9].set(
                        np.round(0 if pars['gamma1'].stderr is None else pars[
                            'gamma1'].stderr, 3))
                self.mininimzed_var_list[10].set(
                    np.round(pars['gamma2'].value, 3))
                if pars['gamma2'].vary:
                    self.error_var_list[10].set(
                        np.round(0 if pars['gamma2'].stderr is None else pars[
                            'gamma2'].stderr, 3))
                
                self.mininimzed_var_list[11].set(np.round(pars['mt'].value, 3))
                if pars['mt'].vary:
                    self.error_var_list[11].set(
                        np.round(0 if pars['mt'].stderr is None else pars[
                            'mt'].stderr, 3))
                
                self.redchisq.set(np.round(self.minresult.redchi, 4))
                self.dof.set(self.minresult.nfree)
                self.rms_rv1.set(np.round(rms_rv1, 4))
                self.rms_rv2.set(np.round(rms_rv2, 4))
                self.rms_as.set(np.round(rms_as, 4))
                self.minimization_run_number += 1
            except ValueError as e:
                print(e)
    
    def set_inferred_params(self):
        self.mprimary.set(np.round(self.system.primary_mass(), 2))
        self.msecondary.set(np.round(self.system.secondary_mass(), 2))
        self.totalmass.set(np.round(self.system.total_mass(), 2))
        self.semimajork1k2.set(
            np.round(self.system.semimajor_axis_from_RV(), 2))
        self.semimajord.set(
            np.round(self.system.semimajor_axis_from_distance(), 2))
    
    def update(self):
        """
        updates the gui, by replotting everything that is selected
        """
        self.datamanager.buildSets()
        if not self.set_system():
            return
        if self.system is not None:
            self.set_inferred_params()
        self.plotter.update_plots()
    
    def save_params(self):
        """
        save minimized parameters to a file
        """
        out = util.getString('name your parameters file',
                             default='fitted_params')
        wd = spl.check_slash(self.wd.get())
        with open(wd + out + '{}.txt'.format(self.minimization_run_number),
                  'w') as f:
            if self.minresult is not None:
                f.write(lm.fit_report(self.minresult))
                f.write('\n')
            #
            # f.write('reduced chisq = {} \n'.format(self.redchisq.get()))
            # f.write('dof = {} \n'.format(self.dof.get()))
        if self.didmcmc:
            np.savetxt(wd + out + '{}_flatchain.txt'.format(
                self.minimization_run_number),
                       self.minresult.flatchain,
                       header='param order: {}'.format(
                           self.minresult.var_names))
        else:
            np.savetxt(
                wd + out + '_covar{}.txt'.format(self.minimization_run_number),
                self.minresult.covar,
                header='param order in this covar matrix: {}'.format(
                    self.minresult.var_names))
    
    def save_guesses(self):
        """
        save guesses parameters to a file
        """
        out = util.getString('name your guesses file', default='guess_save')
        self.set_guess_dict_from_entries()
        spl.guess_saver(self.wd.get(), out, self.guess_dict)


def run(wd):
    """
    Run the main spinOS app tk loop.
    :param wd: working directory to set as root.
    """
    root = tk.Tk()
    wdir = pathlib.PurePath(__file__).parent.parent
    w, h = root.winfo_screenwidth(), root.winfo_screenheight()
    with splash.Splash(root, wdir.joinpath('rsc/spinos100.png'), 2.1, w, h):
        root.geometry("{}x{}+0+0".format(int(0.37 * w),
                                         int(0.95 * h)))  # TODO: on linux this might not scale
        # properly
        root.title('spinOS v{}'.format(cst.VERSION))
        SpinOSGUI(root, wd, w, h)
    
    root.mainloop()
