"""
module that contains the main spinOSApp class and tk loop
"""
import pathlib
import tkinter as tk
from tkinter import ttk
import time

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import EllipseCollection

from modules import spinOSio as spl, spinOSminimizer as spm, spinOSplotter as spp, binary_system as bsys
from modules import spinOSsplash as splash


class SpinOSApp:
    """
    class specifying the main spinOS tk implementation
    """

    def __init__(self, master, wwd, width, heigth):

        mpl.use("TkAgg")  # set the backend

        # define the structure of the frames
        # set the root frame
        tabs = ttk.Notebook(master)
        # set the data frame
        data_frame_wrap = tk.Frame(tabs)
        data_frame = tk.Frame(data_frame_wrap)
        # set the guess frame
        guess_infer_wrap = tk.Frame(tabs)
        guess_infer_frame = tk.Frame(guess_infer_wrap)
        # set the inferation frame
        infer_frame = tk.Frame(guess_infer_frame)
        guess_frame = tk.Frame(guess_infer_frame)
        # make tabs for the bottom two panels
        min_frame_wrap = tk.Frame(tabs)
        min_frame = tk.Frame(min_frame_wrap)
        # set the plot window controls frame
        plt_frame_wrap = tk.Frame(tabs)
        plt_frame = tk.Frame(plt_frame_wrap)

        # initialize some variables
        self.since = 0
        hcolor = '#3399ff'
        self._w, self._h = width, heigth
        # dictionaries
        self.param_dict = None
        self.guess_dict = None
        self.data_dict = None

        # model system
        self.system = None

        # minimization data
        self.minimization_run_number = 0
        self.mcmc_run_number = 0
        self.minresult = None
        self.didmcmc = False

        # figure and line objects
        self.rv_fig = None
        self.as_fig = None
        self.rv_ax = None
        self.as_ax = None
        self.rv1_dot = None
        self.rv2_dot = None
        self.as_dot = None
        self.rv1_line = None
        self.rv2_line = None
        self.as_line = None
        self.as_dist_lines = None
        self.rv1data_line = None
        self.rv2data_line = None
        self.asdata_line = None
        self.peri_dot = None
        self.node_line = None
        self.semi_major = None
        self.as_ellipses = None
        self.as_legend = None
        self.rv_legend = None

        titlesize = 20

        # DATA FRAME #
        firstlabel = tk.Label(data_frame, text='DATA', font=('', titlesize, 'underline'))
        firstlabel.grid(columnspan=5, sticky=tk.N)

        # define inlcusion variables
        self.include_rv1 = tk.BooleanVar()
        self.include_rv2 = tk.BooleanVar()
        self.include_as = tk.BooleanVar()
        self.loading_guesses = False

        # assign to checkbuttons
        rv1check = tk.Checkbutton(data_frame, var=self.include_rv1, command=self.toggle_rv1)
        rv2check = tk.Checkbutton(data_frame, var=self.include_rv2, command=self.toggle_rv2)
        ascheck = tk.Checkbutton(data_frame, var=self.include_as, command=self.toggle_as)

        # put them in a nice grid
        rv1check.grid(row=2)
        rv2check.grid(row=3)
        ascheck.grid(row=4)

        # define labels
        tk.Label(data_frame, text='Working directory').grid(row=1, column=1, sticky=tk.E)
        self.rv1_label = tk.Label(data_frame, text='Primary RV file', state=tk.DISABLED)
        self.rv1_label.grid(row=2, column=1, sticky=tk.E)
        self.rv2_label = tk.Label(data_frame, text='Secondary RV file', state=tk.DISABLED)
        self.rv2_label.grid(row=3, column=1, sticky=tk.E)
        self.as_label = tk.Label(data_frame, text='Astrometric data file', state=tk.DISABLED)
        self.as_label.grid(row=4, column=1, sticky=tk.E)
        tk.Label(data_frame, text='Guess file').grid(row=5, column=1, sticky=tk.E)

        # define entries
        self.wd = tk.Entry(data_frame)
        self.rv1_file = tk.Entry(data_frame)
        self.rv2_file = tk.Entry(data_frame)
        self.as_file = tk.Entry(data_frame)
        self.guess_file = tk.Entry(data_frame)

        # put some mock values
        if wwd:
            self.wd.insert(0, wwd + '/')
        self.rv1_file.insert(0, 'primary_vels.txt')
        self.rv2_file.insert(0, 'secondary_vels.txt')
        self.as_file.insert(0, 'relative_astrometry.txt')
        self.guess_file.insert(0, 'guesses.txt')

        # disable them
        self.rv1_file.config(state=tk.DISABLED)
        self.rv2_file.config(state=tk.DISABLED)
        self.as_file.config(state=tk.DISABLED)

        # put in a nice grid
        self.wd.grid(row=1, column=2)
        self.rv1_file.grid(row=2, column=2)
        self.rv2_file.grid(row=3, column=2)
        self.as_file.grid(row=4, column=2)
        self.guess_file.grid(row=5, column=2)

        self.seppa = tk.BooleanVar(value=True)
        self.seppa.trace_add('write', lambda n, ix, m: self.update())
        self.seppa_but = tk.Radiobutton(data_frame, text='Sep/PA', variable=self.seppa, value=True, state=tk.DISABLED)
        self.seppa_but.grid(row=4, column=3)
        self.en_but = tk.Radiobutton(data_frame, text='E/N', variable=self.seppa, value=False, state=tk.DISABLED)
        self.en_but.grid(row=4, column=4)

        self.data_button = tk.Button(data_frame, text='Load data', command=self.load_data, state=tk.DISABLED)
        self.data_button.grid(row=6, column=2)

        # GUESS FRAME #
        columns = 6
        paramcolumn = 1
        varycheckcolumn = 2
        transfercolumn = 3
        minresultcolumn = 4
        errorcolumn = 5
        numofparams = 12

        # print the labels in the guess frame
        tk.Label(guess_frame, text='MODEL/GUESS PARAMETERS', font=('', titlesize, 'underline')).grid(columnspan=columns)
        tk.Label(guess_frame, text='Vary?').grid(row=1, column=varycheckcolumn)
        tk.Label(guess_frame, text='Result').grid(row=1, column=minresultcolumn)
        tk.Label(guess_frame, text='Error').grid(row=1, column=errorcolumn)

        self.param_name_vars = [tk.StringVar() for _ in range(numofparams)]
        self.param_name_vars[0].set('p (days) =')
        self.param_name_vars[1].set('e =')
        self.param_name_vars[2].set('i (deg) =')
        self.param_name_vars[3].set('omega (deg) =')
        self.param_name_vars[4].set('Omega (deg) =')
        self.param_name_vars[5].set('t0 (JD) =')
        self.param_name_vars[6].set('d (pc) =')
        self.param_name_vars[7].set('k1 (km/s) =')
        self.param_name_vars[8].set('k2 (km/s) =')
        self.param_name_vars[9].set('gamma1 (km/s) =')
        self.param_name_vars[10].set('gamma2 (km/s) =')
        self.param_name_vars[11].set('M_tot (Msun) =')

        self.param_labels = [tk.Label(guess_frame, textvariable=self.param_name_vars[i]) for i in
                             range(numofparams)]

        for i in range(numofparams):
            self.param_labels[i].grid(row=(i + 2), sticky=tk.E)

        # initialize the entry variables
        self.guess_var_list = [tk.StringVar(value='0') for _ in range(numofparams)]

        # define entry boxes
        self.guess_entry_list = [tk.Entry(guess_frame, textvariable=self.guess_var_list[i], width=10) for i in
                                 range(numofparams)]
        self.guess_entry_list[0].insert(0, '1')
        self.guess_entry_list[2].insert(0, '9')
        self.guess_entry_list[6].insert(0, '100')
        # put in a nice grid
        for i in range(numofparams):
            self.guess_entry_list[i].grid(row=(i + 2), column=paramcolumn)

        # add tracers so the model is updated
        for i in range(numofparams):
            self.guess_var_list[i].trace_add('write', lambda n, ix, m: self.update())

        # define the vary state variables
        self.vary_var_list = [tk.BooleanVar() for _ in range(numofparams)]

        # define checkbuttons for vary states
        self.vary_button_list = [tk.Checkbutton(guess_frame, var=self.vary_var_list[i]) for i in
                                 range(numofparams)]

        # put the checkbuttons in a nice grid
        for i in range(numofparams):
            self.vary_button_list[i].grid(row=(i + 2), column=varycheckcolumn)

        # define the transfer buttons
        # for this semantic to work, we need to wrap the lambda function into another one, so that each command
        # references to its own number 'y', rather than the outer 'i' of the list comprehension
        [tk.Button(guess_frame, text='<-', command=(lambda y: (lambda: self.transfer(y)))(i),
                   highlightbackground=hcolor).grid(row=(i + 2), column=transfercolumn)
         for i in range(numofparams)]

        # define the minimized parameter variables
        self.mininimzed_var_list = [tk.StringVar() for _ in range(numofparams)]

        # define the labels the minimized parameters will go in
        self.min_label_list = [tk.Label(guess_frame, textvariable=self.mininimzed_var_list[i], width=10) for i in
                               range(numofparams)]

        for i in range(numofparams):
            self.min_label_list[i].grid(row=(i + 2), column=minresultcolumn)

        # define the error variables
        self.error_var_list = [tk.StringVar() for _ in range(numofparams)]

        # define the labels the errors will go in
        self.error_label_list = [tk.Label(guess_frame, textvariable=self.error_var_list[i], width=10) for i in
                                 range(numofparams)]
        for i in range(numofparams):
            self.error_label_list[i].grid(row=(i + 2), column=errorcolumn)

        # define the buttons in this frame
        tk.Button(guess_frame, text='Load guesses', command=self.load_guesses,
                  highlightbackground=hcolor).grid(row=numofparams + 2)
        tk.Button(guess_frame, text='Save guesses', command=self.save_guesses,
                  highlightbackground=hcolor).grid(row=numofparams + 2, column=1)
        tk.Button(guess_frame, text='Refresh All', command=self.update, highlightbackground=hcolor).grid(
            row=numofparams + 3, column=1)
        tk.Button(guess_frame, text='Save parameters', command=self.save_params, highlightbackground=hcolor).grid(
            row=numofparams + 2, column=4, columnspan=2)

        # INFER FRAME #
        # define variables
        self.mprimary = tk.StringVar()
        self.msecondary = tk.StringVar()
        self.semimajord = tk.StringVar()
        self.semimajork1k2 = tk.StringVar()
        self.totalmass = tk.StringVar()

        # define labels
        tk.Label(infer_frame, text='INFERRED PARAMETERS', font=('', titlesize, 'underline')) \
            .grid(columnspan=4, sticky=tk.N)
        tk.Label(infer_frame, text='From k1/k2', font=('', 13, 'underline')).grid(row=1, columnspan=2)
        tk.Label(infer_frame, text='M1 (M_sun) =').grid(row=3, sticky=tk.E)
        tk.Label(infer_frame, text='M2 (M_sun) =').grid(row=4, sticky=tk.E)
        tk.Label(infer_frame, text='M (M_sun) =').grid(row=5, sticky=tk.E)
        tk.Label(infer_frame, text='a (AU) =').grid(row=2, sticky=tk.E)
        tk.Label(infer_frame, textvariable=self.mprimary).grid(row=3, column=1)
        tk.Label(infer_frame, textvariable=self.msecondary).grid(row=4, column=1)
        tk.Label(infer_frame, textvariable=self.semimajork1k2).grid(row=2, column=1)
        tk.Label(infer_frame, textvariable=self.totalmass).grid(row=5, column=1)
        ttk.Separator(infer_frame).grid(column=2, row=2, rowspan=5, sticky=tk.NS)
        tk.Label(infer_frame, text='From d/M_tot:', font=('', 13, 'underline')).grid(row=1, column=3, columnspan=2)
        tk.Label(infer_frame, text='a (AU) =').grid(row=2, column=3, sticky=tk.E)
        tk.Label(infer_frame, textvariable=self.semimajord).grid(row=2, column=4)

        # MINIMIZATION FRAME #
        # define variables
        self.do_mcmc = tk.BooleanVar()
        self.redchisq = tk.DoubleVar()
        self.dof = tk.IntVar()
        self.steps = tk.IntVar(value=1000)
        self.as_weight = tk.DoubleVar()
        self.rms_rv1 = tk.DoubleVar()
        self.rms_rv2 = tk.DoubleVar()
        self.rms_as = tk.DoubleVar()
        self.do_weight = tk.BooleanVar()
        self.def_weight = tk.DoubleVar()

        # define labels and buttons in a grid
        tk.Label(min_frame, text='MINIMIZATION', font=('', titlesize, 'underline')).grid(columnspan=4)

        tk.Label(min_frame, text='Do MCMC?:').grid(row=1, sticky=tk.E)
        mcmc_check = tk.Checkbutton(min_frame, var=self.do_mcmc, command=self.toggle_mc)
        mcmc_check.grid(row=1, column=1, sticky=tk.W)
        self.steps_label = tk.Label(min_frame, text='# of samples:', state=tk.DISABLED)
        self.steps_label.grid(row=1, column=2, sticky=tk.E)
        self.steps_entry = tk.Entry(min_frame, textvariable=self.steps, width=5, state=tk.DISABLED)
        self.steps_entry.grid(row=1, column=3)
        tk.Label(min_frame, text='astrometric weight from data = ').grid(row=2, column=1, sticky=tk.E)
        self.def_weight_label = tk.Label(min_frame, textvariable=self.def_weight)
        self.def_weight_label.grid(row=2, column=2, sticky=tk.W)
        self.as_weight_button = tk.Checkbutton(min_frame, var=self.do_weight, command=self.toggle_weights)
        self.as_weight_button.grid(row=3, sticky=tk.E)
        self.weight_label = tk.Label(min_frame, text='Custom astrometric weight =', state=tk.DISABLED)
        self.weight_label.grid(row=3, column=1, sticky=tk.E)
        self.weight_slider = tk.Scale(min_frame, variable=self.as_weight, from_=0, to=1, orient=tk.HORIZONTAL,
                                      resolution=0.001, state=tk.DISABLED, length=180)
        self.weight_slider.grid(row=3, column=2, columnspan=2, sticky=tk.W)

        tk.Button(min_frame, text='Minimize!', command=self.minimize, highlightbackground=hcolor).grid(
            row=4, columnspan=4)
        tk.Label(min_frame, text='Results', font=('', titlesize, 'underline')).grid(row=5, columnspan=4)
        tk.Label(min_frame, text='Red. Chi Sqrd =').grid(row=6, sticky=tk.E)
        tk.Label(min_frame, text='Deg. of frdm =').grid(row=7, sticky=tk.E)
        tk.Label(min_frame, textvariable=self.redchisq).grid(row=6, column=1, sticky=tk.W)
        tk.Label(min_frame, textvariable=self.dof).grid(row=7, column=1, sticky=tk.W)
        tk.Label(min_frame, text='RMS Primary =').grid(row=6, column=2, sticky=tk.E)
        tk.Label(min_frame, text='RMS Secondary =').grid(row=7, column=2, sticky=tk.E)
        tk.Label(min_frame, text='RMS Rel. Orb. =').grid(row=8, column=2, sticky=tk.E)
        tk.Label(min_frame, textvariable=self.rms_rv1).grid(row=6, column=3, sticky=tk.W)
        tk.Label(min_frame, textvariable=self.rms_rv2).grid(row=7, column=3, sticky=tk.W)
        tk.Label(min_frame, textvariable=self.rms_as).grid(row=8, column=3, sticky=tk.W)
        self.mcplotbutton = tk.Button(min_frame, text='Make MCMC scatterplot matrix', command=self.plot_corner_diagram,
                                      highlightbackground=hcolor, state=tk.DISABLED)
        self.mcplotbutton.grid(row=9, columnspan=4)

        # PLOT CONTROLS
        tk.Label(plt_frame, text='PLOT CONTROLS', font=('', titlesize, 'underline')).grid(columnspan=6)

        self.do_phasedot = tk.BooleanVar()
        self.do_datarv1 = tk.BooleanVar()
        self.do_datarv2 = tk.BooleanVar()
        self.do_dataas = tk.BooleanVar()
        self.do_modelrv1 = tk.BooleanVar()
        self.do_modelrv2 = tk.BooleanVar()
        self.do_modelas = tk.BooleanVar()
        self.do_nodeline = tk.BooleanVar()
        self.do_semimajor = tk.BooleanVar()
        self.do_peri = tk.BooleanVar()
        self.do_as_dist = tk.BooleanVar()

        self.rv_plot_boolvars = [self.do_datarv1, self.do_datarv2, self.do_modelrv1, self.do_modelrv2]
        self.as_plot_boolvars = [self.do_dataas, self.do_modelas, self.do_nodeline, self.do_semimajor, self.do_peri]

        self.phase = tk.DoubleVar()
        self.phase.trace_add('write', lambda n, ix, m: self.update())

        self.phase_label = tk.Label(plt_frame, text='phase =', state=tk.DISABLED)
        self.phase_label.grid(row=1, column=1, sticky=tk.E)
        self.phase_slider = tk.Scale(plt_frame, variable=self.phase, from_=0, to=1, orient=tk.HORIZONTAL,
                                     resolution=0.005, length=300, state=tk.DISABLED)
        self.phase_slider.grid(row=1, column=2, columnspan=4)
        self.phase_button = tk.Checkbutton(plt_frame, var=self.do_phasedot, command=self.toggle_dots, state=tk.DISABLED)
        self.phase_button.grid(row=1)

        self.plot_rv1data_label = tk.Label(plt_frame, text='Primary RV data', state=tk.DISABLED)
        self.plot_rv1data_label.grid(row=2, column=1)
        self.plot_rv1data_button = tk.Checkbutton(plt_frame, var=self.do_datarv1,
                                                  command=lambda: self.update(self.do_datarv1), state=tk.DISABLED)
        self.plot_rv1data_button.grid(row=2)

        self.plot_rv2data_label = tk.Label(plt_frame, text='Secondary RV data', state=tk.DISABLED)
        self.plot_rv2data_label.grid(row=3, column=1)
        self.plot_rv2data_button = tk.Checkbutton(plt_frame, var=self.do_datarv2,
                                                  command=lambda: self.update(self.do_datarv2), state=tk.DISABLED)
        self.plot_rv2data_button.grid(row=3)

        self.plot_asdata_label = tk.Label(plt_frame, text='Astrometric data', state=tk.DISABLED)
        self.plot_asdata_label.grid(row=4, column=1)
        self.plot_asdata_button = tk.Checkbutton(plt_frame, var=self.do_dataas,
                                                 command=lambda: self.update(self.do_dataas), state=tk.DISABLED)
        self.plot_asdata_button.grid(row=4)

        self.plot_rv1model_label = tk.Label(plt_frame, text='Primary RV model', state=tk.DISABLED)
        self.plot_rv1model_label.grid(row=2, column=3)
        self.plot_rv1model_button = tk.Checkbutton(plt_frame, var=self.do_modelrv1,
                                                   command=lambda: self.update(self.do_modelrv1), state=tk.DISABLED)
        self.plot_rv1model_button.grid(row=2, column=2)

        self.plot_rv2model_label = tk.Label(plt_frame, text='Secondary RV model', state=tk.DISABLED)
        self.plot_rv2model_label.grid(row=3, column=3)
        self.plot_rv2model_button = tk.Checkbutton(plt_frame, var=self.do_modelrv2,
                                                   command=lambda: self.update(self.do_modelrv2), state=tk.DISABLED)
        self.plot_rv2model_button.grid(row=3, column=2)

        self.plot_asmodel_label = tk.Label(plt_frame, text='Model Orbit', state=tk.DISABLED)
        self.plot_asmodel_label.grid(row=4, column=3)
        self.plot_asmodel_button = tk.Checkbutton(plt_frame, var=self.do_modelas,
                                                  command=lambda: self.update(self.do_modelas), state=tk.DISABLED)
        self.plot_asmodel_button.grid(row=4, column=2)

        self.plot_nodeline_label = tk.Label(plt_frame, text='Line of nodes', state=tk.DISABLED)
        self.plot_nodeline_label.grid(row=2, column=5)
        self.plot_nodeline_button = tk.Checkbutton(plt_frame, var=self.do_nodeline,
                                                   command=lambda: self.update(self.do_nodeline), state=tk.DISABLED)
        self.plot_nodeline_button.grid(row=2, column=4)

        self.plot_semimajor_label = tk.Label(plt_frame, text='Semi-major axis', state=tk.DISABLED)
        self.plot_semimajor_label.grid(row=3, column=5)
        self.plot_semimajor_button = tk.Checkbutton(plt_frame, var=self.do_semimajor,
                                                    command=lambda: self.update(self.do_semimajor), state=tk.DISABLED)
        self.plot_semimajor_button.grid(row=3, column=4)

        self.plot_peri_label = tk.Label(plt_frame, text='Periastron', state=tk.DISABLED)
        self.plot_peri_label.grid(row=4, column=5)
        self.plot_peri_button = tk.Checkbutton(plt_frame, var=self.do_peri,
                                               command=lambda: self.update(self.do_peri), state=tk.DISABLED)
        self.plot_peri_button.grid(row=4, column=4)

        self.as_dist_label = tk.Label(plt_frame, text='Astrometric errors', state=tk.DISABLED)
        self.as_dist_label.grid(row=5, column=5)
        self.as_dist_button = tk.Checkbutton(plt_frame, var=self.do_as_dist,
                                             command=lambda: self.update(self.do_as_dist), state=tk.DISABLED)
        self.as_dist_button.grid(row=5, column=4)
        self.modelwidgets = {self.phase_label, self.phase_slider, self.plot_asmodel_label, self.plot_asmodel_button,
                             self.plot_rv1model_button, self.plot_rv2model_button, self.plot_rv1model_label,
                             self.plot_rv2model_label, self.plot_semimajor_button, self.plot_semimajor_label,
                             self.plot_nodeline_button, self.plot_nodeline_label, self.plot_peri_label,
                             self.plot_peri_button, self.phase_button, self.as_dist_button, self.as_dist_label}

        self.do_legend = tk.BooleanVar()
        legend_button = tk.Checkbutton(plt_frame, var=self.do_legend, highlightbackground=hcolor,
                                       command=lambda: self.update(self.do_legend))
        legend_button.grid(row=5)
        tk.Label(plt_frame, text='Legend').grid(row=5, column=1)

        self.plot_phase = tk.BooleanVar(value=True)
        self.plot_phase.trace_add('write', lambda n, ix, m: self.update())
        self.pphase_but = tk.Radiobutton(plt_frame, text='phase', variable=self.plot_phase, value=True)
        self.ptime_but = tk.Radiobutton(plt_frame, text='time', variable=self.plot_phase, value=False)
        self.pphase_but.grid(row=5, column=2)
        self.ptime_but.grid(row=5, column=3)

        # setup the plotting windows
        self.init_plots()

        # pack everything neatly
        data_frame.place(relx=.5, rely=0, anchor="n")
        tabs.add(data_frame_wrap, text='Data Files')

        guess_frame.grid(sticky=tk.N)
        infer_frame.grid(sticky=tk.N)
        guess_infer_frame.place(relx=.5, rely=0, anchor='n')
        tabs.add(guess_infer_wrap, text='System/Parameters')

        min_frame.place(relx=.5, rely=0, anchor="n")
        tabs.add(min_frame_wrap, text='Minimization')

        plt_frame.place(relx=.5, rely=0, anchor="n")
        tabs.add(plt_frame_wrap, text='Plot Controls')

        tabs.pack(fill=tk.BOTH, expand=True)

    def init_plots(self):
        """
        sets up the plot windows
        """

        def move_figure(f, x, y):
            """
            moves window f by x, y pixels
            :param f: window
            :param x: x offset
            :param y: y offset
            """
            f.canvas.manager.window.wm_geometry("+{}+{}".format(x, y))

        if self.rv_fig is not None:
            plt.close(self.rv_fig)
        if self.as_fig is not None:
            plt.close(self.as_fig)
        self.rv_fig = plt.figure(figsize=(10.5, 4.2
                                          ))
        self.as_fig = plt.figure(figsize=(10.5, 4.2))
        move_figure(self.rv_fig, int(0.35 * self._w) + 10, 0)
        move_figure(self.as_fig, int(0.35 * self._w) + 10, int(self._h / 2) + 10)
        self.rv_ax = self.rv_fig.add_subplot(111)
        self.as_ax = self.as_fig.add_subplot(111, aspect=1)
        spp.setup_rvax(self.rv_ax)
        spp.setup_asax(self.as_ax)
        self.rv_fig.tight_layout()
        self.as_fig.tight_layout()
        plt.ion()
        plt.show()

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

    def toggle_dots(self):
        """
        toogles the phase-dot widgets
        """
        for widg in self.phase_slider, self.phase_label:
            self.toggle(widg, self.do_phasedot.get())
        self.update()

    def toggle_mc(self):
        """
        toggles the MCMC widgets
        """
        for widg in self.steps_label, self.steps_entry:
            self.toggle(widg, self.do_mcmc.get())

    def toggle_rv1(self):
        """
        toggles the RV1 widgets
        """
        for widg in self.rv1_file, self.rv1_label, self.plot_rv1data_label, self.plot_rv1data_button:
            self.toggle(widg, self.include_rv1.get())
        if self.do_datarv1.get():
            self.do_datarv1.set(False)
        self.toggle_databutton()
        self.set_RV_or_AS_mode()
        self.update()

    def toggle_rv2(self):
        """
        toggles the RV2 widgets
        """
        if not self.include_rv1.get():
            self.include_rv2.set(False)
            return
        for widg in self.rv2_file, self.rv2_label, self.plot_rv2data_label, self.plot_rv2data_button:
            self.toggle(widg, self.include_rv2.get())
        if self.do_datarv2.get():
            self.do_datarv2.set(False)
        self.toggle_databutton()
        self.set_RV_or_AS_mode()
        self.update()

    def toggle_as(self):
        """
        toggles the AS widgets
        """
        for widg in (self.as_file, self.as_label, self.seppa_but, self.en_but, self.plot_asdata_label,
                     self.plot_asdata_button):
            self.toggle(widg, self.include_as.get())
        if self.do_dataas.get():
            self.do_dataas.set(False)
        self.toggle_databutton()
        self.set_RV_or_AS_mode()
        self.update()

    def toggle_databutton(self):
        """
        toggles the load data button
        """
        if not (self.include_rv1.get() or self.include_rv2.get() or self.include_as.get()):
            self.toggle(self.data_button, False)
        else:
            self.toggle(self.data_button, True)

    def toggle_weights(self):
        """
        toggles the weight widgets
        """
        if not (self.do_weight.get()) or not (
                self.include_as.get() and (self.include_rv1.get() or self.include_rv2.get())):
            self.toggle(self.weight_label, False)
            self.toggle(self.weight_slider, False)
            self.do_weight.set(False)
        else:
            self.toggle(self.weight_label, True)
            self.toggle(self.weight_slider, True)

    def set_RV_or_AS_mode(self):
        """
        sets the parameters in the correct inference mode
        """
        for lst in self.param_labels, self.vary_button_list:
            if self.include_rv1.get() and self.include_rv2.get() and self.include_as.get():
                for i in {2, 4, 6, 7, 8, 9, 10, 11}:
                    lst[i].config(state=tk.NORMAL)
            elif self.include_rv1.get() and self.include_rv2.get():
                for i in {7, 8, 9, 10}:
                    lst[i].config(state=tk.NORMAL)
                for i in {2, 4, 6, 11}:
                    lst[i].config(state=tk.DISABLED)
            elif self.include_rv1.get() and self.include_as.get():
                for i in {2, 4, 6, 7, 9, 11}:
                    lst[i].config(state=tk.NORMAL)
                for i in {8, 10}:
                    lst[i].config(state=tk.DISABLED)
            elif self.include_rv1.get():
                for i in {7, 9, 11}:
                    lst[i].config(state=tk.NORMAL)
                for i in {2, 4, 6, 8, 10, 11}:
                    lst[i].config(state=tk.DISABLED)
            elif self.include_as.get():
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
            self.loading_guesses = True
            self.guess_dict = spl.guess_loader(self.wd.get(), self.guess_file.get())
        except IOError:
            print('cannot find your guess file!')
            self.loading_guesses = False
            self.guess_dict = None
            return
        try:
            self.guess_var_list[0].set(self.guess_dict['p'][0])
            self.guess_var_list[1].set(self.guess_dict['e'][0])
            self.guess_var_list[2].set(self.guess_dict['i'][0])
            self.guess_var_list[3].set(self.guess_dict['omega'][0])
            self.guess_var_list[4].set(self.guess_dict['Omega'][0])
            self.guess_var_list[5].set(self.guess_dict['t0'][0])
            self.guess_var_list[6].set(self.guess_dict['d'][0])
            self.guess_var_list[7].set(self.guess_dict['k1'][0])
            self.guess_var_list[8].set(self.guess_dict['k2'][0])
            self.guess_var_list[9].set(self.guess_dict['gamma1'][0])
            self.guess_var_list[10].set(self.guess_dict['gamma2'][0])
            self.guess_var_list[11].set(self.guess_dict['mt'][0])
            self.vary_var_list[0].set(str(self.guess_dict['p'][1]))
            self.vary_var_list[1].set(str(self.guess_dict['e'][1]))
            self.vary_var_list[2].set(str(self.guess_dict['i'][1]))
            self.vary_var_list[3].set(str(self.guess_dict['omega'][1]))
            self.vary_var_list[4].set(str(self.guess_dict['Omega'][1]))
            self.vary_var_list[5].set(str(self.guess_dict['t0'][1]))
            self.vary_var_list[6].set(str(self.guess_dict['d'][1]))
            self.vary_var_list[7].set(str(self.guess_dict['k1'][1]))
            self.vary_var_list[8].set(str(self.guess_dict['k2'][1]))
            self.vary_var_list[9].set(str(self.guess_dict['gamma1'][1]))
            self.vary_var_list[10].set(str(self.guess_dict['gamma2'][1]))
            self.vary_var_list[11].set(str(self.guess_dict['mt'][1]))
        except (ValueError, KeyError, TypeError):
            print('some parameter has not been set properly')
            self.loading_guesses = False
            self.guess_dict = None
            return
        self.loading_guesses = False
        self.update()

    def set_guess_dict_from_entries(self):
        """
        builds the guess dict from the guess column
        """
        self.guess_dict = dict()
        self.guess_dict = {'p': (float(self.guess_var_list[0].get()), self.vary_var_list[0].get()),
                           'e': (float(self.guess_var_list[1].get()), self.vary_var_list[1].get()),
                           'i': (float(self.guess_var_list[2].get()), self.vary_var_list[2].get()),
                           'omega': (float(self.guess_var_list[3].get()), self.vary_var_list[3].get()),
                           'Omega': (float(self.guess_var_list[4].get()), self.vary_var_list[4].get()),
                           't0': (float(self.guess_var_list[5].get()), self.vary_var_list[5].get()),
                           'd': (float(self.guess_var_list[6].get()), self.vary_var_list[6].get()),
                           'k1': (float(self.guess_var_list[7].get()), self.vary_var_list[7].get()),
                           'k2': (float(self.guess_var_list[8].get()), self.vary_var_list[8].get()),
                           'gamma1': (float(self.guess_var_list[9].get()), self.vary_var_list[9].get()),
                           'gamma2': (float(self.guess_var_list[10].get()), self.vary_var_list[10].get()),
                           'mt': (float(self.guess_var_list[11].get()), self.vary_var_list[11].get())}
        self.param_dict = dict()
        for param, value in self.guess_dict.items():
            self.param_dict[param] = value[0]

    def set_system(self):
        """
        sets the system from the current guess column
        """
        if self.loading_guesses:
            return False
        try:
            self.set_guess_dict_from_entries()
            self.system = bsys.System(self.param_dict)
        except ValueError:
            print('invalid model!')
            self.guess_dict = None
            self.system = None
            for widg in self.modelwidgets:
                self.toggle(widg, False)
            return False
        else:
            for widg in self.modelwidgets:
                self.toggle(widg, True)
            self.mprimary.set(np.round(self.system.primary_mass(), 2))
            self.msecondary.set(np.round(self.system.secondary_mass(), 2))
            self.totalmass.set(np.round(self.system.total_mass(), 2))
            self.semimajork1k2.set(np.round(self.system.semimajor_axis_from_RV(), 2))
            self.semimajord.set(np.round(self.system.semimajor_axis_from_distance(), 2))
            return True

    def load_data(self):
        """
        loads the data from the current selected files
        """
        filetypes = list()
        filenames = list()
        if self.rv1_file.get() != '' and self.include_rv1.get():
            filetypes.append('RV1file')
            filenames.append(self.rv1_file.get())
        if self.rv2_file.get() != '' and self.include_rv2.get():
            filetypes.append('RV2file')
            filenames.append(self.rv2_file.get())
        if self.as_file.get() != '' and self.include_as.get():
            filetypes.append('ASfile')
            filenames.append(self.as_file.get())
        try:
            self.data_dict = dict()
            self.data_dict = spl.data_loader(self.wd.get(), filetypes, filenames, self.seppa.get())
        except (OSError, KeyError) as e:
            self.data_dict = None
            return False
        if not self.include_as.get():
            self.def_weight.set(0)
        elif not self.include_rv1.get() and not self.include_rv2.get():
            self.def_weight.set(1)
        else:
            if self.include_rv1.get() and self.include_rv2.get():
                w = 2 * len(self.data_dict['AS']['hjds']) / (
                        2 * len(self.data_dict['AS']['hjds']) + len(self.data_dict['RV1']['hjds']) +
                        len(self.data_dict['RV2']['hjds']))
            elif self.include_rv1.get():
                w = 2 * len(self.data_dict['AS']['hjds']) / (
                        2 * len(self.data_dict['AS']['hjds']) + len(self.data_dict['RV1']['hjds']))
            else:
                w = 2 * len(self.data_dict['AS']['hjds']) / (
                        2 * len(self.data_dict['AS']['hjds']) + len(self.data_dict['RV2']['hjds']))
            self.def_weight.set(np.round(w, 3))

    def minimize(self):
        """
        launches a minimization run
        """
        self.set_guess_dict_from_entries()
        self.load_data()
        if self.guess_dict is not None and self.data_dict is not None:
            # calculate best parameters
            try:
                if self.do_weight.get():
                    w = self.as_weight.get()
                else:
                    w = None
                self.minresult, rms_rv1, rms_rv2, rms_as = \
                    spm.LMminimizer(self.guess_dict, self.data_dict, self.do_mcmc.get(), self.steps.get(), w)
                if self.do_mcmc.get():
                    self.didmcmc = True
                    self.mcmc_run_number += 1
                    self.toggle(self.mcplotbutton, True)
                else:
                    self.didmcmc = False
                pars = self.minresult.params
                # fill in the entries
                if pars['p'].vary:
                    self.mininimzed_var_list[0].set(np.round(pars['p'].value, 3))
                    self.error_var_list[0].set(np.round(pars['p'].stderr, 3))
                if pars['e'].vary:
                    self.mininimzed_var_list[1].set(np.round(pars['e'].value, 3))
                    self.error_var_list[1].set(np.round(pars['e'].stderr, 3))
                if pars['i'].vary:
                    self.mininimzed_var_list[2].set(np.round(pars['i'].value % 360, 3))
                    self.error_var_list[2].set(np.round(pars['i'].stderr, 3))
                if pars['omega'].vary:
                    self.mininimzed_var_list[3].set(np.round(pars['omega'].value % 360, 3))
                    self.error_var_list[3].set(np.round(pars['omega'].stderr, 3))
                if pars['Omega'].vary:
                    self.mininimzed_var_list[4].set(np.round(pars['Omega'].value % 360, 3))
                    self.error_var_list[4].set(np.round(pars['Omega'].stderr, 3))
                if pars['t0'].vary:
                    self.mininimzed_var_list[5].set(np.round(pars['t0'].value, 3))
                    self.error_var_list[5].set(np.round(pars['t0'].stderr, 3))
                if pars['d'].vary:
                    self.mininimzed_var_list[6].set(np.round(pars['d'].value, 3))
                    self.error_var_list[6].set(np.round(pars['d'].stderr, 3))
                if pars['k1'].vary:
                    self.mininimzed_var_list[7].set(np.round(pars['k1'].value, 3))
                    self.error_var_list[7].set(np.round(pars['k1'].stderr, 3))
                if pars['k2'].vary:
                    self.mininimzed_var_list[8].set(np.round(pars['k2'].value, 3))
                    self.error_var_list[8].set(np.round(pars['k2'].stderr, 3))
                if pars['gamma1'].vary:
                    self.mininimzed_var_list[9].set(np.round(pars['gamma1'].value, 3))
                    self.error_var_list[9].set(np.round(pars['gamma1'].stderr, 3))
                if pars['gamma2'].vary:
                    self.mininimzed_var_list[10].set(np.round(pars['gamma2'].value, 3))
                    self.error_var_list[10].set(np.round(pars['gamma2'].stderr, 3))
                if pars['mt'].vary:
                    self.mininimzed_var_list[11].set(np.round(pars['mt'].value, 3))
                    self.error_var_list[11].set(np.round(pars['mt'].stderr, 3))

                self.redchisq.set(np.round(self.minresult.redchi, 4))
                self.dof.set(self.minresult.nfree)
                self.rms_rv1.set(np.round(rms_rv1, 4))
                self.rms_rv2.set(np.round(rms_rv2, 4))
                self.rms_as.set(np.round(rms_as, 4))
                self.minimization_run_number += 1
            except ValueError as e:
                print(e)

    def update(self, plot_bool=None):
        """
        updates the gui, by replotting everything that is selected
        :param plot_bool: plotting button that was clicked (optional)
        """
        if time.time() - self.since < 1.5:
            return
        else:
            self.since = time.time()
        if self.loading_guesses or not self.set_system():
            if plot_bool:
                plot_bool.set(not plot_bool.get())
            return

        self.load_data()

        # cannot find a way to condense this without messing up references to line objects
        if self.do_dataas.get():
            self.plot_as_data()
        else:
            if self.asdata_line:
                self.asdata_line.remove()
                self.asdata_line = None
            if self.as_ellipses:
                self.as_ellipses.remove()
                self.as_ellipses = None

        if self.do_phasedot.get() and self.plot_phase.get():
            self.plot_dots()
        else:
            if self.as_dot:
                self.as_dot.remove()
                self.as_dot = None
            if self.rv1_dot:
                self.rv1_dot.remove()
                self.rv1_dot = None
            if self.rv2_dot:
                self.rv2_dot.remove()
                self.rv2_dot = None
        if self.do_peri.get():
            self.plot_periastron()
        else:
            if self.peri_dot:
                self.peri_dot.remove()
                self.peri_dot = None
        if self.do_semimajor.get():
            self.plot_semimajor_axis()
        else:
            if self.semi_major:
                self.semi_major.remove()
                self.semi_major = None
        if self.do_nodeline.get():
            self.plot_node_line()
        else:
            if self.node_line:
                self.node_line.remove()
                self.node_line = None
        if self.do_modelas.get():
            self.plot_relative_orbit()
        else:
            if self.as_line:
                self.as_line.remove()
                self.as_line = None
        if self.do_modelrv2.get():
            self.plot_rv2_curve()
        else:
            if self.rv2_line:
                self.rv2_line.remove()
                self.rv2_line = None
        if self.do_modelrv1.get():
            self.plot_rv1_curve()
        else:
            if self.rv1_line:
                self.rv1_line.remove()
                self.rv1_line = None
        if self.do_datarv2.get():
            self.plot_rv2_data()
        else:
            if self.rv2data_line:
                self.rv2data_line.remove()
                self.rv2data_line = None
        if self.do_datarv1.get():
            self.plot_rv1_data()
        else:
            if self.rv1data_line:
                self.rv1data_line.remove()
                self.rv1data_line = None
        if self.do_as_dist.get():
            self.plot_as_dist()
        else:
            if self.as_dist_lines:
                for line in self.as_dist_lines:
                    line.remove()
                self.as_dist_lines = None
        self.plot_legends()
        self.relim_plots()
        self.rv_fig.canvas.draw()
        self.as_fig.canvas.draw()

    def relim_plots(self):
        """
        resizes the plots according to the data limits
        """
        for plot_bool in self.rv_plot_boolvars:
            if plot_bool.get():
                self.rv_ax.relim()
                self.rv_ax.axis('auto')

        for plot_bool in self.as_plot_boolvars:
            if plot_bool.get():
                self.as_ax.relim()
                self.as_ax.axis('image')

    def plot_rv1_data(self):
        """
        plot the rv1 data
        """
        if 'RV1' not in self.data_dict:
            return
        if self.rv1data_line is not None:
            self.rv1data_line.remove()
            self.rv1data_line = None
        if self.plot_phase.get():
            phases, rv, err = self.system.create_phase_extended_RV(self.data_dict['RV1'], 0.15)
            self.rv1data_line = self.rv_ax.errorbar(phases, rv, yerr=err, ls='', capsize=0.1, marker='o',
                                                    ms=5, color='b')
        else:
            self.rv1data_line = self.rv_ax.errorbar(self.data_dict['RV1']['hjds'], self.data_dict['RV1']['RVs'],
                                                    yerr=self.data_dict['RV1']['errors'], ls='', capsize=0.1,
                                                    marker='o', ms=5, color='b')

    def plot_rv2_data(self):
        """
        plot the rv2 data
        """
        if 'RV2' not in self.data_dict:
            return
        if self.rv2data_line is not None:
            self.rv2data_line.remove()
            self.rv2data_line = None
        if self.plot_phase.get():
            phases, rv, err = self.system.create_phase_extended_RV(self.data_dict['RV2'], 0.15)
            self.rv2data_line = self.rv_ax.errorbar(phases, rv, yerr=err, ls='', capsize=0.1, marker='o',
                                                    ms=5, color='r')
        else:
            self.rv2data_line = self.rv_ax.errorbar(self.data_dict['RV2']['hjds'], self.data_dict['RV2']['RVs'],
                                                    yerr=self.data_dict['RV2']['errors'], ls='', capsize=0.1,
                                                    marker='o', ms=5, color='r')

    def plot_as_data(self):
        """
        plot the as data
        """
        if 'AS' not in self.data_dict:
            return
        data = self.data_dict['AS']
        if self.asdata_line is None:
            self.asdata_line, = self.as_ax.plot(data['easts'], data['norths'], 'r.', ls='',
                                                label='Relative position')
        else:
            self.asdata_line.set_xdata(data['easts'])
            self.asdata_line.set_ydata(data['norths'])
        if self.as_ellipses is not None:
            self.as_ellipses.remove()
        self.as_ellipses = EllipseCollection(2 * data['majors'], 2 * data['minors'], data['pas'] - 90,
                                             offsets=np.column_stack((data['easts'], data['norths'])),
                                             transOffset=self.as_ax.transData,
                                             units='x', edgecolors='r', facecolors=(0, 0, 0, 0))
        self.as_ax.add_collection(self.as_ellipses)

    def plot_as_dist(self):
        """
        plot the astrometric distances of each as point
        """
        if 'AS' not in self.data_dict:
            return
        data = self.data_dict['AS']
        if self.as_dist_lines is not None:
            for line in self.as_dist_lines:
                line.remove()
            self.as_dist_lines = None
        self.as_dist_lines = list()
        for i in range(len(data['hjds'])):
            self.as_dist_lines.append(self.as_ax.plot(
                (data['easts'][i], self.system.relative.east_of_hjd(data['hjds'][i])),
                (data['norths'][i], self.system.relative.north_of_hjd(data['hjds'][i])), 'k')[0])

    def plot_rv1_curve(self):
        """
        the the rv1 model curve
        """
        if self.plot_phase.get():
            phases = np.linspace(-0.15, 1.15, num=150)
            vrads1 = self.system.primary.radial_velocity_of_phases(phases)
            if self.rv1_line is None:
                self.rv1_line, = self.rv_ax.plot(phases, vrads1, label=r'primary', color='b', ls='--')
            else:
                self.rv1_line.set_xdata(phases)
                self.rv1_line.set_ydata(vrads1)
            self.rv_ax.set_xlabel(r'orbital phase')
        else:
            m = min(self.data_dict['RV1']['hjds'])
            mm = max(self.data_dict['RV1']['hjds'])
            try:
                m = min(m, min(self.data_dict['RV2']['hjds']))
                mm = max(mm, max(self.data_dict['RV2']['hjds']))
            except KeyError:
                pass
            times = np.linspace(m - 0.01 * (mm - m), m - 0.01 * (mm - m) + self.system.p, endpoint=False, num=100)
            rvs = self.system.primary.radial_velocity_of_phases(self.system.phase_of_hjds(times))
            times, rvs = self.extend_rvs_over_time(times, rvs, mm)
            if self.rv1_line is None:
                self.rv1_line, = self.rv_ax.plot(times, rvs, label=r'primary', color='b', ls='--')
            else:
                self.rv1_line.set_xdata(times)
                self.rv1_line.set_ydata(rvs)
            self.rv_ax.set_xlabel(r'time (JD)')

    def extend_rvs_over_time(self, period, rvs, mm):
        times = np.copy(period)
        n = 0
        while period[-1] <= mm:
            n += 1
            period += self.system.p
            times = np.concatenate((times, period))
        rvs = np.tile(rvs, n + 1)
        return times, rvs

    def plot_rv2_curve(self):
        """
        plot the rv2 model curve
        """
        if self.plot_phase.get():
            phases = np.linspace(-0.15, 1.15, num=150)
            vrads1 = self.system.secondary.radial_velocity_of_phases(phases)
            if self.rv2_line is None:
                self.rv2_line, = self.rv_ax.plot(phases, vrads1, label=r'secondary', color='r', ls='--')
            else:
                self.rv2_line.set_xdata(phases)
                self.rv2_line.set_ydata(vrads1)
            self.rv_ax.set_xlabel(r'orbital phase')
        else:
            m = min(self.data_dict['RV2']['hjds'])
            mm = max(self.data_dict['RV2']['hjds'])
            try:
                m = min(m, min(self.data_dict['RV1']['hjds']))
                mm = max(mm, max(self.data_dict['RV1']['hjds']))
            except KeyError:
                pass
            times = np.linspace(self.system.t0, self.system.t0 + self.system.p, num=100)
            rvs = self.system.secondary.radial_velocity_of_phases(self.system.phase_of_hjds(times))
            times, rvs = self.extend_rvs_over_time(times, rvs, mm)
            if self.rv2_line is None:
                self.rv2_line, = self.rv_ax.plot(times, rvs, label=r'secondary', color='r', ls='--')
            else:
                self.rv2_line.set_xdata(times)
                self.rv2_line.set_ydata(rvs)
            self.rv_ax.set_xlabel(r'time (JD)')

    def plot_relative_orbit(self):
        """
        plot the relative astrometric orbit
        """
        ecc_anoms = np.linspace(0, 2 * np.pi, 200)
        norths = self.system.relative.north_of_ecc(ecc_anoms)
        easts = self.system.relative.east_of_ecc(ecc_anoms)
        if self.as_line is None:
            self.as_line, = self.as_ax.plot(easts, norths, label='relative orbit', color='k')
        else:
            self.as_line.set_xdata(easts)
            self.as_line.set_ydata(norths)

    def plot_node_line(self):
        """
        plot the astrometric node line
        """
        system = self.system.relative
        if self.node_line is None:
            self.node_line, = self.as_ax.plot([system.east_of_true(-system.omega),
                                               system.east_of_true(-system.omega + np.pi)],
                                              [system.north_of_true(-system.omega),
                                               system.north_of_true(-system.omega + np.pi)],
                                              color='0.5', ls='--', label='line of nodes')
        else:
            self.node_line.set_xdata([system.east_of_true(-system.omega),
                                      system.east_of_true(-system.omega + np.pi)])
            self.node_line.set_ydata([system.north_of_true(-system.omega),
                                      system.north_of_true(-system.omega + np.pi)])

    def plot_periastron(self):
        """
        plot the astrometric periastron point
        """
        system = self.system.relative
        if self.peri_dot is None:
            self.peri_dot, = self.as_ax.plot([system.east_of_ecc(0)], [system.north_of_ecc(0)], color='b', marker='s',
                                             ls='', fillstyle='full', label='periastron', markersize=8)
        else:
            self.peri_dot.set_xdata(system.east_of_ecc(0))
            self.peri_dot.set_ydata(system.north_of_ecc(0))

    def plot_semimajor_axis(self):
        """
        plot the astrometric semimajor axis
        """
        system = self.system.relative
        if self.semi_major is None:
            self.semi_major, = self.as_ax.plot([system.east_of_true(0), system.east_of_true(np.pi)],
                                               [system.north_of_true(0), system.north_of_true(np.pi)],
                                               color='0.3', ls='dotted', label='semi-major axis')
        else:
            self.semi_major.set_xdata([system.east_of_true(0), system.east_of_true(np.pi)])
            self.semi_major.set_ydata([system.north_of_true(0), system.north_of_true(np.pi)])

    def plot_dots(self):
        """
        plot the phase dots
        """
        if self.rv1_dot is not None:
            self.rv1_dot.remove()
            self.rv1_dot = None
        if self.do_modelrv1.get() or self.do_datarv1.get():
            rv1 = self.system.primary.radial_velocity_of_phase(self.phase.get())
            self.rv1_dot = self.rv_ax.scatter(self.phase.get(), rv1, s=100, color='b', marker='D',
                                              label=np.round(rv1, 2))
        if self.rv2_dot is not None:
            self.rv2_dot.remove()
            self.rv2_dot = None
        if self.do_modelrv2.get() or self.do_datarv2.get():
            rv2 = self.system.secondary.radial_velocity_of_phase(self.phase.get())
            self.rv2_dot = self.rv_ax.scatter(self.phase.get(), rv2, s=100, color='r', marker='D',
                                              label=np.round(rv2, 2))
        if self.as_dot is not None:
            self.as_dot.remove()
            self.as_dot = None
        if self.do_modelas.get() or self.do_dataas.get():
            N = self.system.relative.north_of_ph(self.phase.get())
            E = self.system.relative.east_of_ph(self.phase.get())
            self.as_dot = self.as_ax.scatter(E, N, s=100, color='r', marker='x',
                                             label='{}E/{}N'.format(np.round(E, 2), np.round(N, 2)))

    def plot_legends(self):
        """
        plot the legends
        """
        try:
            self.rv_ax.get_legend().remove()
        except AttributeError:
            pass
        try:
            self.as_ax.get_legend().remove()
        except AttributeError:
            pass
        if self.do_legend.get():
            if len(self.rv_ax.get_lines()) > 1:
                self.rv_ax.legend()
            if len(self.as_ax.get_lines()) > 1:
                self.as_ax.legend()

    def plot_corner_diagram(self):
        """
        plot a corner diagram of an MCMC run
        """
        if self.didmcmc:
            corner = spp.plot_corner_diagram(self.minresult)
            corner.savefig(self.wd.get() + 'corner{}.png'.format(self.mcmc_run_number))
            plt.close(corner)
        else:
            print('do an mcmc minimization first!')

    def save_params(self):
        """
        save minimized parameters to a file
        """
        with open(self.wd.get() + 'params_run{}.txt'.format(self.minimization_run_number), 'w') as f:
            for i in range(len(self.mininimzed_var_list)):
                f.write(str(self.param_name_vars[i].get()) + ' ' + str(self.mininimzed_var_list[i].get()) + ' ' + str(
                    self.error_var_list[i].get()) + '\n')
            f.write('reduced chisq = {} \n'.format(self.redchisq.get()))
            f.write('dof = {} \n'.format(self.dof.get()))
            f.write('param order: {}'.format(self.minresult.var_names))
        np.savetxt(self.wd.get() + '/covar{}.txt'.format(self.minimization_run_number), self.minresult.covar)

    def save_guesses(self):
        """
        save guesses parameters to a file
        """
        self.set_guess_dict_from_entries()
        spl.guess_saver(self.wd.get(), self.guess_dict)


def run(wd):
    """
    Run the main spinOS app tk loop.
    :param wd: working directory to set as root.
    """
    root = tk.Tk()
    wdir = pathlib.PurePath(__file__).parent.parent
    w, h = root.winfo_screenwidth(), root.winfo_screenheight()
    with splash.Splash(root, wdir.joinpath('rsc/spinos100.png'), 2.1, w, h):
        root.geometry("{}x{}+0+0".format(int(0.35 * w), h))
        root.title('spinOSgui')
        SpinOSApp(root, wd, w, h)

    root.mainloop()
