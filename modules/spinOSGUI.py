"""
module that contains the main spinOSApp class and tk loop
"""
import pathlib
import tkinter as tk
from tkinter import ttk

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import EllipseCollection

import spinOSio as spl
import spinOSminimizer as spm
import spinOSplotter as spp
import binary_system as bsys
import spinOSsplash as splash
import constants as cst

TIME_STR = r'time [day]'
PHASE_STR = r'orbital phase'


class SpinOSGUI:
    """
    class specifying the main spinOS tk implementation
    """

    def __init__(self, master, wwd, width, height):

        mpl.use("TkAgg")  # set the backend

        # FRAME STRUCTURE #
        # set the root frame
        tabs = ttk.Notebook(master)
        # set the data frame
        data_frame_tab = tk.Frame(tabs)
        data_frame = tk.Frame(data_frame_tab)
        # set the guess frame
        guess_infer_tab = tk.Frame(tabs)
        guess_infer_top = tk.Frame(guess_infer_tab)
        # set the inferation frame
        infer_frame = tk.Frame(guess_infer_top)
        guess_frame = tk.Frame(guess_infer_top)
        min_frame_tab = tk.Frame(tabs)
        min_frame = tk.Frame(min_frame_tab)
        # set the plot window controls frame
        plt_frame_tab = tk.Frame(tabs)
        plt_frame_top = tk.Frame(plt_frame_tab)
        plt_frame = tk.Frame(plt_frame_top)
        out_frame_tab = tk.Frame(tabs)
        out_frame = tk.Frame(out_frame_tab)

        # VARS #
        self._w, self._h = width, height
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

        # DATA FRAME #
        firstlabel = tk.Label(data_frame, text='DATA', font=('', cst.TITLESIZE, 'underline'))
        firstlabel.grid(columnspan=5, sticky=tk.N)

        # define inlcusion variables
        self.include_rv1 = tk.BooleanVar()
        self.include_rv2 = tk.BooleanVar()
        self.include_as = tk.BooleanVar()

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
            self.wd.insert(0, wwd)
        self.rv1_file.insert(0, 'primary_vels.txt')
        self.rv2_file.insert(0, 'secondary_vels.txt')
        self.as_file.insert(0, 'relative_astrometry.txt')
        self.guess_file.insert(0, 'guesses.txt')

        # disable them, needed after inserting stuff
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
        self.seppa.trace_add('write', lambda n, ix, m: self.load_data())
        self.seppa_but = tk.Radiobutton(data_frame, text='Sep/PA', variable=self.seppa, value=True, state=tk.DISABLED)
        self.seppa_but.grid(row=4, column=3)
        self.en_but = tk.Radiobutton(data_frame, text='E/N', variable=self.seppa, value=False, state=tk.DISABLED)
        self.en_but.grid(row=4, column=4)

        data_frame.place(relx=.5, rely=0, anchor="n")

        # GUESS FRAME #
        columns = 7
        labelcolumn = 1
        paramcolumn = 2
        varycheckcolumn = 3
        transfercolumn = 4
        minresultcolumn = 5
        errorcolumn = 6

        numofparams = 12
        rparams = range(numofparams)

        # print the labels in the guess frame
        tk.Label(guess_frame, text='System PARAMETERS', font=('', cst.TITLESIZE, 'underline')).grid(columnspan=columns)
        tk.Label(guess_frame, text='Guesses').grid(row=1, column=paramcolumn)
        tk.Label(guess_frame, text='Vary?').grid(row=1, column=varycheckcolumn)
        tk.Label(guess_frame, text='Transfer').grid(row=1, column=transfercolumn)
        tk.Label(guess_frame, text='Result').grid(row=1, column=minresultcolumn)
        tk.Label(guess_frame, text='Error').grid(row=1, column=errorcolumn)

        self.lock_gs = tk.BooleanVar(False)
        self.locked_image = tk.PhotoImage(file="rsc/locked.png")
        self.unlocked_image = tk.PhotoImage(file="rsc/unlocked.png")
        self.lock_gs_button = tk.Button(guess_frame, image=self.locked_image, command=self.toggle_lock)
        self.lock_gs_button.grid(row=12)

        self.lock_q = tk.BooleanVar(False)
        self.lock_q_button = tk.Button(guess_frame, width=1, text='q', command=self.toggle_q)
        self.lock_q_button.grid(row=10)

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

        self.param_label_list = [tk.Label(guess_frame, textvariable=self.param_var_list[i]) for i in
                                 rparams]

        for i in rparams:
            self.param_label_list[i].grid(row=(i + 2), column=labelcolumn, sticky=tk.E)

        # initialize the entry variables
        self.guess_var_list = [tk.StringVar(value='0') for _ in rparams]

        # define entry boxes
        self.guess_entry_list = [tk.Entry(guess_frame, textvariable=self.guess_var_list[i], width=10) for i in
                                 rparams]
        # put in a nice grid
        for i in rparams:
            self.guess_entry_list[i].grid(row=(i + 2), column=paramcolumn)

        # define the vary state variables
        self.vary_var_list = [tk.BooleanVar() for _ in rparams]

        # define checkbuttons for vary states
        self.vary_button_list = [tk.Checkbutton(guess_frame, var=self.vary_var_list[i]) for i in
                                 rparams]

        # put the checkbuttons in a nice grid
        for i in rparams:
            self.vary_button_list[i].grid(row=(i + 2), column=varycheckcolumn)

        # define the transfer buttons
        # for this semantic to work, we need to wrap the lambda function into another one, so that each command
        # references to its own number 'y', rather than the outer 'i' of the list comprehension
        self.transfer_button_list = [
            tk.Button(guess_frame, text='<-', command=(lambda y: (lambda: self.transfer(y)))(i)
                      ).grid(row=(i + 2), column=transfercolumn)
            for i in rparams]

        # define the minimized parameter variables
        self.mininimzed_var_list = [tk.StringVar() for _ in rparams]

        # define the labels the minimized parameters will go in
        self.min_label_list = [tk.Label(guess_frame, textvariable=self.mininimzed_var_list[i], width=8) for i in
                               rparams]

        for i in rparams:
            self.min_label_list[i].grid(row=(i + 2), column=minresultcolumn)

        # define the error variables
        self.error_var_list = [tk.StringVar() for _ in rparams]

        # define the labels the errors will go in
        self.error_label_list = [tk.Label(guess_frame, textvariable=self.error_var_list[i], width=8) for i in
                                 rparams]
        for i in rparams:
            self.error_label_list[i].grid(row=(i + 2), column=errorcolumn)

        # define the buttons in this frame
        tk.Button(guess_frame, text='Load guesses', command=self.load_guesses,
                  highlightbackground=cst.HCOLOR).grid(row=numofparams + 2, column=labelcolumn)
        tk.Button(guess_frame, text='Save guesses', command=self.save_guesses,
                  highlightbackground=cst.HCOLOR).grid(row=numofparams + 2, column=paramcolumn)
        tk.Button(guess_frame, text='Save parameters', command=self.save_params, highlightbackground=cst.HCOLOR).grid(
            row=numofparams + 2, column=minresultcolumn, columnspan=2)

        refreshframe1 = tk.Frame(guess_infer_top)
        tk.Button(refreshframe1, text='Refresh Plots & Inferred Parameters', width=30, height=2, command=self.update,
                  highlightbackground=cst.HCOLOR).pack()

        # INFER FRAME #
        self.mprimary = tk.StringVar()
        self.msecondary = tk.StringVar()
        self.semimajord = tk.StringVar()
        self.semimajork1k2 = tk.StringVar()
        self.totalmass = tk.StringVar()

        # define labels
        tk.Label(infer_frame, text='INFERRED PARAMETERS (from guesses)', font=('', cst.TITLESIZE, 'underline')) \
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

        guess_frame.grid(sticky=tk.N)
        refreshframe1.grid(sticky=tk.N, pady=10)
        infer_frame.grid(sticky=tk.N)
        guess_infer_top.place(relx=1, rely=0, anchor='ne')

        # MINIMIZATION FRAME #
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
        tk.Label(min_frame, text='MINIMIZATION', font=('', cst.TITLESIZE, 'underline')).grid(columnspan=4)
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

        tk.Button(min_frame, text='Minimize!', command=self.minimize, highlightbackground=cst.HCOLOR).grid(
            row=4, columnspan=4)
        tk.Label(min_frame, text='Results', font=('', cst.TITLESIZE, 'underline')).grid(row=5, columnspan=4)
        tk.Label(min_frame, text='Red. Chi Sqrd =').grid(row=6, sticky=tk.E)
        tk.Label(min_frame, text='Deg. of frdm =').grid(row=7, sticky=tk.E)
        tk.Label(min_frame, textvariable=self.redchisq).grid(row=6, column=1, sticky=tk.W)
        tk.Label(min_frame, textvariable=self.dof).grid(row=7, column=1, sticky=tk.W)
        tk.Label(min_frame, text='RMS Primary (km/s) =').grid(row=6, column=2, sticky=tk.E)
        tk.Label(min_frame, text='RMS Secondary (km/s) =').grid(row=7, column=2, sticky=tk.E)
        tk.Label(min_frame, text='RMS Rel. Orbit (mas) =').grid(row=8, column=2, sticky=tk.E)
        tk.Label(min_frame, textvariable=self.rms_rv1).grid(row=6, column=3, sticky=tk.W)
        tk.Label(min_frame, textvariable=self.rms_rv2).grid(row=7, column=3, sticky=tk.W)
        tk.Label(min_frame, textvariable=self.rms_as).grid(row=8, column=3, sticky=tk.W)
        self.mcplotbutton = tk.Button(min_frame, text='Make MCMC scatterplot matrix', command=self.make_corner_diagram,
                                      highlightbackground=cst.HCOLOR, state=tk.DISABLED)
        self.mcplotbutton.grid(row=9, columnspan=4)

        min_frame.place(relx=.5, rely=0, anchor="n")

        # PLOT CONTROLS #
        tk.Label(plt_frame, text='PLOT CONTROLS', font=('', cst.TITLESIZE, 'underline')).grid(columnspan=6)
        # vars
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
        # UI elements
        self.phase_label = tk.Label(plt_frame, text='phase =', state=tk.DISABLED)
        self.phase_label.grid(row=1, column=1, sticky=tk.E)
        self.phase_slider = tk.Scale(plt_frame, variable=self.phase, from_=0, to=1, orient=tk.HORIZONTAL,
                                     resolution=0.005, length=300, state=tk.DISABLED)
        self.phase_slider.grid(row=1, column=2, columnspan=4)
        self.phase_button = tk.Checkbutton(plt_frame, var=self.do_phasedot, command=self.toggle_dot, state=tk.DISABLED)
        self.phase_button.grid(row=1)

        self.plot_rv1data_label = tk.Label(plt_frame, text='Primary RV data', state=tk.DISABLED)
        self.plot_rv1data_label.grid(row=2, column=1)
        self.plot_rv1data_button = tk.Checkbutton(plt_frame, var=self.do_datarv1, state=tk.DISABLED)
        self.plot_rv1data_button.grid(row=2)

        self.plot_rv2data_label = tk.Label(plt_frame, text='Secondary RV data', state=tk.DISABLED)
        self.plot_rv2data_label.grid(row=3, column=1)
        self.plot_rv2data_button = tk.Checkbutton(plt_frame, var=self.do_datarv2, state=tk.DISABLED)
        self.plot_rv2data_button.grid(row=3)

        self.plot_asdata_label = tk.Label(plt_frame, text='Astrometric data', state=tk.DISABLED)
        self.plot_asdata_label.grid(row=4, column=1)
        self.plot_asdata_button = tk.Checkbutton(plt_frame, var=self.do_dataas, state=tk.DISABLED)
        self.plot_asdata_button.grid(row=4)

        self.plot_rv1model_label = tk.Label(plt_frame, text='Primary RV model', state=tk.DISABLED)
        self.plot_rv1model_label.grid(row=2, column=3)
        self.plot_rv1model_button = tk.Checkbutton(plt_frame, var=self.do_modelrv1, state=tk.DISABLED)
        self.plot_rv1model_button.grid(row=2, column=2)

        self.plot_rv2model_label = tk.Label(plt_frame, text='Secondary RV model', state=tk.DISABLED)
        self.plot_rv2model_label.grid(row=3, column=3)
        self.plot_rv2model_button = tk.Checkbutton(plt_frame, var=self.do_modelrv2, state=tk.DISABLED)
        self.plot_rv2model_button.grid(row=3, column=2)

        self.plot_asmodel_label = tk.Label(plt_frame, text='Model Orbit', state=tk.DISABLED)
        self.plot_asmodel_label.grid(row=4, column=3)
        self.plot_asmodel_button = tk.Checkbutton(plt_frame, var=self.do_modelas, state=tk.DISABLED)
        self.plot_asmodel_button.grid(row=4, column=2)

        self.plot_nodeline_label = tk.Label(plt_frame, text='Line of nodes', state=tk.DISABLED)
        self.plot_nodeline_label.grid(row=2, column=5)
        self.plot_nodeline_button = tk.Checkbutton(plt_frame, var=self.do_nodeline, state=tk.DISABLED)
        self.plot_nodeline_button.grid(row=2, column=4)

        self.plot_semimajor_label = tk.Label(plt_frame, text='Semi-major axis', state=tk.DISABLED)
        self.plot_semimajor_label.grid(row=3, column=5)
        self.plot_semimajor_button = tk.Checkbutton(plt_frame, var=self.do_semimajor, state=tk.DISABLED)
        self.plot_semimajor_button.grid(row=3, column=4)

        self.plot_peri_label = tk.Label(plt_frame, text='Periastron', state=tk.DISABLED)
        self.plot_peri_label.grid(row=4, column=5)
        self.plot_peri_button = tk.Checkbutton(plt_frame, var=self.do_peri, state=tk.DISABLED)
        self.plot_peri_button.grid(row=4, column=4)

        self.as_dist_label = tk.Label(plt_frame, text='Astrometric errors', state=tk.DISABLED)
        self.as_dist_label.grid(row=5, column=5)
        self.as_dist_button = tk.Checkbutton(plt_frame, var=self.do_as_dist, state=tk.DISABLED)
        self.as_dist_button.grid(row=5, column=4)

        self.do_legend = tk.BooleanVar()
        legend_button = tk.Checkbutton(plt_frame, var=self.do_legend, highlightbackground=cst.HCOLOR)
        legend_button.grid(row=5)
        tk.Label(plt_frame, text='Legend').grid(row=5, column=1)

        self.plot_vs_phase = tk.BooleanVar(value=False)
        self.pphase_but = tk.Radiobutton(plt_frame, text='phase', command=self.toggle_phase_time,
                                         variable=self.plot_vs_phase, value=True, state=tk.DISABLED)
        self.ptime_but = tk.Radiobutton(plt_frame, text='time', command=self.toggle_phase_time,
                                        variable=self.plot_vs_phase, value=False, state=tk.DISABLED)
        self.pphase_but.grid(row=5, column=2)
        self.ptime_but.grid(row=5, column=3)
        self.modelwidgets = {self.plot_asmodel_label, self.plot_asmodel_button,
                             self.plot_rv1model_button, self.plot_rv2model_button, self.plot_rv1model_label,
                             self.plot_rv2model_label, self.plot_semimajor_button, self.plot_semimajor_label,
                             self.plot_nodeline_button, self.plot_nodeline_label, self.plot_peri_label,
                             self.plot_peri_button, self.as_dist_button, self.as_dist_label,
                             self.pphase_but, self.ptime_but}
        refreshframe2 = tk.Frame(plt_frame_top)
        tk.Button(refreshframe2, text='Refresh Plots', width=20, height=2, command=self.update,
                  highlightbackground=cst.HCOLOR).pack()

        plt_frame.pack()
        refreshframe2.pack(pady=10)
        plt_frame_top.place(relx=.5, rely=0, anchor="n")

        # OUTPUT TAB #
        tk.Label(out_frame, text='Output file names').grid(columnspan=2)
        tk.Label(out_frame, text='Guessed parameters').grid(row=1)
        tk.Label(out_frame, text='Fitted parameters').grid(row=2)
        tk.Label(out_frame, text='Corner plot').grid(row=3)

        self.guess_out = tk.StringVar()
        self.fit_out = tk.StringVar()
        self.corner_out = tk.StringVar()

        tk.Entry(out_frame, textvariable=self.guess_out).grid(row=1, column=1)
        tk.Entry(out_frame, textvariable=self.fit_out).grid(row=2, column=1)
        tk.Entry(out_frame, textvariable=self.corner_out).grid(row=3, column=1)

        out_frame.pack()

        # setup the plotting windows
        self.init_plots()

        # pack everything neatly
        tabs.add(data_frame_tab, text='Data Files')
        tabs.add(guess_infer_tab, text='System/Parameters')
        tabs.add(min_frame_tab, text='Minimization')
        tabs.add(plt_frame_tab, text='Plot Controls')
        tabs.add(out_frame_tab, text='Output names')
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
        self.rv_fig = plt.figure(figsize=(10.5, 4.2))
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

    def toggle_dot(self):
        for widg in self.phase_slider, self.phase_label:
            self.toggle(widg, self.do_phasedot.get())

    def toggle_q(self):
        self.lock_q.set(not self.lock_q.get())
        if self.lock_q.get():
            self.param_var_list[8].set('q = k1/k2 =')
            self.toggle(self.vary_button_list[8], False)
            self.lock_q_button.config(text='k')
            if float(self.guess_var_list[8].get()) != 0:
                self.guess_var_list[8].set(
                    np.round(float(self.guess_var_list[7].get()) / float(self.guess_var_list[8].get()), 6))

        else:
            self.param_var_list[8].set('k2 (km/s) =')
            self.toggle(self.vary_button_list[8], True)
            self.lock_q_button.config(text='q')
            if float(self.guess_var_list[8].get()) != 0:
                self.guess_var_list[8].set(
                    np.round(float(self.guess_var_list[7].get()) / float(self.guess_var_list[8].get()), 6))

    def toggle_lock(self):
        self.lock_gs.set(not self.lock_gs.get())
        if self.lock_gs.get():
            self.lock_gs_button.config(image=self.unlocked_image)
            self.guess_var_list[10].set('')
            for widgset in self.param_label_list, self.vary_button_list, self.guess_entry_list:
                self.toggle(widgset[10], False)

        else:
            self.lock_gs_button.config(image=self.locked_image)
            for widgset in self.param_label_list, self.vary_button_list, self.guess_entry_list:
                self.toggle(widgset[10], True)
            self.guess_var_list[10].set(str(self.guess_var_list[9].get()))

    def toggle_phase_time(self):
        if not self.plot_vs_phase.get() and self.do_phasedot.get():
            self.do_phasedot.set(False)
        self.toggle(self.phase_button, self.plot_vs_phase.get())

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
        self.set_RV_or_AS_mode()

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
        self.set_RV_or_AS_mode()

    def toggle_as(self):
        """
        toggles the AS widgets
        """
        for widg in (self.as_file, self.as_label, self.seppa_but, self.en_but, self.plot_asdata_label,
                     self.plot_asdata_button):
            self.toggle(widg, self.include_as.get())
        if self.do_dataas.get():
            self.do_dataas.set(False)
        self.set_RV_or_AS_mode()

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
        for lst in self.param_label_list, self.vary_button_list:
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
            self.guess_dict = spl.guess_loader(self.wd.get(), self.guess_file.get())
        except IOError:
            print('cannot find your guess file!')
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
            self.guess_dict = None
            return
        self.set_system()

    def set_guess_dict_from_entries(self):
        """
        builds the guess dict from the guess column
        """
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
                           'gamma2': (float(self.guess_var_list[10].get()), self.vary_var_list[10].get())
                           if not self.lock_gs.get() else (
                               float(self.guess_var_list[9].get()), self.vary_var_list[9].get()),
                           'mt': (float(self.guess_var_list[11].get()), self.vary_var_list[11].get())
                           }

        self.param_dict = dict()
        for param, value in self.guess_dict.items():
            self.param_dict[param] = value[0]

    def set_system(self):
        """
        sets the system from the current guess column
        """
        try:
            self.set_guess_dict_from_entries()
            self.system = bsys.System(self.param_dict)
        except ValueError:
            print('invalid model!')
            self.guess_dict = None
            self.system = None
            self.plot_vs_phase.set(False)
            self.toggle_phase_time()
            for widg in self.modelwidgets:
                self.toggle(widg, False)
            return False
        else:
            for widg in self.modelwidgets:
                self.toggle(widg, True)
            return True

    def load_data(self):
        """
        loads the data from the currently selected files
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
            self.data_dict = spl.data_loader(self.wd.get(), filetypes, filenames, self.seppa.get())
        except (OSError, KeyError):
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
                    spm.LMminimizer(self.guess_dict, self.data_dict, self.do_mcmc.get(), self.steps.get(), w,
                                    self.lock_gs.get(), self.lock_q.get())
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
                    self.error_var_list[0].set(np.round(0 if pars['p'].stderr is None else pars['p'].stderr, 3))
                if pars['e'].vary:
                    self.mininimzed_var_list[1].set(np.round(pars['e'].value, 3))
                    self.error_var_list[1].set(np.round(0 if pars['e'].stderr is None else pars['e'].stderr, 3))
                if pars['i'].vary:
                    self.mininimzed_var_list[2].set(np.round(pars['i'].value % 360, 3))
                    self.error_var_list[2].set(np.round(0 if pars['i'].stderr is None else pars['i'].stderr, 3))
                if pars['omega'].vary:
                    self.mininimzed_var_list[3].set(np.round(pars['omega'].value % 360, 3))
                    self.error_var_list[3].set(np.round(0 if pars['omega'].stderr is None else pars['omega'].stderr, 3))
                if pars['Omega'].vary:
                    self.mininimzed_var_list[4].set(np.round(pars['Omega'].value % 360, 3))
                    self.error_var_list[4].set(np.round(0 if pars['Omega'].stderr is None else pars['Omega'].stderr, 3))
                if pars['t0'].vary:
                    self.mininimzed_var_list[5].set(np.round(pars['t0'].value, 3))
                    self.error_var_list[5].set(np.round(0 if pars['t0'].stderr is None else pars['t0'].stderr, 3))
                if pars['d'].vary:
                    self.mininimzed_var_list[6].set(np.round(pars['d'].value, 3))
                    self.error_var_list[6].set(np.round(0 if pars['d'].stderr is None else pars['d'].stderr, 3))
                if pars['k1'].vary:
                    self.mininimzed_var_list[7].set(np.round(pars['k1'].value, 3))
                    self.error_var_list[7].set(np.round(0 if pars['k1'].stderr is None else pars['k1'].stderr, 3))
                if pars['k2'].vary and not self.lock_q.get():
                    self.mininimzed_var_list[8].set(np.round(pars['k2'].value, 3))
                    self.error_var_list[8].set(np.round(0 if pars['k2'].stderr is None else pars['k2'].stderr, 3))
                if pars['gamma1'].vary:
                    self.mininimzed_var_list[9].set(np.round(pars['gamma1'].value, 3))
                    self.error_var_list[9].set(
                        np.round(0 if pars['gamma1'].stderr is None else pars['gamma1'].stderr, 3))
                if pars['gamma2'].vary:
                    self.mininimzed_var_list[10].set(np.round(pars['gamma2'].value, 3))
                    self.error_var_list[10].set(
                        np.round(0 if pars['gamma2'].stderr is None else pars['gamma2'].stderr, 3))
                if pars['mt'].vary:
                    self.mininimzed_var_list[11].set(np.round(pars['mt'].value, 3))
                    self.error_var_list[11].set(np.round(0 if pars['mt'].stderr is None else pars['mt'].stderr, 3))

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
        self.semimajork1k2.set(np.round(self.system.semimajor_axis_from_RV(), 2))
        self.semimajord.set(np.round(self.system.semimajor_axis_from_distance(), 2))

    def update(self):
        """
        updates the gui, by replotting everything that is selected
        """
        self.load_data()
        if self.plot_vs_phase.get() and not self.set_system():
            return
        if self.system is not None:
            self.set_inferred_params()
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

        if self.do_phasedot.get() and self.plot_vs_phase.get():
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
        if self.plot_vs_phase.get():
            phases, rv, err = self.system.create_phase_extended_RV(self.data_dict['RV1'], 0.15)
            self.rv1data_line = self.rv_ax.errorbar(phases, rv, yerr=err, ls='', capsize=0.1, marker='o',
                                                    ms=5, color='b')
            self.rv_ax.set_xlabel(PHASE_STR)
        else:
            self.rv1data_line = self.rv_ax.errorbar(self.data_dict['RV1']['hjds'], self.data_dict['RV1']['RVs'],
                                                    yerr=self.data_dict['RV1']['errors'], ls='', capsize=0.1,
                                                    marker='o', ms=5, color='b')
            self.rv_ax.set_xlabel(TIME_STR)

    def plot_rv2_data(self):
        """
        plot the rv2 data
        """
        if 'RV2' not in self.data_dict:
            return
        if self.rv2data_line is not None:
            self.rv2data_line.remove()
            self.rv2data_line = None
        if self.plot_vs_phase.get():
            phases, rv, err = self.system.create_phase_extended_RV(self.data_dict['RV2'], 0.15)
            self.rv2data_line = self.rv_ax.errorbar(phases, rv, yerr=err, ls='', capsize=0.1, marker='o',
                                                    ms=5, color='r')
            self.rv_ax.set_xlabel(PHASE_STR)
        else:
            self.rv2data_line = self.rv_ax.errorbar(self.data_dict['RV2']['hjds'], self.data_dict['RV2']['RVs'],
                                                    yerr=self.data_dict['RV2']['errors'], ls='', capsize=0.1,
                                                    marker='o', ms=5, color='r')
            self.rv_ax.set_xlabel(TIME_STR)

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
                (data['norths'][i], self.system.relative.north_of_hjd(data['hjds'][i])),
                c=(0.75, 0.25, 0.0, 0.9))[0])

    def plot_rv1_curve(self):
        """
        the the rv1 model curve
        """
        if self.plot_vs_phase.get():
            phases = np.linspace(-0.15, 1.15, num=150)
            vrads1 = self.system.primary.radial_velocity_of_phases(phases)
            if self.rv1_line is None:
                self.rv1_line, = self.rv_ax.plot(phases, vrads1, label=r'primary', color='b', ls='--')
            else:
                self.rv1_line.set_xdata(phases)
                self.rv1_line.set_ydata(vrads1)
            self.rv_ax.set_xlabel(PHASE_STR)
        else:
            m = np.infty
            mm = -np.infty
            if self.include_rv1.get():
                m = min(m, min(self.data_dict['RV1']['hjds']))
                mm = max(mm, max(self.data_dict['RV1']['hjds']))
            if self.include_rv2.get():
                m = min(m, min(self.data_dict['RV2']['hjds']))
                mm = max(mm, max(self.data_dict['RV2']['hjds']))
            times = np.linspace(m - 0.01 * (mm - m), m - 0.01 * (mm - m) + self.system.p, endpoint=False, num=100)
            rvs = self.system.primary.radial_velocity_of_phases(self.system.phase_of_hjds(times))
            times, rvs = self.system.extend_rvs_until_time(times, rvs, mm)
            if self.rv1_line is None:
                self.rv1_line, = self.rv_ax.plot(times, rvs, label=r'primary', color='b', ls='--')
            else:
                self.rv1_line.set_xdata(times)
                self.rv1_line.set_ydata(rvs)
            self.rv_ax.set_xlabel(TIME_STR)

    def plot_rv2_curve(self):
        """
        plot the rv2 model curve
        """
        if self.plot_vs_phase.get():
            phases = np.linspace(-0.15, 1.15, num=150)
            vrads1 = self.system.secondary.radial_velocity_of_phases(phases)
            if self.rv2_line is None:
                self.rv2_line, = self.rv_ax.plot(phases, vrads1, label=r'secondary', color='r', ls='--')
            else:
                self.rv2_line.set_xdata(phases)
                self.rv2_line.set_ydata(vrads1)
            self.rv_ax.set_xlabel(PHASE_STR)
        else:
            m = np.infty
            mm = -np.infty
            if self.include_rv2.get():
                m = min(m, min(self.data_dict['RV2']['hjds']))
                mm = max(mm, max(self.data_dict['RV2']['hjds']))
            if self.include_rv1.get():
                m = min(m, min(self.data_dict['RV1']['hjds']))
                mm = max(mm, max(self.data_dict['RV1']['hjds']))
            times = np.linspace(m, m + self.system.p, num=100)
            rvs = self.system.secondary.radial_velocity_of_phases(self.system.phase_of_hjds(times))
            times, rvs = self.system.extend_rvs_until_time(times, rvs, mm)
            if self.rv2_line is None:
                self.rv2_line, = self.rv_ax.plot(times, rvs, label=r'secondary', color='r', ls='--')
            else:
                self.rv2_line.set_xdata(times)
                self.rv2_line.set_ydata(rvs)
            self.rv_ax.set_xlabel(TIME_STR)

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

    def make_corner_diagram(self):
        """
        plot a corner diagram of an MCMC run
        """
        if self.didmcmc:
            corner = spp.plot_corner_diagram(self.minresult)
            out = self.corner_out.get()
            if out == '':
                out = 'corner_diagram'
            corner.savefig(self.wd.get() + '/' + out + '{}.png'.format(self.mcmc_run_number))
            plt.close(corner)
        else:
            print('do an mcmc minimization first!')

    def save_params(self):
        """
        save minimized parameters to a file
        """
        out = self.fit_out.get()
        if out == '':
            out = 'fitted_params'
        with open(self.wd.get() + '/' + out + '{}.txt'.format(self.minimization_run_number), 'w') as f:
            for i in range(len(self.mininimzed_var_list)):
                f.write(str(self.param_var_list[i].get()) + ' ' + str(self.mininimzed_var_list[i].get()) + ' ' + str(
                    self.error_var_list[i].get()) + '\n')
            f.write('reduced chisq = {} \n'.format(self.redchisq.get()))
            f.write('dof = {} \n'.format(self.dof.get()))
            f.write('param order: {}'.format(self.minresult.var_names))
        np.savetxt(self.wd.get() + '/' + out + '_covar{}.txt'.format(self.minimization_run_number),
                   self.minresult.covar)

    def save_guesses(self):
        """
        save guesses parameters to a file
        """
        out = self.guess_out.get()
        if out == '':
            out = 'guessed_params'
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
        root.geometry("{}x{}+0+0".format(int(0.35 * w), h))
        root.title('spinOSgui')
        SpinOSGUI(root, wd, w, h)

    root.mainloop()
