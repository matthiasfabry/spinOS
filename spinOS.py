"""
This is spinOS, the SPectroscopic and INterferometric Orbiral Solution finder, put in a graphical user interface

Goal:
This is a graphical user interface implementation of the command line version spinOS. It uses (spectroscopic or
otherwise) Radial Velocity (RV) measurements of either of the two components of the binary, as well as (
interferometric or otherwise) relative astrometry (AS) measurements of the binary. You then need to supply a guess
for the parameters that define the binary orbit:
    -) e:       the eccentricity of the orbit
    -) i:       the inclination of the orbit (with respect to the plane of the sky)
    -) omega:   the argument of periastron of the secondary, with respect to its ascending node
    -) Omega:   the longitude of the ascending node of the secondary, measured East from North
    -) t0:      the time of periastron passage (this number should be between 0 and the supposed period)
    -) p:       the period of the binary
    -) d:       the distance to the system
and:
    -) k1:      the semiamplitude of the RV curve of the primary
    -) k2:      the semiamplitude of the RV curve of the secondary
    -) gamma1:  the peculiar velocity of the primary
    -) gamma2:  the peculiar velocity of the secondary
if you have an SB2 with or without astrometric data, or:
    -) k1:      the semiamplitude of the RV curve of the primary
    -) gamma1:  the peculiar velocity of the primary
    -) mt:      the total dynamical mass of the system (which sets the apparent size of the orbit)
if you an SB1 with or without astrometric data, or:
    -) mt:      the total dynamical mass of the system (which sets the apparent size of the orbit)
if you have an astrometric orbit only.

This application allows for easy plotting of data and models, as well as minimization of the model to your supplied
data. The program then gives a best fit value for the parameters itemized above, as well as the component masses.
Error are estimated as the diagonal elements of the correlation matrix, or as half of the difference
between the 15.87 and 84.13 percentiles found in an Markov Chain Monte Carlo sampling.

Usage:
To start spinOS, simply run:
 python3 spinOS.py

At the top, put your working directory, where all your data files are (typically <objectname/> don't forget the slash!)

The application expects the data to be in the following format: All data files should be plain text files, formatted as:
for RV data:
 JD(days) RV(km/s) error_on_RV(km/s)
eg:
 45000 25.1 2.1
 45860 -4.2 1.1
 etc...

for AS data:
either:
 JD(days) E_separation(mas) N_separation(mas) semimajor_ax_errorellipse(mas)
                                                            semiminor_ax_errorellipse(mas) angle_E_of_N_of_major_ax(deg)
eg:
 48000 -2.5 2.4 0.1 0.8 60
 48050 2.1 8.4 0.4 0.5 90
 etc...

or:
 JD(days) separation(mas) PA(deg) semimajor_ax_errorellipse(mas)
                                                            semiminor_ax_errorellipse(mas) angle_E_of_N_of_major_ax(deg)
eg:
 48000 3.5 316 0.1 0.8 60
 48050 8.7 76 0.4 0.5 90
 etc...

for the guess file, format should be eg:
 e 0.648 True
 i 86.53 True
 omega 211.0 True
 Omega 67.3 True
 t0 56547.1 True
 k1 31.0 False
 k2 52.0 True
 p 3252.0 True
 gamma1 15.8 False
 gamma2 5.6 False
 mt 30.0 True
All eleven parameters should be guessed.

Use the provided buttons to load data and guesses from the files designated. The Plot control buttons allow plotting of
the relevant data and models.
You can minimize the model to the selected data with the minimize button, with or without an mcmc error estimation.
The save buttons save either the guesses to guesses.txt or minimized parameters to params_runi.txt (i is number of
minimization runs performed in this session).
If the last minimization run contained an MCMC analysis, you can create a corner plot with the button provided. It will
be saved at corneri.png (i is run number).


Dependencies:
    python 3.7.7
    numpy 1.18.1
    scipy 1.4.1
    lmfit 0.9.14
    matplotlib 3.1.3
    emcee 3.0.0 (if MCMC error calculation is performed)
    corner 2.0.1 (if MCMC error calculation is performed)

Author:
    Matthias Fabry
    Instituut voor Sterrekunde, KU Leuven, Belgium

Date:
    1 Apr 2020

Version:
    2.2

Acknowledgements:
    This python3 implementation is heavily based on an earlier IDL implementation by Hugues Sana.
    We thank the authors of lmfit for the development of their package.
"""

import sys
import tkinter as tk

import lmfit as lm
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import EllipseCollection

from modules import spinOSio as spl, spinOSminimizer as spm, spinOSplotter as spp, binary_system as bsys


class SpinOSApp:
    def __init__(self, master, wwd):
        # set the root frame
        self.frame = tk.Frame(master)
        # set the data frame
        data_frame = tk.Frame(self.frame)
        data_frame.grid(row=0, columnspan=2)
        # set the guess frame
        guess_frame = tk.Frame(self.frame)
        guess_frame.grid(row=1, columnspan=2)
        # set the inferation frame
        infer_frame = tk.Frame(self.frame)
        infer_frame.grid(row=2)
        # set the minimization frame
        min_frame = tk.Frame(self.frame)
        min_frame.grid(row=2, column=1)
        # set the p[ot window controls frame
        plt_frame = tk.Frame(self.frame)
        plt_frame.grid(row=3, columnspan=2)

        # initialize some variables
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
        tk.Label(data_frame, text='DATA', font=('', titlesize)).grid(columnspan=5, sticky=tk.N)

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

        self.seppa = tk.BooleanVar()
        self.seppa.set(True)
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
        tk.Label(guess_frame, text='MODEL/GUESS PARAMETERS', font=('', titlesize)).grid(row=0, columnspan=columns)
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
        self.guess_var_list = [tk.StringVar() for _ in range(numofparams)]
        for i in range(numofparams):
            self.guess_var_list[i].set('1.0')

        # define entry boxes
        self.guess_entry_list = [tk.Entry(guess_frame, textvariable=self.guess_var_list[i], width=10) for i in
                                 range(numofparams)]
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
        tk.Button(guess_frame, text='Load guesses', command=self.set_guesses_from_file,
                  highlightbackground=hcolor).grid(row=numofparams + 2)
        tk.Button(guess_frame, text='Save guesses', command=self.save_guesses,
                  highlightbackground=hcolor).grid(row=numofparams + 2, column=1)
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
        tk.Label(infer_frame, text='INFERRED PARAMETERS', font=('', titlesize)).grid(columnspan=4, sticky=tk.N)
        tk.Label(infer_frame, text='From k_1/k_2').grid(row=1, columnspan=2)
        tk.Label(infer_frame, text='M1 (M_sun) =').grid(row=3, sticky=tk.E)
        tk.Label(infer_frame, text='M2 (M_sun) =').grid(row=4, sticky=tk.E)
        tk.Label(infer_frame, text='a (AU) =').grid(row=2, sticky=tk.E)
        tk.Label(infer_frame, textvariable=self.mprimary).grid(row=3, column=1)
        tk.Label(infer_frame, textvariable=self.msecondary).grid(row=4, column=1)
        tk.Label(infer_frame, textvariable=self.semimajork1k2).grid(row=2, column=1)
        tk.Label(infer_frame, text='From d/M_tot:').grid(row=1, column=3, columnspan=2)
        tk.Label(infer_frame, text='a (AU) =').grid(row=2, column=3, sticky=tk.E)
        tk.Label(infer_frame, textvariable=self.semimajord).grid(row=2, column=4)

        # MINIMIZATION FRAME #
        # define variables
        self.do_mcmc = tk.BooleanVar()
        self.redchisq = tk.DoubleVar()
        self.dof = tk.IntVar()
        self.steps = tk.IntVar(value=1000)

        # define labels and buttons in a grid
        tk.Label(min_frame, text='MINIMIZATION', font=('', titlesize)).grid(row=0, columnspan=4)

        tk.Label(min_frame, text='Do MCMC?:').grid(row=1, sticky=tk.E)
        mcmc_check = tk.Checkbutton(min_frame, var=self.do_mcmc, command=self.toggle_mc)
        mcmc_check.grid(row=1, column=1)
        self.steps_label = tk.Label(min_frame, text='# of samples:', state=tk.DISABLED)
        self.steps_label.grid(row=1, column=2, sticky=tk.E)
        self.steps_entry = tk.Entry(min_frame, textvariable=self.steps, width=5, state=tk.DISABLED)
        self.steps_entry.grid(row=1, column=3)
        tk.Label(min_frame, text='Red. Chi Sqrd =').grid(row=2, columnspan=2, sticky=tk.E)
        tk.Label(min_frame, text='Deg. of frdm =').grid(row=3, columnspan=2, sticky=tk.E)
        tk.Label(min_frame, textvariable=self.redchisq).grid(row=2, column=2, columnspan=2, sticky=tk.W)
        tk.Label(min_frame, textvariable=self.dof).grid(row=3, column=2, columnspan=2, sticky=tk.W)
        tk.Button(min_frame, text='Minimize!', command=self.minimize, highlightbackground=hcolor).grid(
            row=4, columnspan=2)
        tk.Button(min_frame, text='Make MCMC plot', command=self.plot_corner_diagram,
                  highlightbackground=hcolor).grid(row=4, column=2, columnspan=2)

        # PLOT WINDOW CONTROLS
        tk.Label(plt_frame, text='PLOT CONTROLS', font=('', titlesize)).grid(columnspan=5)

        self.do_phasedot = tk.BooleanVar(False)
        self.do_datarv1 = tk.BooleanVar(False)
        self.do_datarv2 = tk.BooleanVar(False)
        self.do_dataas = tk.BooleanVar(False)
        self.do_modelrv1 = tk.BooleanVar(False)
        self.do_modelrv2 = tk.BooleanVar(False)
        self.do_modelas = tk.BooleanVar(False)
        self.do_nodeline = tk.BooleanVar(False)
        self.do_semimajor = tk.BooleanVar(False)
        self.do_peri = tk.BooleanVar(False)

        self.common_plot_bools = [self.do_phasedot]
        self.common_plot_funcs = [self.plot_dots]
        self.common_plot_line_groups = [[self.as_dot, self.rv1_dot, self.rv2_dot]]
        self.rv_plot_bools = [self.do_datarv1, self.do_datarv2, self.do_modelrv1, self.do_modelrv2]
        self.rv_plot_funcs = [self.plot_rv1_data, self.plot_rv2_data, self.plot_rv1_curve, self.plot_rv2_curve]
        self.rv_plot_line_groups = [[self.rv1data_line], [self.rv2data_line], [self.rv1_line], [self.rv2_line]]
        self.as_plot_bools = [self.do_dataas, self.do_modelas, self.do_nodeline, self.do_semimajor, self.do_peri]
        self.as_plot_funcs = [self.plot_as_data, self.plot_relative_orbit, self.plot_node_line,
                              self.plot_semimajor_axis, self.plot_periastron]
        self.as_plot_line_groups = [[self.asdata_line, self.as_ellipses], [self.as_line], [self.node_line],
                                    [self.semi_major], [self.peri_dot]]

        self.phase: tk.DoubleVar = tk.DoubleVar()
        self.phase.trace_add('write', lambda n, ix, m: self.update())

        self.phase_label = tk.Label(plt_frame, text='phase =', state=tk.DISABLED)
        self.phase_label.grid(row=1, column=1, sticky=tk.E)
        self.phase_slider = tk.Scale(plt_frame, variable=self.phase, from_=0, to=1, orient=tk.HORIZONTAL,
                                     resolution=0.005, length=300, state=tk.DISABLED)
        self.phase_slider.grid(row=1, column=2, columnspan=4)
        self.phase_button = tk.Checkbutton(plt_frame, var=self.do_phasedot, command=self.toggle_dots,
                                           state=tk.DISABLED)
        self.phase_button.grid(row=1)

        self.plot_rv1data_label = tk.Label(plt_frame, text='Primary RV data', state=tk.DISABLED)
        self.plot_rv1data_label.grid(row=2, column=1)
        self.plot_rv1data_button = tk.Checkbutton(plt_frame, var=self.do_datarv1, command=self.update,
                                                  state=tk.DISABLED)
        self.plot_rv1data_button.grid(row=2)

        self.plot_rv2data_label = tk.Label(plt_frame, text='Secondary RV data', state=tk.DISABLED)
        self.plot_rv2data_label.grid(row=3, column=1)
        self.plot_rv2data_button = tk.Checkbutton(plt_frame, var=self.do_datarv2, command=self.update,
                                                  state=tk.DISABLED)
        self.plot_rv2data_button.grid(row=3)

        self.plot_asdata_label = tk.Label(plt_frame, text='Astrometric data', state=tk.DISABLED)
        self.plot_asdata_label.grid(row=4, column=1)
        self.plot_asdata_button = tk.Checkbutton(plt_frame, var=self.do_dataas, command=self.update,
                                                 state=tk.DISABLED)
        self.plot_asdata_button.grid(row=4)

        self.plot_rv1model_label = tk.Label(plt_frame, text='Primary RV model', state=tk.DISABLED)
        self.plot_rv1model_label.grid(row=2, column=3)
        self.plot_rv1model_button = tk.Checkbutton(plt_frame, var=self.do_modelrv1, command=self.update,
                                                   state=tk.DISABLED)
        self.plot_rv1model_button.grid(row=2, column=2)

        self.plot_rv2model_label = tk.Label(plt_frame, text='Secondary RV model', state=tk.DISABLED)
        self.plot_rv2model_label.grid(row=3, column=3)
        self.plot_rv2model_button = tk.Checkbutton(plt_frame, var=self.do_modelrv2, command=self.update,
                                                   state=tk.DISABLED)
        self.plot_rv2model_button.grid(row=3, column=2)

        self.plot_asmodel_label = tk.Label(plt_frame, text='Model Orbit', state=tk.DISABLED)
        self.plot_asmodel_label.grid(row=4, column=3)
        self.plot_asmodel_button = tk.Checkbutton(plt_frame, var=self.do_modelas, command=self.update,
                                                  state=tk.DISABLED)
        self.plot_asmodel_button.grid(row=4, column=2)

        self.plot_nodeline_label = tk.Label(plt_frame, text='Line of nodes', state=tk.DISABLED)
        self.plot_nodeline_label.grid(row=2, column=5)
        self.plot_nodeline_button = tk.Checkbutton(plt_frame, var=self.do_nodeline, command=self.update,
                                                   state=tk.DISABLED)
        self.plot_nodeline_button.grid(row=2, column=4)

        self.plot_semimajor_label = tk.Label(plt_frame, text='Semi-major axis', state=tk.DISABLED)
        self.plot_semimajor_label.grid(row=3, column=5)
        self.plot_semimajor_button = tk.Checkbutton(plt_frame, var=self.do_semimajor, command=self.update,
                                                    state=tk.DISABLED)
        self.plot_semimajor_button.grid(row=3, column=4)

        self.plot_peri_label = tk.Label(plt_frame, text='Periastron', state=tk.DISABLED)
        self.plot_peri_label.grid(row=4, column=5)
        self.plot_peri_button = tk.Checkbutton(plt_frame, var=self.do_peri, command=self.update,
                                               state=tk.DISABLED)
        self.plot_peri_button.grid(row=4, column=4)

        self.modelwidgets = {self.phase_label, self.phase_slider, self.plot_asmodel_label, self.plot_asmodel_button,
                             self.plot_rv1model_button, self.plot_rv2model_button, self.plot_rv1model_label,
                             self.plot_rv2model_label, self.plot_semimajor_button, self.plot_semimajor_label,
                             self.plot_nodeline_button, self.plot_nodeline_label, self.plot_peri_label,
                             self.plot_peri_button}
        self.legend_button_bool = True
        legend_button = tk.Button(plt_frame, text='Legend', highlightbackground=hcolor,
                                  command=lambda: self.toggle_legends(legend_button))
        legend_button.grid(row=5, column=3)

        # BACK TO ROOT FRAME #
        # display the root frame
        self.init_plots()
        self.frame.pack()

    @staticmethod
    def toggle(widg, boolvalue):
        if boolvalue:
            widg.config(state=tk.NORMAL)
        else:
            widg.config(state=tk.DISABLED)

    def toggle_dots(self):
        for widg in self.phase_slider, self.phase_label:
            self.toggle(widg, self.do_phasedot.get())
        self.update()

    def toggle_mc(self):
        for widg in self.steps_label, self.steps_entry:
            self.toggle(widg, self.do_mcmc.get())

    def toggle_rv1(self):
        for widg in self.rv1_file, self.rv1_label, self.plot_rv1data_label, self.plot_rv1data_button:
            self.toggle(widg, self.include_rv1.get())
        self.toggle_databutton()
        self.set_RV_or_AS_mode()

    def toggle_rv2(self):
        if not self.include_rv1.get():
            self.include_rv2.set(False)
            return
        for widg in self.rv2_file, self.rv2_label, self.plot_rv2data_label, self.plot_rv2data_button:
            self.toggle(widg, self.include_rv2.get())
        self.toggle_databutton()
        self.set_RV_or_AS_mode()

    def toggle_as(self):
        for widg in self.as_file, self.as_label, self.seppa_but, self.en_but, \
                    self.plot_asdata_label, self.plot_asdata_button:
            self.toggle(widg, self.include_as.get())
        self.toggle_databutton()
        self.set_RV_or_AS_mode()

    def toggle_databutton(self):
        if not (self.include_rv1.get() or self.include_rv2.get() or self.include_as.get()):
            self.toggle(self.data_button, False)
        else:
            self.toggle(self.data_button, True)

    def set_RV_or_AS_mode(self):
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
                for i in {2, 4, 11}:
                    lst[i].config(state=tk.NORMAL)
                for i in {6, 7, 8, 9, 10}:
                    lst[i].config(state=tk.DISABLED)
            else:
                for i in {2, 4, 6, 7, 8, 9, 10, 11}:
                    lst[i].config(state=tk.NORMAL)

    def toggle_legends(self, button):
        if self.legend_button_bool:
            self.legend_button_bool = False
            button.config(relief=tk.SUNKEN)
            self.plot_legends()
        else:
            try:
                self.as_ax.get_legend().remove()
            except AttributeError:
                pass
            try:
                self.rv_ax.get_legend().remove()
            except AttributeError:
                pass
            self.legend_button_bool = True
            self.rv_fig.canvas.draw_idle()
            self.as_fig.canvas.draw_idle()
            button.config(relief=tk.RAISED)

    def init_plots(self):
        if self.rv_fig is not None:
            plt.close(self.rv_fig)
        if self.as_fig is not None:
            plt.close(self.as_fig)
        self.rv_fig = plt.figure(figsize=(10.5, 4.5))
        self.as_fig = plt.figure(figsize=(10.5, 4.5))
        move_figure(self.rv_fig, int(w / 3) + 20, 0)
        move_figure(self.as_fig, int(w / 3) + 20, int(h / 2) + 20)
        self.rv_ax = self.rv_fig.add_subplot(111)
        self.as_ax = self.as_fig.add_subplot(111, aspect=1)
        spp.setup_rvax(self.rv_ax)
        spp.setup_asax(self.as_ax)
        self.rv_fig.tight_layout()
        self.as_fig.tight_layout()
        plt.ion()
        plt.show()

    def transfer(self, varno):
        self.guess_var_list[varno].set(self.mininimzed_var_list[varno].get())

    def set_guesses_from_file(self):
        try:
            self.loading_guesses = True
            self.guess_dict = spl.guess_loader(self.wd.get(), self.guess_file.get())
            self.fill_guess_entries_from_dict()
        except IOError:
            print('cannot find your guess file!')
            self.guess_dict = None
        except (ValueError, KeyError):
            self.guess_dict = None
            print('some parameter has not been set properly')
        finally:
            self.loading_guesses = False
            self.update()

    def set_guess_dict_from_entries(self):
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

    def fill_guess_entries_from_dict(self):
        # here we must convert from dict to list, no way to write this faster
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

    def set_system(self):
        try:
            self.set_guess_dict_from_entries()
            self.system = bsys.System(self.param_dict)
        except ValueError:
            print('invalid model!')
            self.guess_dict = None
            self.system = None
            for widg in self.modelwidgets:
                self.toggle(widg, False)
        else:
            for widg in self.modelwidgets:
                self.toggle(widg, True)
            self.mprimary.set(np.round(self.system.primary_mass(), 2))
            self.msecondary.set(np.round(self.system.secondary_mass(), 2))
            self.semimajork1k2.set(np.round(self.system.semimajor_axis_from_RV(), 2))
            self.semimajord.set(np.round(self.system.semimajor_axis_from_distance(), 2))

    def load_data(self):
        filetypes = list()
        filenames = list()
        if self.rv1_file.get() is not '' and self.include_rv1.get():
            filetypes.append('RV1file')
            filenames.append(self.rv1_file.get())
        if self.rv2_file.get() is not '' and self.include_rv2.get():
            filetypes.append('RV2file')
            filenames.append(self.rv2_file.get())
        if self.as_file.get() is not '' and self.include_as.get():
            filetypes.append('ASfile')
            filenames.append(self.as_file.get())
        if len(filenames) > 0:
            try:
                self.data_dict = dict()
                self.data_dict = spl.data_loader(self.wd.get(), filetypes, filenames, self.seppa.get())
            except OSError as e:
                print('Cannot find datafile: ', e)
                self.data_dict = None
        else:
            print('no data files entered, or no data included!')
            self.data_dict = None

    def minimize(self):
        self.set_guess_dict_from_entries()
        self.load_data()
        if self.guess_dict is not None and self.data_dict is not None:
            # calculate best parameters
            try:
                self.minresult = spm.LMminimizer(self.guess_dict, self.data_dict, self.do_mcmc.get(), self.steps.get())
                if self.do_mcmc.get():
                    self.didmcmc = True
                    self.mcmc_run_number += 1
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
                self.minimization_run_number += 1
            except ValueError as e:
                print(e)

    def update(self):
        if self.loading_guesses:
            return
        self.load_data()
        self.set_system()

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
        if self.system:
            if self.do_phasedot.get():
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
        self.plot_legends()
        self.relim_plots()
        self.rv_fig.canvas.draw()
        self.rv_fig.canvas.flush_events()
        self.as_fig.canvas.draw()
        self.as_fig.canvas.flush_events()

    def relim_plots(self):
        for plot_bool in self.rv_plot_bools:
            if plot_bool.get():
                self.rv_ax.relim()
                self.rv_ax.axis('auto')

        for plot_bool in self.as_plot_bools:
            if plot_bool.get():
                self.as_ax.relim()
                self.as_ax.axis('image')

    def plot_rv1_data(self):
        if 'RV1' not in self.data_dict:
            return
        phases, rv, err = self.system.create_phase_extended_RV(self.data_dict['RV1'], 0.15)
        if self.rv1data_line is None:
            self.rv1data_line = self.rv_ax.errorbar(phases, rv, yerr=err, ls='', capsize=0.1, marker='o',
                                                    ms=5, color='b', label='Primary RV')
        else:
            self.rv1data_line[0].set_xdata(phases)
            self.rv1data_line[0].set_ydata(rv)

    def plot_rv2_data(self):
        if 'RV2' not in self.data_dict:
            return
        phases, rv, err = self.system.create_phase_extended_RV(self.data_dict['RV2'], 0.15)
        if self.rv2data_line is None:
            self.rv2data_line = self.rv_ax.errorbar(phases, rv, yerr=err, ls='', capsize=0.1, marker='o',
                                                    ms=5, color='r', label='Secondary RV')
        else:
            self.rv2data_line[0].set_xdata(phases)
            self.rv2data_line[0].set_ydata(rv)

    def plot_as_data(self):
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

    def plot_rv1_curve(self):
        phases = np.linspace(-0.15, 1.15, num=200)
        vrads1 = self.system.primary.radial_velocity_of_phases(phases)
        if self.rv1_line is None:
            self.rv1_line, = self.rv_ax.plot(phases, vrads1, label='primary', color='b', ls='--')
        else:
            self.rv1_line.set_ydata(vrads1)

    def plot_rv2_curve(self):
        phases = np.linspace(-0.15, 1.15, num=200)
        vrads2 = self.system.secondary.radial_velocity_of_phases(phases)
        if self.rv2_line is None:
            self.rv2_line, = self.rv_ax.plot(phases, vrads2, label='secondary', color='r', ls='--')
        else:
            self.rv2_line.set_ydata(vrads2)

    def plot_relative_orbit(self):
        ecc_anoms = np.linspace(0, 2 * np.pi, 200)
        norths = self.system.relative.north_of_ecc(ecc_anoms)
        easts = self.system.relative.east_of_ecc(ecc_anoms)
        if self.as_line is None:
            self.as_line, = self.as_ax.plot(easts, norths, label='relative orbit', color='k')
        else:
            self.as_line.set_xdata(easts)
            self.as_line.set_ydata(norths)

    def plot_node_line(self):
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
        system = self.system.relative
        if self.peri_dot is None:
            self.peri_dot, = self.as_ax.plot([system.east_of_ecc(0)], [system.north_of_ecc(0)], color='b', marker='s',
                                             ls='', fillstyle='full', label='periastron', markersize=8)
        else:
            self.peri_dot.set_xdata(system.east_of_ecc(0))
            self.peri_dot.set_ydata(system.north_of_ecc(0))

    def plot_semimajor_axis(self):
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
        if self.do_modelrv1.get() or self.do_datarv1.get():
            rv1 = self.system.primary.radial_velocity_of_phase(self.phase.get())
            self.rv1_dot = self.rv_ax.scatter(self.phase.get(), rv1, s=100, color='b', marker='x',
                                              label=np.round(rv1, 2))
        if self.rv2_dot is not None:
            self.rv2_dot.remove()
        if self.do_modelrv2.get() or self.do_datarv2.get():
            rv2 = self.system.secondary.radial_velocity_of_phase(self.phase.get())
            self.rv2_dot = self.rv_ax.scatter(self.phase.get(), rv2, s=100, color='r', marker='x',
                                              label=np.round(rv2, 2))
        if self.as_dot is not None:
            self.as_dot.remove()
        if self.do_modelas.get() or self.do_dataas.get():
            N = self.system.relative.north_of_ph(self.phase.get())
            E = self.system.relative.east_of_ph(self.phase.get())
            self.as_dot = self.as_ax.scatter(E, N, s=100, color='r', marker='x',
                                             label='{}E/{}N'.format(np.round(E, 2), np.round(N, 2)))

    def plot_legends(self):
        if self.legend_button_bool:
            return
        try:
            self.rv_ax.get_legend().remove()
        except AttributeError:
            pass
        try:
            self.as_ax.get_legend().remove()
        except AttributeError:
            pass
        if len(self.rv_ax.get_lines()) > 1:
            self.rv_ax.legend()
        if len(self.as_ax.get_lines()) > 1:
            self.as_ax.legend()

    def plot_corner_diagram(self):
        if self.didmcmc:
            corner = spp.plot_corner_diagram(self.minresult)
            corner.savefig(self.wd.get() + 'corner{}.png'.format(self.mcmc_run_number))
            plt.close(corner)
        else:
            print('do an mcmc minimization first!')

    def save_params(self):
        with open(self.wd.get() + 'params_run{}.txt'.format(self.minimization_run_number), 'w') as f:
            for i in range(len(self.mininimzed_var_list)):
                f.write(str(self.param_name_vars[i].get()) + ' ' + str(self.mininimzed_var_list[i].get()) + ' ' + str(
                    self.error_var_list[i].get()) + '\n')
            f.write('reduced chisq = {} \n'.format(self.redchisq.get()))
            f.write('dof = {} \n'.format(self.dof.get()))
            f.write(lm.report_fit(self.minresult.params))

    def save_guesses(self):
        self.set_guess_dict_from_entries()
        spl.guess_saver(self.wd.get(), self.guess_dict)


def move_figure(f, x, y):
    f.canvas.manager.window.wm_geometry("+{}+{}".format(x, y))


mpl.use("TkAgg")  # set the backend

hcolor = '#3399ff'
try:
    wd = sys.argv[1]
except IndexError:
    wd = ''

root = tk.Tk()
w, h = root.winfo_screenwidth(), root.winfo_screenheight()
root.geometry("{}x{}".format(int(w / 3), h))
root.title('spinOSgui')
app = SpinOSApp(root, wd)
root.mainloop()
