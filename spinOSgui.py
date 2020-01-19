"""
This is spinOSgui, the SPectroscopic and INterferometric Orbiral Solution finder, put in a graphical user interface

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
    -) k1:      the semiamplitude of the RV curve of the primary
    -) k2:      the semiamplitude of the RV curve of the secondary
    -) p:       the period of the binary
    -) gamma1:  the peculiar velocity of the primary
    -) gamma2:  the peculiar velocity of the secondary
    -) d:       the distance to the system

This application allows for easy plotting of data and models, as well as minimization of the model to your supplied
data. The program then gives a best fit value for the parameters itemized above, as well as the component masses.
Errors are calculated using a Markov Chain Monte Carlo (MCMC) method, the reported errors are half of the difference
between the 15.87 and 84.13 percentiles found in the MCMC sampling.

Usage:
To start spinOSgui, simply run:
 python3 spinOSgui.py

The application expects the data to be in the following format: All data files should be plain text files, with:
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
 JD(days) separation(mas) PA semimajor_ax_errorellipse(mas)
                                                            semiminor_ax_errorellipse(mas) angle_E_of_N_of_major_ax(deg)
eg:
 48000 3.5 316 0.1 0.8 60
 48050 8.7 76 0.4 0.5 90
 etc...

Dependencies:
    python 3.7
    numpy 1.17.2
    scipy 1.3.1
    lmfit 0.9.14
    matplotlib 3.1.1
    emcee 3.0.0 (if MCMC error calculation is performed)

Author:
    Matthias Fabry
    Instituut voor Sterrekunde, KU Leuven, Belgium

Date:
    21 Nov 2019

Version:
    1.5

Acknowledgements:
    This python3 implementation is heavily based on an earlier IDL implementation by Hugues Sana.
    We thank the authors of lmfit for the development of their package.

"""

import tkinter as tk
import sys
import lmfit as lm
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

import binarySystem as bsys
import spinOSio as spl
import spinOSminimizer as spm
import spinOSplotter as spp
import plotsave as pls


class SpinOSApp:
    def __init__(self, master, wd):
        # set the root frame
        self.frame = tk.Frame(master)

        # set the guess frame
        guess_frame = tk.Frame(self.frame)
        guess_frame.grid(row=1, column=0, sticky=tk.W)
        # set the inferation frame
        infer_frame = tk.Frame(self.frame)
        infer_frame.grid(row=2, column=0, sticky=tk.W)
        # set the data frame
        data_frame = tk.Frame(self.frame)
        data_frame.grid(row=0, column=0, sticky=tk.W)
        # set the minimization frame
        min_frame = tk.Frame(self.frame)
        min_frame.grid(row=3, column=0, sticky=tk.W)
        # set the plotting frame
        plot_frame = tk.Frame(self.frame)
        plot_frame.grid(padx=200, row=0, rowspan=4, column=1, sticky=tk.E)

        master.grid_columnconfigure(0, weight=0)
        master.grid_columnconfigure(1, weight=1)

        # initialize some variables that will be set later
        self.minimization_run_number = 0
        self.params = None
        self.guess_dict = None
        self.data_dict = None
        self.minresult = None
        self.system = None
        self.primaryline = None
        self.secondaryline = None
        self.relativeline = None
        self.nodeline = None
        self.axisline = None
        self.periastron = None
        self.rv1line = None
        self.rv2line = None
        self.astrometryellipses = None
        self.rv1errorlines = None
        self.rv2errorlines = None
        self.rv_window = None
        self.as_window = None
        self.rv_fig = None
        self.as_fig = None
        self.rv_ax = None
        self.as_ax = None
        self.didmcmc = False

        hcolor = '#3399ff'

        # PLOT FRAME #
        # initialize the plotting frame
        self.init_plots(plot_frame)
        # define the buttons
        tk.Button(plot_frame, text='Save RV figure', command=self.save_RV_plot, highlightbackground=hcolor).grid(row=1)
        tk.Button(plot_frame, text='Save AS figure', command=self.save_AS_plot, highlightbackground=hcolor).grid(row=3)
        tk.Button(plot_frame, text='Reset plots', command=lambda: self.init_plots(plot_frame),
                  highlightbackground=hcolor).grid(row=3, column=2, sticky=tk.SE)

        # GUESS FRAME #
        columns = 6
        paramcolumn = 1
        varycheckcolumn = 2
        transfercolumn = 3
        minresultcolumn = 4
        errorcolumn = 5

        # print the labels in the guess frame
        tk.Label(guess_frame, text='MODEL/GUESS PARAMETERS', font=('', 24)).grid(row=0, columnspan=columns, sticky=tk.N)
        tk.Label(guess_frame, text='Vary?').grid(row=1, column=varycheckcolumn)
        tk.Label(guess_frame, text='Result').grid(row=1, column=minresultcolumn)
        tk.Label(guess_frame, text='Error').grid(row=1, column=errorcolumn)

        self.param_name_vars = [tk.StringVar() for _ in range(11)]
        self.param_name_vars[0].set('e =')
        self.param_name_vars[1].set('i (deg) =')
        self.param_name_vars[2].set('omega (deg) =')
        self.param_name_vars[3].set('Omega (deg) =')
        self.param_name_vars[4].set('t0 (jd) =')
        self.param_name_vars[5].set('k1 (km/s) =')
        self.param_name_vars[6].set('k2 (km/s) =')
        self.param_name_vars[7].set('p (days) =')
        self.param_name_vars[8].set('gamma1 (km/s) =')
        self.param_name_vars[9].set('gamma2 (km/s) =')
        self.param_name_vars[10].set('d (pc) =')

        for i in range(len(self.param_name_vars)):
            tk.Label(guess_frame, textvariable=self.param_name_vars[i]).grid(row=(i + 2), sticky=tk.E)

        # initialize the entry variables
        self.guess_var_list = [tk.DoubleVar() for _ in range(11)]
        # define entry boxes
        self.guess_entry_list = [tk.Entry(guess_frame, textvariable=self.guess_var_list[i], width=10) for i in
                                 range(len(self.guess_var_list))]

        for i in range(len(self.guess_entry_list)):
            self.guess_entry_list[i].grid(row=(i + 2), column=paramcolumn)

        # define the vary state variables
        self.vary_var_list = [tk.BooleanVar() for _ in range(11)]

        # define checkbuttons for vary states
        self.vary_button_list = [tk.Checkbutton(guess_frame, var=self.vary_var_list[i]) for i in
                                 range(len(self.vary_var_list))]

        # put the checkbuttons in a nice grid
        for i in range(len(self.vary_button_list)):
            self.vary_button_list[i].grid(row=(i + 2), column=varycheckcolumn)

        # define the transfer buttons
        # for this semantic to work, we need to wrap the lambda function into another one, so that each command
        # references to its own number 'y', rather than the outer 'i' of the list comprehension
        self.transfer_button_list = [tk.Button(guess_frame, text='<-',
                                               command=(lambda y: (lambda: self.transfer(y)))(i),
                                               highlightbackground=hcolor).grid(row=(i + 2), column=transfercolumn)
                                     for i in range(len(self.guess_var_list))]

        # define the minimized parameter variables
        self.mininimzed_var_list = [tk.StringVar() for _ in range(11)]

        # define the labels the minimized parameters will go in
        self.minimized_label_list = [tk.Label(guess_frame, textvariable=self.mininimzed_var_list[i], width=10) for i in
                                     range(len(self.mininimzed_var_list))]

        # put the labels in a nice grid
        for i in range(len(self.minimized_label_list)):
            self.minimized_label_list[i].grid(row=(i + 2), column=minresultcolumn)

        # define the error variables
        self.error_var_list = [tk.StringVar() for _ in range(11)]

        # define the labels the errors will go in
        self.error_label_list = [tk.Label(guess_frame, textvariable=self.error_var_list[i], width=10) for i in
                                 range(len(self.error_var_list))]

        # put the labels in a nice grid
        for i in range(len(self.error_label_list)):
            self.error_label_list[i].grid(row=(i + 2), column=errorcolumn)

        # define the buttons in this frame
        tk.Button(guess_frame, text='plot model', command=self.plot_model,
                  highlightbackground=hcolor).grid(row=13)
        tk.Button(guess_frame, text='save guesses', command=self.save_guesses,
                  highlightbackground=hcolor).grid(row=13, column=1)
        tk.Button(guess_frame, text='Save parameters', command=self.save_params, highlightbackground=hcolor).grid(
            row=13, column=4, columnspan=2)

        # INFER FRAME #
        # define variables
        self.mprimary = tk.StringVar()
        self.msecondary = tk.StringVar()
        self.semimajor = tk.StringVar()

        # define labels
        tk.Label(infer_frame, text='INFERRED PARAMETERS', font=('', 24)).grid(columnspan=2, sticky=tk.N)
        tk.Label(infer_frame, text='Primary (M_sun) =').grid(row=1, sticky=tk.E)
        tk.Label(infer_frame, text='Secondary (M_sun) =').grid(row=2, sticky=tk.E)
        tk.Label(infer_frame, text='Physical semi major axis (R_sun) =').grid(row=3, sticky=tk.E)
        tk.Label(infer_frame, textvariable=self.mprimary).grid(row=1, column=1)
        tk.Label(infer_frame, textvariable=self.msecondary).grid(row=2, column=1)
        tk.Label(infer_frame, textvariable=self.semimajor).grid(row=3, column=1)

        # DATA FRAME #
        # define labels
        tk.Label(data_frame, text='DATA', font=('', 24)).grid(columnspan=4, sticky=tk.N)
        tk.Label(data_frame, text='Include?').grid(row=1, column=2)
        tk.Label(data_frame, text='Working directory').grid(row=2, sticky=tk.E)
        tk.Label(data_frame, text='Primary RV file').grid(row=3, sticky=tk.E)
        tk.Label(data_frame, text='Secondary RV file').grid(row=4, sticky=tk.E)
        tk.Label(data_frame, text='Astrometric data file').grid(row=5, sticky=tk.E)
        tk.Label(data_frame, text='Guess file').grid(row=6, sticky=tk.E)

        # define entries
        self.wd = tk.Entry(data_frame)
        self.rv1_file = tk.Entry(data_frame)
        self.rv2_file = tk.Entry(data_frame)
        self.as_file = tk.Entry(data_frame)
        self.guess_file = tk.Entry(data_frame)

        # put some mock values
        self.wd.insert(0, wd + '/')
        self.rv1_file.insert(0, 'primary_vels.txt')
        self.rv2_file.insert(0, 'secondary_vels.txt')
        self.as_file.insert(0, 'relative_astrometry.txt')
        self.guess_file.insert(0, 'guesses.txt')

        # put in a nice grid
        self.wd.grid(row=2, column=1)
        self.rv1_file.grid(row=3, column=1)
        self.rv2_file.grid(row=4, column=1)
        self.as_file.grid(row=5, column=1)
        self.guess_file.grid(row=6, column=1)

        # define inlcusion variables
        self.include_rv1 = tk.BooleanVar()
        self.include_rv2 = tk.BooleanVar()
        self.include_as = tk.BooleanVar()

        # assign to checkbuttons
        rv1check = tk.Checkbutton(data_frame, var=self.include_rv1)
        rv2check = tk.Checkbutton(data_frame, var=self.include_rv2)
        ascheck = tk.Checkbutton(data_frame, var=self.include_as)

        # put them in a nice grid
        rv1check.grid(row=3, column=2)
        rv2check.grid(row=4, column=2)
        ascheck.grid(row=5, column=2)

        # define the buttons
        tk.Button(data_frame, text='load data', command=self.load_data, highlightbackground=hcolor).grid(row=7,
                                                                                                         columnspan=2)

        tk.Button(data_frame, text='load guesses', command=self.set_guesses_from_file,
                  highlightbackground=hcolor).grid(row=8, columnspan=2)
        tk.Button(data_frame, text='plot data', command=self.plot_data, highlightbackground=hcolor).grid(row=9,
                                                                                                         columnspan=2)
        self.seppa = tk.BooleanVar()
        self.seppa.set(True)
        tk.Radiobutton(data_frame, text='Sep/PA', variable=self.seppa, value=True).grid(row=5, column=3)
        tk.Radiobutton(data_frame, text='E/N', variable=self.seppa, value=False).grid(row=5, column=4)

        # MINIMIZATION FRAME #
        # define variables
        self.mcmc = tk.BooleanVar()
        self.redchisq = tk.DoubleVar()
        self.dof = tk.IntVar()
        self.steps = tk.IntVar(value=1000)

        # define labels and buttons in a grid
        tk.Label(min_frame, text='MINIMIZATION', font=('', 24)).grid(row=0, columnspan=4)
        tk.Button(min_frame, text='Minimize model to data', command=self.minimize, highlightbackground=hcolor).grid(
            row=2, columnspan=4)
        tk.Label(min_frame, text='Do MCMC?:').grid(row=1, sticky=tk.E)
        mcmccheck = tk.Checkbutton(min_frame, var=self.mcmc)
        mcmccheck.grid(row=1, column=1)
        tk.Label(min_frame, text='# of samples:').grid(row=1, column=2, sticky=tk.E)
        self.stepsentry = tk.Entry(min_frame, textvariable=self.steps, width=5)
        self.stepsentry.grid(row=1, column=3)
        tk.Label(min_frame, text='Reduced Chi Squared').grid(row=3, columnspan=2, sticky=tk.E)
        tk.Label(min_frame, text='Degrees of freedom').grid(row=4, columnspan=2, sticky=tk.E)
        tk.Label(min_frame, textvariable=self.redchisq).grid(row=3, column=2, columnspan=2)
        tk.Label(min_frame, textvariable=self.dof).grid(row=4, column=2, columnspan=2)
        tk.Button(min_frame, text='make corner diagram', command=self.plot_corner_diagram,
                  highlightbackground=hcolor).grid(row=5, columnspan=4)

        # BACK TO ROOT FRAME #
        # display the root frame
        self.frame.pack(fill=tk.BOTH, expand=1)

    def make_plots(self):
        self.rv_fig = plt.figure()
        self.as_fig = plt.figure()
        self.rv_ax = self.rv_fig.add_subplot(111)
        self.as_ax = self.as_fig.add_subplot(111, aspect=1)
        spp.setup_rvax(self.rv_ax)
        spp.setup_asax(self.as_ax)
        self.rv_fig.tight_layout()
        self.as_fig.tight_layout()
        self.data_dict = None
        self.system = None

    def init_plots(self, plot_window):
        if self.rv_fig is not None:
            plt.close(self.rv_fig)
        if self.as_fig is not None:
            plt.close(self.as_fig)
        self.make_plots()
        # set the rv figure
        self.rv_window = FigureCanvasTkAgg(self.rv_fig, master=plot_window)
        self.rv_window.draw()
        self.rv_window.get_tk_widget().grid(row=0, sticky=tk.E)
        # set the as figure
        self.as_window = FigureCanvasTkAgg(self.as_fig, master=plot_window)
        self.as_window.draw()
        self.as_window.get_tk_widget().grid(row=2, sticky=tk.E)

    def transfer(self, varno):
        self.guess_var_list[varno].set(self.mininimzed_var_list[varno].get())

    def set_guesses_from_file(self):
        try:
            self.guess_dict = spl.guess_loader(self.wd.get(), self.guess_file.get())
            self.fill_guess_entries_from_dict()
        except IOError:
            print('cannot find your guess file!')
            self.guess_dict = None
        except (ValueError, KeyError):
            self.guess_dict = None
            print('some parameter has not been set properly')

    def set_guess_dict_from_entries(self):
        try:
            # here we must convert form list to dict, no way to write this faster
            self.guess_dict = dict()
            self.guess_dict = {'e': (self.guess_var_list[0].get(), self.vary_var_list[0].get()),
                               'i': (self.guess_var_list[1].get(), self.vary_var_list[1].get()),
                               'omega': (self.guess_var_list[2].get(), self.vary_var_list[2].get()),
                               'Omega': (self.guess_var_list[3].get(), self.vary_var_list[3].get()),
                               't0': (self.guess_var_list[4].get(), self.vary_var_list[4].get()),
                               'k1': (self.guess_var_list[5].get(), self.vary_var_list[5].get()),
                               'k2': (self.guess_var_list[6].get(), self.vary_var_list[6].get()),
                               'p': (self.guess_var_list[7].get(), self.vary_var_list[7].get()),
                               'gamma1': (self.guess_var_list[8].get(), self.vary_var_list[8].get()),
                               'gamma2': (self.guess_var_list[9].get(), self.vary_var_list[9].get()),
                               'd': (self.guess_var_list[10].get(), self.vary_var_list[10].get())}
        except ValueError or AttributeError:
            self.guess_dict = None
            self.system = None
            print('some parameter has not been set correctly!')

    def set_params_from_entries(self):
        try:
            # here we must convert form list to dict, no way to write this faster
            self.params = dict()
            self.params = {'e': self.guess_var_list[0].get(),
                           'i': self.guess_var_list[1].get(),
                           'omega': self.guess_var_list[2].get(),
                           'Omega': self.guess_var_list[3].get(),
                           't0': self.guess_var_list[4].get(),
                           'k1': self.guess_var_list[5].get(),
                           'k2': self.guess_var_list[6].get(),
                           'p': self.guess_var_list[7].get(),
                           'gamma1': self.guess_var_list[8].get(),
                           'gamma2': self.guess_var_list[9].get(),
                           'd': self.guess_var_list[10].get()}
        except ValueError or AttributeError:
            self.guess_dict = None
            self.system = None
            print('some parameter has not been set correctly!')

    def fill_guess_entries_from_dict(self):
        # here we must convert from dict to list, no way to write this faster
        self.guess_var_list[0].set(self.guess_dict['e'][0])
        self.guess_var_list[1].set(self.guess_dict['i'][0])
        self.guess_var_list[2].set(self.guess_dict['omega'][0])
        self.guess_var_list[3].set(self.guess_dict['Omega'][0])
        self.guess_var_list[4].set(self.guess_dict['t0'][0])
        self.guess_var_list[5].set(self.guess_dict['k1'][0])
        self.guess_var_list[6].set(self.guess_dict['k2'][0])
        self.guess_var_list[7].set(self.guess_dict['p'][0])
        self.guess_var_list[8].set(self.guess_dict['gamma1'][0])
        self.guess_var_list[9].set(self.guess_dict['gamma2'][0])
        self.guess_var_list[10].set(self.guess_dict['d'][0])

        self.vary_var_list[0].set(str(self.guess_dict['e'][1]))
        self.vary_var_list[1].set(str(self.guess_dict['i'][1]))
        self.vary_var_list[2].set(str(self.guess_dict['omega'][1]))
        self.vary_var_list[3].set(str(self.guess_dict['Omega'][1]))
        self.vary_var_list[4].set(str(self.guess_dict['t0'][1]))
        self.vary_var_list[5].set(str(self.guess_dict['k1'][1]))
        self.vary_var_list[6].set(str(self.guess_dict['k2'][1]))
        self.vary_var_list[7].set(str(self.guess_dict['p'][1]))
        self.vary_var_list[8].set(str(self.guess_dict['gamma1'][1]))
        self.vary_var_list[9].set(str(self.guess_dict['gamma2'][1]))
        self.vary_var_list[10].set(str(self.guess_dict['d'][1]))

    def set_system(self):
        try:
            self.set_params_from_entries()
            self.system = bsys.System(self.params)
            self.mprimary.set(self.system.primary_mass())
            self.msecondary.set(self.system.secondary_mass())
            self.semimajor.set(self.system.semimajor_axis())
        except ValueError as e:
            print(e)

    def plot_model(self):
        self.set_system()
        if self.system is not None:
            spp.plot_rv_curves(self.rv_ax, self.system)
            spp.plot_relative_orbit(self.as_ax, self.system)
            self.rv_fig.tight_layout()
            self.as_fig.tight_layout()
            self.rv_window.draw()
            self.as_window.draw()

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
            except OSError:
                print('Some file has not been found! Check your file paths!')
                self.data_dict = None
        else:
            print('no data files entered, or no data included!')
            self.data_dict = None

    def plot_data(self):
        self.load_data()
        self.set_system()
        try:
            spp.plot_as_data(self.as_ax, self.data_dict)
            self.as_fig.tight_layout()
            self.as_window.draw_idle()
        except (AttributeError, KeyError) as e:
            print(e)
        try:
            spp.plot_rv_data(self.rv_ax, self.data_dict, self.system)
            self.rv_fig.tight_layout()
            self.rv_window.draw_idle()
        except (KeyError, AttributeError) as e:
            print(e)
            print('set a model to plot RV data')

    def minimize(self):
        self.set_guess_dict_from_entries()
        self.load_data()
        if self.guess_dict is not None and self.data_dict is not None:
            # calculate best parameters
            self.minresult = spm.LMminimizer(self.guess_dict, self.data_dict, self.mcmc.get(), self.steps.get())
            if self.mcmc.get():
                self.didmcmc = True
            else:
                self.didmcmc = False
            pars = self.minresult.params
            # fill in the entries
            if self.guess_dict['e'][1]:
                self.mininimzed_var_list[0].set(np.round(pars['e'].value, 3))
                self.error_var_list[0].set(np.round(pars['e'].stderr, 3))
            if self.guess_dict['i'][1]:
                self.mininimzed_var_list[1].set(np.round(pars['i'].value, 3))
                self.error_var_list[1].set(np.round(pars['i'].stderr, 3))
            if self.guess_dict['omega'][1]:
                self.mininimzed_var_list[2].set(np.round(pars['omega'].value, 3))
                self.error_var_list[2].set(np.round(pars['omega'].stderr, 3))
            if self.guess_dict['Omega'][1]:
                self.mininimzed_var_list[3].set(np.round(pars['Omega'].value, 3))
                self.error_var_list[3].set(np.round(pars['Omega'].stderr, 3))
            if self.guess_dict['t0'][1]:
                self.mininimzed_var_list[4].set(np.round(pars['t0'].value, 3))
                self.error_var_list[4].set(np.round(pars['t0'].stderr, 3))
            if self.guess_dict['k1'][1]:
                self.mininimzed_var_list[5].set(np.round(pars['k1'].value, 3))
                self.error_var_list[5].set(np.round(pars['k1'].stderr, 3))
            if self.guess_dict['k2'][1]:
                self.mininimzed_var_list[6].set(np.round(pars['k2'].value, 3))
                self.error_var_list[6].set(np.round(pars['k2'].stderr, 3))
            if self.guess_dict['p'][1]:
                self.mininimzed_var_list[7].set(np.round(pars['p'].value, 3))
                self.error_var_list[7].set(np.round(pars['p'].stderr, 3))
            if self.guess_dict['gamma1'][1]:
                self.mininimzed_var_list[8].set(np.round(pars['gamma1'].value, 3))
                self.error_var_list[8].set(np.round(pars['gamma1'].stderr, 3))
            if self.guess_dict['gamma2'][1]:
                self.mininimzed_var_list[9].set(np.round(pars['gamma2'].value, 3))
                self.error_var_list[9].set(np.round(pars['gamma2'].stderr, 3))
            if self.guess_dict['d'][1]:
                self.mininimzed_var_list[10].set(np.round(pars['d'].value, 3))
                self.error_var_list[10].set(np.round(pars['d'].stderr, 3))

            self.redchisq.set(np.round(self.minresult.redchi, 4))
            self.dof.set(self.minresult.nfree)
            self.minimization_run_number += 1

    def plot_corner_diagram(self):
        if self.didmcmc:
            corner = spp.plot_corner_diagram(self.minresult)
            corner.savefig(self.wd.get() + 'corner{}.png'.format(self.minimization_run_number))
        else:
            print('do an mcmc minimization first!')

    def save_RV_plot(self):
        self.rv_fig.tight_layout()
        pls.plotsave(self.wd.get() + 'rv_plot', self.rv_fig, dpi=200)

    def save_AS_plot(self):
        self.as_fig.tight_layout()
        pls.plotsave(self.wd.get() + 'as_plot', self.as_fig, dpi=200)


wd = ''
try:
    wd = sys.argv[1]
except IndexError:
    pass

root = tk.Tk()
w, h = root.winfo_screenwidth(), root.winfo_screenheight()
root.geometry("{}x{}".format(w, h))
root.title('spinOSgui')
app = SpinOSApp(root, wd)
root.mainloop()
