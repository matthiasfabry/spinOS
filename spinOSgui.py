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
 JD(days) E_separation(mas) N_separation(mas) semimajor_ax_errorellipse(mas)
                                                            semiminor_ax_errorellipse(mas) angle_E_of_N_of_major_ax(deg)
eg:
 48000 -2.5 2.4 0.1 0.8 60
 48050 2.1 8.4 0.4 0.5 90
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
    13 Nov 2019

Version:
    1.3

Acknowledgements:
    This python3 implementation is heavily based on an earlier IDL implementation by Hugues Sana.
    We thank the authors of lmfit for the development of their package.

"""

import tkinter as tk

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

import binarySystem as bsys
import spinOSloader as spl
import spinOSplotter as spp
import spinOSminimizer as spm


def make_plots():
    fig1 = plt.figure(figsize=(10, 4.5), dpi=100)
    fig2 = plt.figure(figsize=(10, 4.5), dpi=100)
    ax1 = fig1.add_subplot(111)
    ax2 = fig2.add_subplot(111, aspect=1)
    spp.setup_rvax(ax1)
    spp.setup_asax(ax2)
    fig1.tight_layout()
    fig2.tight_layout()
    return fig1, fig2, ax1, ax2


class SpinOSApp:
    def __init__(self, master):
        # set the root frame
        self.frame = tk.Frame(master)

        # set the guess frame
        guess_frame = tk.Frame(self.frame)
        guess_frame.grid(row=1, column=0)
        # set the mass frame
        mass_frame = tk.Frame(self.frame)
        mass_frame.grid(row=2, column=0)
        # set the data frame
        data_frame = tk.Frame(self.frame)
        data_frame.grid(row=0, column=0)
        # set the minimization frame
        min_frame = tk.Frame(self.frame)
        min_frame.grid(row=3, column=0)
        # set the plotting frame
        plot_frame = tk.Frame(self.frame)
        plot_frame.grid(row=0, rowspan=4, column=1)

        # initialize some variables that will be set later
        self.guess_dict = None
        self.data_dict = None
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

        hcolor = '#3399ff'
        # PLOT FRAME #
        # initialize the plotting frame
        self.init_plots(plot_frame)
        # define the buttons
        tk.Button(plot_frame, text='Save RV figure', command=self.save_RV_plot, highlightbackground=hcolor).grid(row=1)
        tk.Button(plot_frame, text='Save AS figure', command=self.save_AS_plot, highlightbackground=hcolor).grid(row=3)
        tk.Button(plot_frame, text='Reset plots', command=lambda: self.init_plots(plot_frame),
                  highlightbackground=hcolor).grid(row=3, column=2)

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

        tk.Label(guess_frame, text='e =').grid(row=2, sticky=tk.E)
        tk.Label(guess_frame, text='i (deg)=').grid(row=3, sticky=tk.E)
        tk.Label(guess_frame, text='omega (deg)=').grid(row=4, sticky=tk.E)
        tk.Label(guess_frame, text='Omega (deg)=').grid(row=5, sticky=tk.E)
        tk.Label(guess_frame, text='t0 (hjd)=').grid(row=6, sticky=tk.E)
        tk.Label(guess_frame, text='k1 (km/s)=').grid(row=7, sticky=tk.E)
        tk.Label(guess_frame, text='k2 (km/s)=').grid(row=8, sticky=tk.E)
        tk.Label(guess_frame, text='p (days)=').grid(row=9, sticky=tk.E)
        tk.Label(guess_frame, text='gamma1 (km/s)=').grid(row=10, sticky=tk.E)
        tk.Label(guess_frame, text='gamma2 (km/s)=').grid(row=11, sticky=tk.E)
        tk.Label(guess_frame, text='d (pc)=').grid(row=12, sticky=tk.E)

        # initialize the entry variables
        self.entryvarlist = [tk.DoubleVar() for _ in range(11)]
        # define entry boxes
        self.entrylist = [tk.Entry(guess_frame, textvariable=self.entryvarlist[i], width=10) for i in
                          range(len(self.entryvarlist))]

        for i in range(len(self.entrylist)):
            self.entrylist[i].grid(row=(i + 2), column=paramcolumn)

        # define the vary state variables
        self.varyvarlist = [tk.BooleanVar() for _ in range(11)]

        # define checkbuttons for vary states
        self.varybutlist = [tk.Checkbutton(guess_frame, var=self.varyvarlist[i]) for i in range(len(self.varyvarlist))]

        # put the checkbuttons in a nice grid
        for i in range(len(self.varybutlist)):
            self.varybutlist[i].grid(row=(i + 2), column=varycheckcolumn)

        # define the transfer buttons
        # for this semantic to work, we need to wrap the lambda function into another one, so that each command
        # references to its own number 'y', rather than the outer 'i' of the list comprehension
        self.transferbutlist = [tk.Button(guess_frame, text='<-',
                                          command=(lambda y: (lambda: self.transfer(y)))(i),
                                          highlightbackground=hcolor).grid(row=(i + 2), column=transfercolumn)
                                for i in range(len(self.entryvarlist))]

        # define the minimized parameter variables
        self.minvarlist = [tk.StringVar() for _ in range(11)]

        # define the labels the minimized parameters will go in
        self.minlabellist = [tk.Label(guess_frame, textvariable=self.minvarlist[i], width=10) for i in
                             range(len(self.minvarlist))]

        # put the labels in a nice grid
        for i in range(len(self.minlabellist)):
            self.minlabellist[i].grid(row=(i + 2), column=minresultcolumn)

        # define the error variables
        self.errorvarlist = [tk.StringVar() for _ in range(11)]

        # define the labels the errors will go in
        self.errorlabellist = [tk.Label(guess_frame, textvariable=self.errorvarlist[i], width=10) for i in
                               range(len(self.errorvarlist))]

        # put the labels in a nice grid
        for i in range(len(self.errorlabellist)):
            self.errorlabellist[i].grid(row=(i + 2), column=errorcolumn)

        # define the buttons in this frame
        tk.Button(guess_frame, text='plot model', command=self.plot_guesses,
                  highlightbackground=hcolor).grid(row=13, columnspan=2)

        # MASS FRAME #
        # define variables
        self.mprimary = tk.StringVar()
        self.msecondary = tk.StringVar()

        # define labels
        tk.Label(mass_frame, text='INFERRED MASSES', font=('', 24)).grid(row=0, columnspan=2, sticky=tk.N)
        tk.Label(mass_frame, text='Primary (M_sun) =').grid(row=1, column=0, sticky=tk.E)
        tk.Label(mass_frame, text='Secondary (M_sun) =').grid(row=2, column=0, sticky=tk.E)
        tk.Label(mass_frame, textvariable=self.mprimary).grid(row=1, column=1)
        tk.Label(mass_frame, textvariable=self.msecondary).grid(row=2, column=1)

        # DATA FRAME #
        # define labels
        tk.Label(data_frame, text='DATA', font=('', 24)).grid(row=0, columnspan=3, sticky=tk.N)
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
        self.wd.insert(0, '9 Sgr/')
        self.rv1_file.insert(0, 'Ostarvels.txt')
        self.rv2_file.insert(0, 'WRstarvels.txt')
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

        # MINIMIZATION FRAME #
        # define variables
        self.mcmc = tk.BooleanVar()
        self.redchisq = tk.DoubleVar()
        self.dof = tk.IntVar()

        # define labels and buttons in a grid
        tk.Label(min_frame, text='MINIMIZATION', font=('', 24)).grid(row=0, columnspan=2)
        tk.Button(min_frame, text='Minimize model to data', command=self.minimize, highlightbackground=hcolor).grid(
            row=2)
        tk.Label(min_frame, text='Do MCMC?:').grid(row=1, column=0)
        mcmccheck = tk.Checkbutton(min_frame, var=self.mcmc)
        mcmccheck.grid(row=1, column=1)
        tk.Label(min_frame, text='Reduced Chi Squared').grid(row=3, column=0, sticky=tk.E)
        tk.Label(min_frame, text='Degrees of freedom').grid(row=4, column=0, sticky=tk.E)
        tk.Label(min_frame, textvariable=self.redchisq).grid(row=3, column=1)
        tk.Label(min_frame, textvariable=self.dof).grid(row=4, column=1)

        # BACK TO ROOT FRAME #
        # display the root frame
        self.frame.pack()

    def init_plots(self, plot_window):
        if self.rv_fig is not None:
            self.rv_fig.close()
        if self.as_fig is not None:
            self.as_fig.close()
        self.rv_fig, self.as_fig, self.rv_ax, self.as_ax = make_plots()
        # set the rv figure
        self.rv_window = FigureCanvasTkAgg(self.rv_fig, master=plot_window)
        self.rv_window.draw()
        self.rv_window.get_tk_widget().grid(row=0)
        # set the as figure
        self.as_window = FigureCanvasTkAgg(self.as_fig, master=plot_window)
        self.as_window.draw()
        self.as_window.get_tk_widget().grid(row=2)

    def transfer(self, varno):
        self.entryvarlist[varno].set(self.minvarlist[varno].get())

    def set_guesses_from_file(self):
        try:
            self.guess_dict = spl.guess_loader(self.wd.get(), self.guess_file.get())
            self.set_system()
            self.fill_guess_entries_from_dict()
        except IOError:
            print('cannot find your guess file!')
            self.guess_dict = None
            self.system = None
        except (ValueError, KeyError):
            self.guess_dict = None
            self.system = None
            print('some parameter has not been set properly')

    def set_guesses_from_entries(self):
        try:
            # here we must convert form list to dict, no way to write this faster
            self.guess_dict = dict()
            self.guess_dict['guesses'] = {'e': self.entryvarlist[0].get(),
                                          'i': self.entryvarlist[1].get() * np.pi / 180,
                                          'omega': self.entryvarlist[2].get() * np.pi / 180,
                                          'Omega': self.entryvarlist[3].get() * np.pi / 180,
                                          't0': self.entryvarlist[4].get(),
                                          'k1': self.entryvarlist[5].get(),
                                          'k2': self.entryvarlist[6].get(),
                                          'p': self.entryvarlist[7].get(),
                                          'gamma1': self.entryvarlist[8].get(),
                                          'gamma2': self.entryvarlist[9].get(),
                                          'd': self.entryvarlist[10].get()}
            self.guess_dict['varying'] = {'e': self.varyvarlist[0].get(),
                                          'i': self.varyvarlist[1].get(),
                                          'omega': self.varyvarlist[2].get(),
                                          'Omega': self.varyvarlist[3].get(),
                                          't0': self.varyvarlist[4].get(),
                                          'k1': self.varyvarlist[5].get(),
                                          'k2': self.varyvarlist[6].get(),
                                          'p': self.varyvarlist[7].get(),
                                          'gamma1': self.varyvarlist[8].get(),
                                          'gamma2': self.varyvarlist[9].get(),
                                          'd': self.varyvarlist[10].get()}
            self.set_system()
        except ValueError or AttributeError:
            self.guess_dict = None
            self.system = None
            print('some parameter has not been set correctly!')

    def fill_guess_entries_from_dict(self):
        # here we must convert from dict to list, no way to write this faster
        self.entryvarlist[0].set(self.guess_dict['guesses']['e'])
        self.entryvarlist[1].set(self.guess_dict['guesses']['i'] * 180 / np.pi)
        self.entryvarlist[2].set(self.guess_dict['guesses']['omega'] * 180 / np.pi)
        self.entryvarlist[3].set(self.guess_dict['guesses']['Omega'] * 180 / np.pi)
        self.entryvarlist[4].set(self.guess_dict['guesses']['t0'] % self.guess_dict['guesses']['p'])
        self.entryvarlist[5].set(self.guess_dict['guesses']['k1'])
        self.entryvarlist[6].set(self.guess_dict['guesses']['k2'])
        self.entryvarlist[7].set(self.guess_dict['guesses']['p'])
        self.entryvarlist[8].set(self.guess_dict['guesses']['gamma1'])
        self.entryvarlist[9].set(self.guess_dict['guesses']['gamma2'])
        self.entryvarlist[10].set(self.guess_dict['guesses']['d'])

        for entry in self.entryvarlist:
            entry.set(np.round(entry.get(), 3))

        self.varyvarlist[0].set(str(self.guess_dict['varying']['e']))
        self.varyvarlist[1].set(str(self.guess_dict['varying']['i']))
        self.varyvarlist[2].set(str(self.guess_dict['varying']['omega']))
        self.varyvarlist[3].set(str(self.guess_dict['varying']['Omega']))
        self.varyvarlist[4].set(str(self.guess_dict['varying']['t0']))
        self.varyvarlist[5].set(str(self.guess_dict['varying']['k1']))
        self.varyvarlist[6].set(str(self.guess_dict['varying']['k2']))
        self.varyvarlist[7].set(str(self.guess_dict['varying']['p']))
        self.varyvarlist[8].set(str(self.guess_dict['varying']['gamma1']))
        self.varyvarlist[9].set(str(self.guess_dict['varying']['gamma2']))
        self.varyvarlist[10].set(str(self.guess_dict['varying']['d']))

    def set_system(self):
        self.system = bsys.System(self.guess_dict['guesses'])
        self.mprimary.set(self.system.primary_mass())
        self.msecondary.set(self.system.secondary_mass())

    def plot_guesses(self):
        self.set_guesses_from_entries()
        if self.guess_dict is not None:
            spp.plot_rv_curves(self.rv_ax, self.system)
            spp.plot_relative_orbit(self.as_ax, self.system)
            self.rv_window.draw()
            self.as_window.draw()

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
                self.data_dict = spl.data_loader(self.wd.get(), filetypes, filenames)
            except OSError:
                print('Some file has not been found! Check your file paths!')
                self.data_dict = None
        else:
            print('no data files entered, or no data included!')
            self.data_dict = None

    def plot_data(self):
        if self.data_dict is None:
            self.load_data()
        if self.data_dict is not None:
            spp.plot_as_data(self.as_ax, self.data_dict)
            self.as_window.draw_idle()
        if self.data_dict is not None and self.system is not None:
            spp.plot_rv_data(self.rv_ax, self.data_dict, self.system)
        else:
            print('set a model first!')

    def minimize(self):
        self.set_guesses_from_entries()
        self.load_data()
        if self.guess_dict is not None and self.data_dict is not None:
            # calculate new best parameters
            minresult = spm.LMminimizer(self.guess_dict, self.data_dict, self.mcmc.get())
            pars = minresult.params
            # fill in the entries
            if self.guess_dict['varying']['e']:
                self.minvarlist[0].set(pars['e'].value)
                self.errorvarlist[0].set(pars['e'].stderr)
            if self.guess_dict['varying']['i']:
                self.minvarlist[1].set(pars['i'].value * 180 / np.pi)
                self.errorvarlist[1].set(pars['i'].stderr * 180 / np.pi)
            if self.guess_dict['varying']['omega']:
                self.minvarlist[2].set(pars['omega'].value * 180 / np.pi)
                self.errorvarlist[2].set(pars['omega'].stderr * 180 / np.pi)
            if self.guess_dict['varying']['Omega']:
                self.minvarlist[3].set(pars['Omega'].value * 180 / np.pi)
                self.errorvarlist[3].set(pars['Omega'].stderr * 180 / np.pi)
            if self.guess_dict['varying']['t0']:
                self.minvarlist[4].set(pars['t0'].value % pars['p'].value)
                self.errorvarlist[4].set(pars['t0'].stderr)
            if self.guess_dict['varying']['k1']:
                self.minvarlist[5].set(pars['k1'].value)
                self.errorvarlist[5].set(pars['k1'].stderr)
            if self.guess_dict['varying']['k2']:
                self.minvarlist[6].set(pars['k2'].value)
                self.errorvarlist[6].set(pars['k2'].stderr)
            if self.guess_dict['varying']['p']:
                self.minvarlist[7].set(pars['p'].value)
                self.errorvarlist[7].set(pars['p'].stderr)
            if self.guess_dict['varying']['gamma1']:
                self.minvarlist[8].set(pars['gamma1'].value)
                self.errorvarlist[8].set(pars['gamma1'].stderr)
            if self.guess_dict['varying']['gamma2']:
                self.minvarlist[9].set(pars['gamma2'].value)
                self.errorvarlist[9].set(pars['gamma2'].stderr)
            if self.guess_dict['varying']['d']:
                self.minvarlist[10].set(pars['d'].value)
                self.errorvarlist[10].set(pars['d'].stderr)

            for minvar in self.minvarlist:
                minvar.set(np.round(minvar.get(), 3))

            # set new guessdict and system masses
            self.set_guesses_from_entries()
            self.redchisq.set(minresult.redchi)
            self.dof.set(minresult.nfree)

    def save_RV_plot(self):
        self.rv_fig.savefig('rvplot')

    def save_AS_plot(self):
        self.as_fig.savefig('asplot')


root = tk.Tk()
w, h = root.winfo_screenwidth(), root.winfo_screenheight()
root.geometry("%dx%d" % (w, h))
root.title('spinOSgui')
app = SpinOSApp(root)
root.mainloop()
