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
 JD(days) E_separation(mas) N_separation(mas) major_ax_errorellipse(mas)
                                                            minor_ax_errorellipse(mas) angle_E_of_N_of_major_ax(deg)
eg:
 48000 -2.5 2.4 0.1 0.8 60
 48050 2.1 8.4 0.4 0.5 90
 etc...

Version:
    1.0

Author:
    Matthias Fabry
    Instituut voor Sterrekunde, KU Leuven, Belgium

Date:
    6 November 2019
"""

import tkinter as tk

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

import orbit
import spinOSloader
import spinOSplotter
import spinOSminimizer


def save_fig(fig, title):
    fig.savefig(title)


def make_plots():
    fig1 = plt.figure(figsize=(10, 4.5))
    fig2 = plt.figure(figsize=(10, 5))
    ax1 = fig1.add_subplot(111)
    ax2 = fig2.add_subplot(111, aspect=1)
    spinOSplotter.setup_rvax(ax1)
    spinOSplotter.setup_asax(ax2)
    fig1.tight_layout()
    fig2.tight_layout()
    return fig1, fig2, ax1, ax2


class SpinOSApp:
    def __init__(self, master):
        # set the root frame
        self.frame = tk.Frame(master)

        # set the guess frame
        guess_frame = tk.Frame(self.frame)
        guess_frame.grid(row=0, column=0)
        # set the mass frame
        mass_frame = tk.Frame(self.frame)
        mass_frame.grid(row=1, column=0)
        # set the data frame
        data_frame = tk.Frame(self.frame)
        data_frame.grid(row=2, column=0)
        # set the minimization frame
        min_frame = tk.Frame(self.frame)
        min_frame.grid(row=3, column=0)
        # set the button frame
        button_frame = tk.Frame(self.frame)
        button_frame.grid(row=4, column=0)
        # set the plotting frame
        plot_window = tk.Frame(self.frame)
        plot_window.grid(row=0, rowspan=5, column=1)

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

        # initialize the plotting frame
        self.init_plots(plot_window)
        # define the save buttons
        tk.Button(plot_window, text='Save RV figure', command=self.save_RV_plot).grid(row=1)
        tk.Button(plot_window, text='Save AS figure', command=self.save_AS_plot).grid(row=3)

        # GUESS FRAME #
        columns = 4
        paramcolumn = 1
        errorcolumn = 2
        varycheckcolumn = 3

        # print the labels in the guess frame
        tk.Label(guess_frame, text='MODEL/GUESS PARAMETERS', font=('', 24)).grid(row=0, columnspan=columns, sticky=tk.N)
        tk.Label(guess_frame, text='Vary?').grid(row=1, column=varycheckcolumn)
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
        tk.Label(guess_frame, text='Fractional search interval =').grid(row=14, sticky=tk.E)

        # initialize the entry variables
        self.e = tk.DoubleVar()
        self.i = tk.DoubleVar()
        self.omega = tk.DoubleVar()
        self.Omega = tk.DoubleVar()
        self.t0 = tk.DoubleVar()
        self.k1 = tk.DoubleVar()
        self.k2 = tk.DoubleVar()
        self.p = tk.DoubleVar()
        self.gamma1 = tk.DoubleVar()
        self.gamma2 = tk.DoubleVar()
        self.d = tk.DoubleVar()
        self.search_interval = tk.DoubleVar()

        # define entry boxes
        eentry = tk.Entry(guess_frame, textvariable=self.e)
        ientry = tk.Entry(guess_frame, textvariable=self.i)
        omegaentry = tk.Entry(guess_frame, textvariable=self.omega)
        Omegaentry = tk.Entry(guess_frame, textvariable=self.Omega)
        t0entry = tk.Entry(guess_frame, textvariable=self.t0)
        k1entry = tk.Entry(guess_frame, textvariable=self.k1)
        k2entry = tk.Entry(guess_frame, textvariable=self.k2)
        pentry = tk.Entry(guess_frame, textvariable=self.p)
        gamma1entry = tk.Entry(guess_frame, textvariable=self.gamma1)
        gamma2entry = tk.Entry(guess_frame, textvariable=self.gamma2)
        dentry = tk.Entry(guess_frame, textvariable=self.d)
        search_intervalentry = tk.Entry(guess_frame, textvariable=self.search_interval)

        # put the entries in a nice grid
        eentry.grid(row=2, column=paramcolumn)
        ientry.grid(row=3, column=paramcolumn)
        omegaentry.grid(row=4, column=paramcolumn)
        Omegaentry.grid(row=5, column=paramcolumn)
        t0entry.grid(row=6, column=paramcolumn)
        k1entry.grid(row=7, column=paramcolumn)
        k2entry.grid(row=8, column=paramcolumn)
        pentry.grid(row=9, column=paramcolumn)
        gamma1entry.grid(row=10, column=paramcolumn)
        gamma2entry.grid(row=11, column=paramcolumn)
        dentry.grid(row=12, column=paramcolumn)
        search_intervalentry.grid(row=14, column=paramcolumn)

        # initialize the entries with mock values
        self.e.set(0.89)
        self.i.set(120.)
        self.omega.set(225.)
        self.Omega.set(173.)
        self.t0.set(2600.)
        self.k1.set(29.)
        self.k2.set(78.)
        self.p.set(2850.)
        self.gamma1.set(1.1)
        self.gamma2.set(4.2)
        self.d.set(1640.)
        self.search_interval.set(0.5)

        # define the error variables
        self.eerror = tk.StringVar()
        self.ierror = tk.StringVar()
        self.omegaerror = tk.StringVar()
        self.Omegaerror = tk.StringVar()
        self.t0error = tk.StringVar()
        self.k1error = tk.StringVar()
        self.k2error = tk.StringVar()
        self.perror = tk.StringVar()
        self.gamma1error = tk.StringVar()
        self.gamma2error = tk.StringVar()
        self.derror = tk.StringVar()

        # define the labels the errors will go in
        eerrorlabel = tk.Label(guess_frame, textvariable=self.eerror)
        ierrorlabel = tk.Label(guess_frame, textvariable=self.ierror)
        omegaerrorlabel = tk.Label(guess_frame, textvariable=self.omegaerror)
        Omegaerrorlabel = tk.Label(guess_frame, textvariable=self.Omegaerror)
        t0errorlabel = tk.Label(guess_frame, textvariable=self.t0error)
        k1errorlabel = tk.Label(guess_frame, textvariable=self.k1error)
        k2errorlabel = tk.Label(guess_frame, textvariable=self.k2error)
        perrorlabel = tk.Label(guess_frame, textvariable=self.perror)
        gamma1errorlabel = tk.Label(guess_frame, textvariable=self.gamma1error)
        gamma2errorlabel = tk.Label(guess_frame, textvariable=self.gamma2error)
        derrorlabel = tk.Label(guess_frame, textvariable=self.derror)

        # put the labels in a nice grid
        eerrorlabel.grid(row=2, column=errorcolumn)
        ierrorlabel.grid(row=3, column=errorcolumn)
        omegaerrorlabel.grid(row=4, column=errorcolumn)
        Omegaerrorlabel.grid(row=5, column=errorcolumn)
        t0errorlabel.grid(row=6, column=errorcolumn)
        k1errorlabel.grid(row=7, column=errorcolumn)
        k2errorlabel.grid(row=8, column=errorcolumn)
        perrorlabel.grid(row=9, column=errorcolumn)
        gamma1errorlabel.grid(row=10, column=errorcolumn)
        gamma2errorlabel.grid(row=11, column=errorcolumn)
        derrorlabel.grid(row=12, column=errorcolumn)

        # define the vary state variables
        self.vary_e = tk.BooleanVar()
        self.vary_i = tk.BooleanVar()
        self.vary_omega = tk.BooleanVar()
        self.vary_Omega = tk.BooleanVar()
        self.vary_t0 = tk.BooleanVar()
        self.vary_k1 = tk.BooleanVar()
        self.vary_k2 = tk.BooleanVar()
        self.vary_p = tk.BooleanVar()
        self.vary_gamma1 = tk.BooleanVar()
        self.vary_gamma2 = tk.BooleanVar()
        self.vary_d = tk.BooleanVar()

        # define checkbuttons for vary states
        echeck = tk.Checkbutton(guess_frame, var=self.vary_e)
        icheck = tk.Checkbutton(guess_frame, var=self.vary_i)
        ocheck = tk.Checkbutton(guess_frame, var=self.vary_omega)
        Ocheck = tk.Checkbutton(guess_frame, var=self.vary_Omega)
        t0check = tk.Checkbutton(guess_frame, var=self.vary_t0)
        k1check = tk.Checkbutton(guess_frame, var=self.vary_k1)
        k2check = tk.Checkbutton(guess_frame, var=self.vary_k2)
        pcheck = tk.Checkbutton(guess_frame, var=self.vary_p)
        g1check = tk.Checkbutton(guess_frame, var=self.vary_gamma1)
        g2check = tk.Checkbutton(guess_frame, var=self.vary_gamma2)
        dcheck = tk.Checkbutton(guess_frame, var=self.vary_d)

        # put the checkbuttons in a nice grid
        echeck.grid(row=2, column=varycheckcolumn)
        icheck.grid(row=3, column=varycheckcolumn)
        ocheck.grid(row=4, column=varycheckcolumn)
        Ocheck.grid(row=5, column=varycheckcolumn)
        t0check.grid(row=6, column=varycheckcolumn)
        k1check.grid(row=7, column=varycheckcolumn)
        k2check.grid(row=8, column=varycheckcolumn)
        pcheck.grid(row=9, column=varycheckcolumn)
        g1check.grid(row=10, column=varycheckcolumn)
        g2check.grid(row=11, column=varycheckcolumn)
        dcheck.grid(row=12, column=varycheckcolumn)

        # MASS FRAME #
        # define labels
        tk.Label(mass_frame, text='INFERRED MASSES', font=('', 24)).grid(row=0, columnspan=2, sticky=tk.N)
        tk.Label(mass_frame, text='Primary (M_sun) =').grid(row=1, column=0, sticky=tk.E)
        tk.Label(mass_frame, text='Secondary (M_sun) =').grid(row=2, column=0, sticky=tk.E)

        self.mprimary = tk.StringVar()
        self.msecondary = tk.StringVar()

        tk.Label(mass_frame, textvariable=self.mprimary).grid(row=1, column=1)
        tk.Label(mass_frame, textvariable=self.msecondary).grid(row=2, column=1)

        # DATA FRAME #
        # define labels
        tk.Label(data_frame, text='DATA', font=('', 24)).grid(row=0, columnspan=2, sticky=tk.N)
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
        self.wd.insert(0, 'testcase/')
        self.rv1_file.insert(0, 'Ostarvels.txt')
        self.rv2_file.insert(0, 'WRstarvels.txt')
        self.as_file.insert(0, 'relative_astrometry.txt')
        self.guess_file.insert(0, 'guessfile.txt')

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

        # MINIMIZATION FRAME #
        # define some labels
        tk.Label(min_frame, text='MINIMIZATION STATISTICS', font=('', 24)).grid(row=0, columnspan=2)
        tk.Label(min_frame, text='Reduced Chi Squared').grid(row=1, column=0)
        tk.Label(min_frame, text='Degrees of freedom').grid(row=2, column=0)

        self.redchisq = tk.DoubleVar()
        self.dof = tk.IntVar()

        tk.Label(min_frame, textvariable=self.redchisq).grid(row=1, column=1)
        tk.Label(min_frame, textvariable=self.dof).grid(row=2, column=1)

        # BUTTON FRAME #
        # define buttons
        tk.Button(button_frame, text='load data', command=self.load_data).grid(row=0, column=0)
        tk.Button(button_frame, text='plot data', command=self.plot_data).grid(row=0, column=1)
        tk.Button(button_frame, text='load guesses', command=self.load_guesses_from_file).grid(row=1, column=0)
        tk.Button(button_frame, text='plot model', command=self.plot_guesses).grid(row=1, column=1)
        tk.Button(button_frame, text='Minimize model to data', command=self.minimize).grid(row=2, columnspan=2)
        tk.Button(button_frame, text='Reset plots', command=lambda: self.init_plots(plot_window)).grid(row=3,
                                                                                                       columnspan=2)

        # display the root frame
        self.frame.pack()

    def init_plots(self, plot_window):
        self.rv_fig, self.as_fig, self.rv_ax, self.as_ax = make_plots()
        # set the rv figure
        self.rv_window = FigureCanvasTkAgg(self.rv_fig, master=plot_window)
        self.rv_window.draw()
        self.rv_window.get_tk_widget().grid(row=0)
        # set the as figure
        self.as_window = FigureCanvasTkAgg(self.as_fig, master=plot_window)
        self.as_window.draw()
        self.as_window.get_tk_widget().grid(row=2)

    def load_guesses_from_file(self):
        try:
            self.guess_dict = spinOSloader.guess_loader(self.wd.get(), self.guess_file.get())
            self.set_guess_entries()
        except IOError:
            print('cannot find your guess file!')
            self.guess_dict = dict()

    def set_guess_entries(self):
        self.e.set(self.guess_dict['guesses']['e'])
        self.i.set(self.guess_dict['guesses']['i'] * 180 / np.pi)
        self.omega.set(self.guess_dict['guesses']['omega'] * 180 / np.pi)
        self.Omega.set(self.guess_dict['guesses']['Omega'] * 180 / np.pi)
        self.t0.set(self.guess_dict['guesses']['t0'] % self.guess_dict['guesses']['p'])
        self.k1.set(self.guess_dict['guesses']['k1'])
        self.k2.set(self.guess_dict['guesses']['k2'])
        self.p.set(self.guess_dict['guesses']['p'])
        self.gamma1.set(self.guess_dict['guesses']['gamma1'])
        self.gamma2.set(self.guess_dict['guesses']['gamma2'])
        self.d.set(self.guess_dict['guesses']['d'])

        self.vary_e.set(str(self.guess_dict['varying']['e']))
        self.vary_i.set(str(self.guess_dict['varying']['i']))
        self.vary_omega.set(str(self.guess_dict['varying']['omega']))
        self.vary_Omega.set(str(self.guess_dict['varying']['Omega']))
        self.vary_t0.set(str(self.guess_dict['varying']['t0']))
        self.vary_k1.set(str(self.guess_dict['varying']['k1']))
        self.vary_k2.set(str(self.guess_dict['varying']['k2']))
        self.vary_p.set(str(self.guess_dict['varying']['p']))
        self.vary_gamma1.set(str(self.guess_dict['varying']['gamma1']))
        self.vary_gamma2.set(str(self.guess_dict['varying']['gamma2']))
        self.vary_d.set(str(self.guess_dict['varying']['d']))

    def set_guesses(self):
        try:
            self.guess_dict = dict()
            self.guess_dict['guesses'] = {'e': self.e.get(), 'i': self.i.get() * np.pi / 180,
                                          'omega': self.omega.get() * np.pi / 180,
                                          'Omega': self.Omega.get() * np.pi / 180, 't0': self.t0.get(),
                                          'k1': self.k1.get(), 'k2': self.k2.get(), 'p': self.p.get(),
                                          'gamma1': self.gamma1.get(), 'gamma2': self.gamma2.get(), 'd': self.d.get()}
            self.guess_dict['varying'] = {'e': self.vary_e.get(), 'i': self.vary_i.get(),
                                          'omega': self.vary_omega.get(),
                                          'Omega': self.vary_Omega.get(), 't0': self.vary_t0.get(),
                                          'k1': self.vary_k1.get(), 'k2': self.vary_k2.get(), 'p': self.vary_p.get(),
                                          'gamma1': self.vary_gamma1.get(), 'gamma2': self.vary_gamma2.get(),
                                          'd': self.vary_d.get()}
            self.system = orbit.System(self.guess_dict['guesses'])
            self.mprimary.set(self.system.primary_mass())
            self.msecondary.set(self.system.secondary_mass())
        except ValueError or AttributeError:
            self.guess_dict = None
            self.system = None
            print('some parameter has not been set correctly!')

    def plot_guesses(self):
        self.set_guesses()
        if self.guess_dict is not None:
            spinOSplotter.plot_rv_curves(self.rv_ax, self.system)
            spinOSplotter.plot_relative_orbit(self.as_ax, self.system)
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
                self.data_dict = spinOSloader.data_loader(self.wd.get(), filetypes, filenames)
            except OSError:
                print('Some file has not been found! Check your file paths!')
                self.data_dict = None
        else:
            print('no data files entered, or no data included!')
            self.data_dict = None

    def plot_data(self):
        self.set_guesses()
        if self.data_dict is None:
            self.load_data()
        if self.data_dict is not None and self.system is not None:
            spinOSplotter.plot_data(self.rv_ax, self.as_ax, self.data_dict, self.system)
            self.rv_window.draw_idle()
            self.as_window.draw()

    def minimize(self):
        self.set_guesses()
        self.load_data()
        if self.guess_dict is not None and self.data_dict is not None:
            # calculate new best parameters
            minimizationresult = spinOSminimizer.LMminimizer(self.guess_dict, self.data_dict,
                                                             float(self.search_interval.get()))
            # fill in best pars
            self.guess_dict['guesses'] = minimizationresult.params.valuesdict()
            # fill in the entries
            self.set_guess_entries()
            # set new guessdict and system masses
            self.set_guesses()
            self.redchisq.set(minimizationresult.redchi)
            self.dof.set(minimizationresult.nfree)

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
