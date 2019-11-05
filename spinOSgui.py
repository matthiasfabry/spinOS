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
    spinOSplotter.setup_relax(ax2)
    fig1.tight_layout()
    fig2.tight_layout()
    return fig1, fig2, ax1, ax2


class SpinOSApp:
    def __init__(self, master):
        # set the root frame
        self.frame = tk.Frame(master)

        # set the guess frame
        guess_frame = tk.Frame(self.frame, borderwidth=2)
        guess_frame.grid(row=0, column=0)
        # set the plotting frame
        plot_window = tk.Frame(self.frame)
        plot_window.grid(row=0, rowspan=3, column=1)
        # set the data frame
        data_frame = tk.Frame(self.frame, borderwidth=2)
        data_frame.grid(row=1, column=0)
        # set the button frame
        button_frame = tk.Frame(self.frame)
        button_frame.grid(row=2, column=0)

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

        # print the labels in the guess frame
        tk.Label(guess_frame, text='MODEL/GUESS PARAMETERS', font=('', 24)).grid(row=0, columnspan=2, sticky=tk.N)
        tk.Label(guess_frame, text='Vary?').grid(row=1, column=2)
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

        # initialize the entries values
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
        eentry.grid(row=2, column=1)
        ientry.grid(row=3, column=1)
        omegaentry.grid(row=4, column=1)
        Omegaentry.grid(row=5, column=1)
        t0entry.grid(row=6, column=1)
        k1entry.grid(row=7, column=1)
        k2entry.grid(row=8, column=1)
        pentry.grid(row=9, column=1)
        gamma1entry.grid(row=10, column=1)
        gamma2entry.grid(row=11, column=1)
        dentry.grid(row=12, column=1)
        search_intervalentry.grid(row=14, column=1)

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
        echeck.grid(row=2, column=2)
        icheck.grid(row=3, column=2)
        ocheck.grid(row=4, column=2)
        Ocheck.grid(row=5, column=2)
        t0check.grid(row=6, column=2)
        k1check.grid(row=7, column=2)
        k2check.grid(row=8, column=2)
        pcheck.grid(row=9, column=2)
        g1check.grid(row=10, column=2)
        g2check.grid(row=11, column=2)
        dcheck.grid(row=12, column=2)

        # The data frame
        datatitle = tk.Label(data_frame, text='DATA', font=('', 24))
        datatitle.grid(row=0, columnspan=2, sticky=tk.N)
        tk.Label(data_frame, text='Working directory').grid(row=1, sticky=tk.E)
        tk.Label(data_frame, text='Primary RV file').grid(row=2, sticky=tk.E)
        tk.Label(data_frame, text='Secondary RV file').grid(row=3, sticky=tk.E)
        tk.Label(data_frame, text='Astrometric data file').grid(row=4, sticky=tk.E)
        self.wd = tk.Entry(data_frame)
        self.wd.insert(0, 'testcase/')
        self.rv1_file = tk.Entry(data_frame)
        self.rv1_file.insert(0, 'Ostarvels.txt')
        self.rv2_file = tk.Entry(data_frame)
        self.rv2_file.insert(0, 'WRstarvels.txt')
        self.as_file = tk.Entry(data_frame)
        self.as_file.insert(0, 'relative_astrometry.txt')
        self.wd.grid(row=1, column=1)
        self.rv1_file.grid(row=2, column=1)
        self.rv2_file.grid(row=3, column=1)
        self.as_file.grid(row=4, column=1)

        # define buttons
        tk.Button(button_frame, text='load data', command=self.load_data, bg='blue').grid(row=0, column=0)
        tk.Button(button_frame, text='plot data', command=self.plot_data).grid(row=0, column=1)
        tk.Button(button_frame, text='plot model', command=self.plot_guesses).grid(row=1, column=0)
        tk.Button(button_frame, text='Minimize model to data', command=self.minimize).grid(row=1, column=1)
        tk.Button(button_frame, text='Reset plots', command=lambda: self.init_plots(plot_window)).grid(row=2,
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
        if self.rv1_file.get() is not '':
            filetypes.append('RV1file')
            filenames.append(self.rv1_file.get())
        if self.rv2_file.get() is not '':
            filetypes.append('RV2file')
            filenames.append(self.rv2_file.get())
        if self.as_file.get() is not '':
            filetypes.append('ASfile')
            filenames.append(self.as_file.get())
        if len(filenames) > 0:
            try:
                self.data_dict = spinOSloader.data_loader(self.wd.get(), filetypes, filenames)
            except OSError:
                print('Some file has not been found! Check your file paths!')
                self.data_dict = None
        else:
            print('no data files entered!')
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
            bestpars = spinOSminimizer.LMminimizer(self.guess_dict, float(self.search_interval.get()), self.data_dict)
            self.e.set(bestpars['e'])
            self.i.set(bestpars['i'])
            self.omega.set(bestpars['omega'])
            self.Omega.set(bestpars['Omega'])
            self.t0.set(bestpars['t0'])
            self.k1.set(bestpars['k1'])
            self.k2.set(bestpars['k2'])
            self.p.set(bestpars['p'])
            self.gamma1.set(bestpars['gamma1'])
            self.gamma2.set(bestpars['gamma2'])
            self.d.set(bestpars['d'])

    def save_RV_plot(self):
        self.rv_fig.savefig('rvplot')

    def save_AS_plot(self):
        self.as_fig.savefig('asplot')


root = tk.Tk()
w, h = root.winfo_screenwidth(), root.winfo_screenheight()
root.geometry("%dx%d+0+0" % (w, h))
root.title('spinOSgui')
app = SpinOSApp(root)
root.mainloop()
