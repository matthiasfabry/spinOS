import tkinter as tk
import spinOSplotter
import orbit
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
import numpy as np


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
        frame = tk.Frame(master)

        # set the guess frame
        guess_frame = tk.Frame(frame, borderwidth=2)
        guess_frame.grid(row=0, column=0)
        # set the plotting frame
        plot_window = tk.Frame(frame)
        plot_window.grid(row=0, rowspan=3, column=1)
        # set the data frame
        data_frame = tk.Frame(frame, borderwidth=2)
        data_frame.grid(row=1, column=0)
        # set the button frame
        button_frame = tk.Frame(frame)
        button_frame.grid(row=2, column=0)

        # initialize the plotting frame
        self.rvfig, self.asfig, self.rvax, self.relax = make_plots()
        # set the rv figure
        self.rv_plot = FigureCanvasTkAgg(self.rvfig, master=plot_window)
        self.rv_plot.draw()
        self.rv_plot.get_tk_widget().grid(row=0)
        # define the save button
        tk.Button(plot_window, text='Save RV figure', command=lambda: save_fig(self.rvfig, 'rvplot')).grid(row=1)
        # set the as figure
        self.astrometry_plot = FigureCanvasTkAgg(self.asfig, master=plot_window)
        self.astrometry_plot.draw()
        self.astrometry_plot.get_tk_widget().grid(row=2)
        # define the save button
        tk.Button(plot_window, text='Save AS figure', command=lambda: save_fig(self.asfig, 'asplot')).grid(row=3)

        # print the labels in the guess frame
        tk.Label(guess_frame, text='MODEL/GUESS PARAMETERS', font=('', 24)).grid(row=0, columnspan=2, sticky=tk.N)
        tk.Label(guess_frame, text='e =').grid(row=1, sticky=tk.E)
        tk.Label(guess_frame, text='i (deg)=').grid(row=2, sticky=tk.E)
        tk.Label(guess_frame, text='omega (deg)=').grid(row=3, sticky=tk.E)
        tk.Label(guess_frame, text='Omega (deg)=').grid(row=4, sticky=tk.E)
        tk.Label(guess_frame, text='t0 (hjd)=').grid(row=5, sticky=tk.E)
        tk.Label(guess_frame, text='k1 (km/s)=').grid(row=6, sticky=tk.E)
        tk.Label(guess_frame, text='k2 (km/s)=').grid(row=7, sticky=tk.E)
        tk.Label(guess_frame, text='p (days)=').grid(row=8, sticky=tk.E)
        tk.Label(guess_frame, text='gamma1 (km/s)=').grid(row=9, sticky=tk.E)
        tk.Label(guess_frame, text='gamma2 (km/s)=').grid(row=10, sticky=tk.E)
        tk.Label(guess_frame, text='d (pc)=').grid(row=11, sticky=tk.E)
        # initialize the entries with some mock values
        self.e = tk.Entry(guess_frame)
        self.e.insert(0, 0.5)
        self.e.grid(row=1, column=1)
        self.i = tk.Entry(guess_frame)
        self.i.insert(0, 20)
        self.i.grid(row=2, column=1)
        self.omega = tk.Entry(guess_frame)
        self.omega.insert(0, 20)
        self.omega.grid(row=3, column=1)
        self.Omega = tk.Entry(guess_frame)
        self.Omega.insert(0, 20)
        self.Omega.grid(row=4, column=1)
        self.t0 = tk.Entry(guess_frame)
        self.t0.insert(0, 0)
        self.t0.grid(row=5, column=1)
        self.k1 = tk.Entry(guess_frame)
        self.k1.insert(0, 10)
        self.k1.grid(row=6, column=1)
        self.k2 = tk.Entry(guess_frame)
        self.k2.insert(0, 10)
        self.k2.grid(row=7, column=1)
        self.p = tk.Entry(guess_frame)
        self.p.insert(0, 500)
        self.p.grid(row=8, column=1)
        self.gamma1 = tk.Entry(guess_frame)
        self.gamma1.insert(0, 5)
        self.gamma1.grid(row=9, column=1)
        self.gamma2 = tk.Entry(guess_frame)
        self.gamma2.insert(0, 5)
        self.gamma2.grid(row=10, column=1)
        self.d = tk.Entry(guess_frame)
        self.d.insert(0, 1000)
        self.d.grid(row=11, column=1)
        datatitle = tk.Label(data_frame, text='DATA', font=('', 24))
        datatitle.grid(row=0, columnspan=2, sticky=tk.N)
        tk.Label(data_frame, text='Primary RV file').grid(row=1, sticky=tk.E)
        tk.Label(data_frame, text='Secondary RV file').grid(row=2, sticky=tk.E)
        tk.Label(data_frame, text='Astrometric data file').grid(row=3, sticky=tk.E)
        self.rv1_file = tk.Entry(data_frame)
        self.rv2_file = tk.Entry(data_frame)
        self.as_file = tk.Entry(data_frame)
        self.rv1_file.grid(row=1, column=1)
        self.rv2_file.grid(row=2, column=1)
        self.as_file.grid(row=3, column=1)

        # define buttons
        tk.Button(button_frame, text='load data', command=self.load_data).grid(row=0, column=0)
        tk.Button(button_frame, text='plot data', command=self.plot_data).grid(row=0, column=1)
        tk.Button(button_frame, text='plot model', command=self.plot_guesses).grid(row=1, column=0)
        tk.Button(button_frame, text='Minimize model to data', command=self.minimize).grid(row=1, column=1)

        # display the root frame
        frame.pack()

    def plot_guesses(self):
        try:
            system = orbit.System(
                {'e': np.float64(self.e.get()), 'i': np.float64(self.i.get()) * np.pi / 180,
                 'omega': np.float64(self.omega.get()) * np.pi / 180,
                 'Omega': np.float64(self.Omega.get()) * np.pi / 180, 't0': np.float64(self.t0.get()),
                 'k1': np.float64(self.k1.get()),
                 'k2': np.float64(self.k2.get()), 'p': np.float64(self.p.get()),
                 'gamma1': np.float64(self.gamma1.get()),
                 'gamma2': np.float64(self.gamma2.get()), 'd': np.float64(self.d.get())})
            spinOSplotter.plot_rv_curves(self.rvax, system)
            spinOSplotter.plot_relative_orbit(self.relax, system)
            self.rv_plot.draw()
            self.astrometry_plot.draw()
        except AttributeError or ValueError:
            print('some parameter has not been set!')

    def load_data(self):
        pass

    def plot_data(self):
        pass

    def minimize(self):
        pass

    def save_RV_plot(self):
        self.rvfig.savefig('rvplot')

    def save_AS_plot(self):
        self.asfig.savefig('asplot')


root = tk.Tk()
w, h = root.winfo_screenwidth(), root.winfo_screenheight()
root.geometry("%dx%d+0+0" % (w, h))
root.title('spinOSgui')
app = SpinOSApp(root)
root.mainloop()
