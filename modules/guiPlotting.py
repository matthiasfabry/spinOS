import tkinter as tk

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.collections import EllipseCollection

import constants as cst


class Plotting:

    def __init__(self, gui, tab):
        self.gui = gui
        plt_frame_top = tk.Frame(tab)
        plt_frame = tk.Frame(plt_frame_top)
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
        self.phase_button = tk.Checkbutton(plt_frame, var=self.do_phasedot, command=self.gui.toggle_dot,
                                           state=tk.DISABLED)
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
        self.pphase_but = tk.Radiobutton(plt_frame, text='phase', command=self.gui.toggle_phase_time,
                                         variable=self.plot_vs_phase, value=True, state=tk.DISABLED)
        self.ptime_but = tk.Radiobutton(plt_frame, text='time', command=self.gui.toggle_phase_time,
                                        variable=self.plot_vs_phase, value=False, state=tk.DISABLED)
        self.pphase_but.grid(row=5, column=2)
        self.ptime_but.grid(row=5, column=3)
        self.modelwidgets = {self.plot_asmodel_label, self.plot_asmodel_button,
                             self.plot_rv1model_button, self.plot_rv2model_button, self.plot_rv1model_label,
                             self.plot_rv2model_label, self.plot_semimajor_button, self.plot_semimajor_label,
                             self.plot_nodeline_button, self.plot_nodeline_label, self.plot_peri_label,
                             self.plot_peri_button, self.as_dist_button, self.as_dist_label,
                             self.pphase_but, self.ptime_but}
        plt_frame.pack()

        refreshframe = tk.Frame(plt_frame_top)
        tk.Button(refreshframe, text='Refresh Plots', width=20, height=2, command=self.gui.update,
                  highlightbackground=cst.HCOLOR).pack()
        refreshframe.pack(pady=10)

        settings_frame = tk.Frame(plt_frame_top)
        entrycol = 1
        tk.Label(settings_frame, text='PLOT SETTINGS', font=('', cst.TITLESIZE, 'underline')).grid(columnspan=2)
        tk.Label(settings_frame, text='Axis label size').grid(row=1)
        self.axeslabelsize = tk.DoubleVar(value=plt.rcParams['font.size'])
        tk.Entry(settings_frame, textvariable=self.axeslabelsize).grid(row=1, column=entrycol)
        tk.Label(settings_frame, text='Tick label size').grid(row=2)
        self.ticklabelsize = tk.DoubleVar(value=plt.rcParams['font.size'])
        tk.Entry(settings_frame, textvariable=self.ticklabelsize).grid(row=2, column=entrycol)

        settings_frame.pack(pady=10)

        refreshframe = tk.Frame(plt_frame_top)
        tk.Button(refreshframe, text='Refresh Plots', width=20, height=2, command=self.gui.update,
                  highlightbackground=cst.HCOLOR).pack()
        refreshframe.pack(pady=10)

        plt_frame_top.place(relx=.5, rely=0, anchor="n")

    def update_plots(self):
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
        self.setup_rv_ax()
        self.setup_as_ax()
        self.relim_plots()
        self.rv_fig.canvas.draw()
        self.as_fig.canvas.draw()

    def init_plots(self):
        """
        sets up the plot windows
        """
        plt.rc('text', usetex=True)
        plt.rc('font', size=20)
        plt.rc('font', family='serif')

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
        move_figure(self.rv_fig, int(0.35 * self.gui.w) + 10, 0)
        move_figure(self.as_fig, int(0.35 * self.gui.w) + 10, int(self.gui.h / 2) + 10)
        self.rv_ax = self.rv_fig.add_subplot(111)
        self.as_ax = self.as_fig.add_subplot(111, aspect=1)
        self.rv_ax.axhline(linestyle=':', color='black')
        self.rv_ax.grid()
        self.setup_rv_ax()
        self.as_ax.invert_xaxis()
        self.as_ax.axhline(linestyle=':', color='black')
        self.as_ax.axvline(linestyle=':', color='black')
        self.as_ax.grid()
        self.setup_as_ax()
        self.rv_fig.tight_layout()
        self.as_fig.tight_layout()
        plt.ion()
        plt.show()

    def setup_rv_ax(self):
        self.rv_ax.tick_params(axis='both', which='major', direction='in', labelsize=self.ticklabelsize.get())
        if self.plot_vs_phase.get():
            self.rv_ax.set_xlabel(cst.PHASE_STR, fontdict={'size': self.axeslabelsize.get()})
        else:
            self.rv_ax.set_xlabel(cst.TIME_STR, fontdict={'size': self.axeslabelsize.get()})
        self.rv_ax.set_ylabel(r'$RV$ (km s$^{-1}$)', fontdict={'size': self.axeslabelsize.get()})
        self.rv_ax.set_xlim((-0.18, 1.18))
        self.rv_ax.set_ylim((-50, 50))

    def setup_as_ax(self):
        self.as_ax.tick_params(axis='both', which='major', direction='in', labelsize=self.ticklabelsize.get())
        self.as_ax.set_xlabel(r'East (mas)', fontdict={'size': self.axeslabelsize.get()})
        self.as_ax.set_ylabel(r'North (mas)', fontdict={'size': self.axeslabelsize.get()})
        self.as_ax.set_xlim((-10, 10))
        self.as_ax.set_ylim((-10, 10))
        # asax.xaxis.set_major_locator(MultipleLocator(2.5))
        # asax.yaxis.set_major_locator(MultipleLocator(2.5))

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
        if 'RV1' not in self.gui.data_dict:
            return
        if self.rv1data_line is not None:
            self.rv1data_line.remove()
            self.rv1data_line = None
        if self.plot_vs_phase.get():
            phases, rv, err = self.gui.system.create_phase_extended_RV(self.gui.data_dict['RV1'], 0.15)
            self.rv1data_line = self.rv_ax.errorbar(phases, rv, yerr=err, ls='', capsize=0.1, marker='o',
                                                    ms=5, color='b')
        else:
            self.rv1data_line = self.rv_ax.errorbar(self.gui.data_dict['RV1']['hjds'], self.gui.data_dict['RV1']['RVs'],
                                                    yerr=self.gui.data_dict['RV1']['errors'], ls='', capsize=0.1,
                                                    marker='o', ms=5, color='b')

    def plot_rv2_data(self):
        """
        plot the rv2 data
        """
        if 'RV2' not in self.gui.data_dict:
            return
        if self.rv2data_line is not None:
            self.rv2data_line.remove()
            self.rv2data_line = None
        if self.plot_vs_phase.get():
            phases, rv, err = self.gui.system.create_phase_extended_RV(self.gui.data_dict['RV2'], 0.15)
            self.rv2data_line = self.rv_ax.errorbar(phases, rv, yerr=err, ls='', capsize=0.1, marker='o',
                                                    ms=5, color='r')
        else:
            self.rv2data_line = self.rv_ax.errorbar(self.gui.data_dict['RV2']['hjds'], self.gui.data_dict['RV2']['RVs'],
                                                    yerr=self.gui.data_dict['RV2']['errors'], ls='', capsize=0.1,
                                                    marker='o', ms=5, color='r')

    def plot_as_data(self):
        """
        plot the as data
        """
        if 'AS' not in self.gui.data_dict:
            return
        data = self.gui.data_dict['AS']
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
        if 'AS' not in self.gui.data_dict:
            return
        data = self.gui.data_dict['AS']
        if self.as_dist_lines is not None:
            for line in self.as_dist_lines:
                line.remove()
            self.as_dist_lines = None
        self.as_dist_lines = list()
        for i in range(len(data['hjds'])):
            self.as_dist_lines.append(self.as_ax.plot(
                (data['easts'][i], self.gui.system.relative.east_of_hjd(data['hjds'][i])),
                (data['norths'][i], self.gui.system.relative.north_of_hjd(data['hjds'][i])),
                c=(0.75, 0.25, 0.0, 0.9))[0])

    def plot_rv1_curve(self):
        """
        the the rv1 model curve
        """
        if self.plot_vs_phase.get():
            phases = np.linspace(-0.15, 1.15, num=150)
            vrads1 = self.gui.system.primary.radial_velocity_of_phases(phases)
            if self.rv1_line is None:
                self.rv1_line, = self.rv_ax.plot(phases, vrads1, label=r'primary', color='b', ls='--')
            else:
                self.rv1_line.set_xdata(phases)
                self.rv1_line.set_ydata(vrads1)
        else:
            m = np.infty
            mm = -np.infty
            if self.gui.include_rv1.get():
                m = min(m, min(self.gui.data_dict['RV1']['hjds']))
                mm = max(mm, max(self.gui.data_dict['RV1']['hjds']))
            if self.gui.include_rv2.get():
                m = min(m, min(self.gui.data_dict['RV2']['hjds']))
                mm = max(mm, max(self.gui.data_dict['RV2']['hjds']))
            times = np.linspace(m - 0.01 * (mm - m), m - 0.01 * (mm - m) + self.gui.system.p, endpoint=False, num=100)
            rvs = self.gui.system.primary.radial_velocity_of_phases(self.gui.system.phase_of_hjds(times))
            times, rvs = self.gui.system.extend_rvs_until_time(times, rvs, mm)
            if self.rv1_line is None:
                self.rv1_line, = self.rv_ax.plot(times, rvs, label=r'primary', color='b', ls='--')
            else:
                self.rv1_line.set_xdata(times)
                self.rv1_line.set_ydata(rvs)

    def plot_rv2_curve(self):
        """
        plot the rv2 model curve
        """
        if self.plot_vs_phase.get():
            phases = np.linspace(-0.15, 1.15, num=150)
            vrads1 = self.gui.system.secondary.radial_velocity_of_phases(phases)
            if self.rv2_line is None:
                self.rv2_line, = self.rv_ax.plot(phases, vrads1, label=r'secondary', color='r', ls='--')
            else:
                self.rv2_line.set_xdata(phases)
                self.rv2_line.set_ydata(vrads1)
        else:
            m = np.infty
            mm = -np.infty
            if self.gui.include_rv2.get():
                m = min(m, min(self.gui.data_dict['RV2']['hjds']))
                mm = max(mm, max(self.gui.data_dict['RV2']['hjds']))
            if self.gui.include_rv1.get():
                m = min(m, min(self.gui.data_dict['RV1']['hjds']))
                mm = max(mm, max(self.gui.data_dict['RV1']['hjds']))
            times = np.linspace(m, m + self.gui.system.p, num=100)
            rvs = self.gui.system.secondary.radial_velocity_of_phases(self.gui.system.phase_of_hjds(times))
            times, rvs = self.gui.system.extend_rvs_until_time(times, rvs, mm)
            if self.rv2_line is None:
                self.rv2_line, = self.rv_ax.plot(times, rvs, label=r'secondary', color='r', ls='--')
            else:
                self.rv2_line.set_xdata(times)
                self.rv2_line.set_ydata(rvs)

    def plot_relative_orbit(self):
        """
        (re)plot the relative astrometric orbit
        """
        ecc_anoms = np.linspace(0, 2 * np.pi, 200)
        norths = self.gui.system.relative.north_of_ecc(ecc_anoms)
        easts = self.gui.system.relative.east_of_ecc(ecc_anoms)
        if self.as_line is None:
            self.as_line, = self.as_ax.plot(easts, norths, label='relative orbit', color='k')
        else:
            self.as_line.set_xdata(easts)
            self.as_line.set_ydata(norths)

    def plot_node_line(self):
        """
        (re)plot the astrometric node line
        """
        system = self.gui.system.relative
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
        (re)plot the astrometric periastron point
        """
        system = self.gui.system.relative
        if self.peri_dot is None:
            self.peri_dot, = self.as_ax.plot([system.east_of_ecc(0)], [system.north_of_ecc(0)], color='b', marker='s',
                                             ls='', fillstyle='full', label='periastron', markersize=8)
        else:
            self.peri_dot.set_xdata(system.east_of_ecc(0))
            self.peri_dot.set_ydata(system.north_of_ecc(0))

    def plot_semimajor_axis(self):
        """
        (re)plot the astrometric semimajor axis
        """
        system = self.gui.system.relative
        if self.semi_major is None:
            self.semi_major, = self.as_ax.plot([system.east_of_true(0), system.east_of_true(np.pi)],
                                               [system.north_of_true(0), system.north_of_true(np.pi)],
                                               color='0.3', ls='dotted', label='semi-major axis')
        else:
            self.semi_major.set_xdata([system.east_of_true(0), system.east_of_true(np.pi)])
            self.semi_major.set_ydata([system.north_of_true(0), system.north_of_true(np.pi)])

    def plot_dots(self):
        """
        (re)plot diamond shapes at the specified phase
        """
        if self.rv1_dot is not None:
            self.rv1_dot.remove()
            self.rv1_dot = None
        if self.do_modelrv1.get() or self.do_datarv1.get():
            rv1 = self.gui.system.primary.radial_velocity_of_phase(self.phase.get())
            self.rv1_dot = self.rv_ax.scatter(self.phase.get(), rv1, s=100, color='b', marker='D',
                                              label=np.round(rv1, 2))
        if self.rv2_dot is not None:
            self.rv2_dot.remove()
            self.rv2_dot = None
        if self.do_modelrv2.get() or self.do_datarv2.get():
            rv2 = self.gui.system.secondary.radial_velocity_of_phase(self.phase.get())
            self.rv2_dot = self.rv_ax.scatter(self.phase.get(), rv2, s=100, color='r', marker='D',
                                              label=np.round(rv2, 2))
        if self.as_dot is not None:
            self.as_dot.remove()
            self.as_dot = None
        if self.do_modelas.get() or self.do_dataas.get():
            N = self.gui.system.relative.north_of_ph(self.phase.get())
            E = self.gui.system.relative.east_of_ph(self.phase.get())
            self.as_dot = self.as_ax.scatter(E, N, s=100, color='r', marker='x',
                                             label='{}E/{}N'.format(np.round(E, 2), np.round(N, 2)))

    def plot_legends(self):
        """
        (re)plot the legends
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
