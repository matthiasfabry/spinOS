from matplotlib.widgets import Slider, Button
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d, Axes3D
import numpy as np
import matplotlib.pyplot as plt


def plot_3d(orbit):
    phases = np.linspace(0, 1, 500)
    vrads = np.zeros(len(phases))
    thetas = np.zeros(len(phases))
    ecc_anoms = np.zeros(len(phases))
    for i in range(len(vrads)):
        vrads[i], thetas[i], ecc_anoms[i] = orbit.radial_velocity(phases[i], getAngles=True)
    xs = orbit.x(thetas)
    ys = orbit.y(thetas)
    zs = orbit.z(thetas)
    minz, maxz = min(zs), max(zs)
    fig = plt.figure()

    # make the orbital plane graph
    ax2 = fig.add_subplot(111, projection='3d')
    marker_point2, = ax2.plot([xs[0]], [ys[0]], [zs[0]], 'ro')
    x_line, = ax2.plot([0, xs[0]], [ys[0], ys[0]], [zs[0], zs[0]], 'r')
    y_line, = ax2.plot([xs[0], xs[0]], [0, ys[0]], [zs[0], zs[0]], 'r')
    z_line, = ax2.plot([xs[0], xs[0]], [ys[0], ys[0]], [0, zs[0]], 'r')
    # define the axes
    x_line = np.array([min(xs), max(xs)])
    y_line = np.array([min(ys), max(ys)])
    zeros = np.zeros(2)
    # plot x_y_plane
    x_axis, y_axis = np.meshgrid(x_line, y_line)
    x_y_plane = ax2.plot_surface(x_axis, y_axis, np.zeros([2, 2]), color='r', alpha=0.2)
    # plot x_z_plane
    x_axis, y_axis = np.meshgrid(x_line, zeros)
    x_z_plane = ax2.plot_surface(x_axis, y_axis, np.array([[minz, minz], [maxz, maxz]]), color='r', alpha=0.2)
    # plot y_z_plane
    x_axis, y_axis = np.meshgrid(zeros, y_line)
    y_z_plane = ax2.plot_surface(x_axis, y_axis, np.array([[minz, maxz], [minz, maxz]]), color='r', alpha=0.2)

    ax2.plot(xs, ys, zs, '.')
    ax2.plot([0, 0], [0, 0], [minz, maxz], label='line of sight', color='orange')
    arrow = Arrow3D([0, 0], [0, 0], [minz - 0.25 * (maxz - minz),
                                     maxz + 0.25 * (maxz - minz)],
                    mutation_scale=20,
                    lw=2, arrowstyle="-|>", label='line of sight', color='orange')
    ax2.add_artist(arrow)
    ax2.plot([orbit.x(orbit.omega), orbit.x(orbit.omega + np.pi)],
             [orbit.y(orbit.omega), orbit.y(orbit.omega + np.pi)],
             [0, 0], label='line of nodes')
    ax2.plot([orbit.x(0.), orbit.x(np.pi)],
             [orbit.y(0.), orbit.y(np.pi)],
             [orbit.z(0.), orbit.z(np.pi)], label='major axis')
    ax2.set_xlabel(r'$x/a$')
    ax2.set_ylabel(r'$y/a$')
    ax2.set_zlabel(r'$z/a$')
    ax2.legend()

    # make a Button to toggle axial planes
    # noinspection PyTypeChecker
    buttonaxes = plt.axes([0.85, 0.05, 0.15, 0.05])
    label = str('Toggle axes')
    but = Button(buttonaxes, label)

    class Toggle:
        def __init__(self):
            self._tog = True

        def toggle(self):
            self._tog = not self._tog

        def get_toggle_state(self):
            return self._tog

    tog = Toggle()

    def func(dummy):  # dummy argument is required by matplotlib
        x_y_plane.set_visible(not x_y_plane.get_visible())
        x_z_plane.set_visible(not x_z_plane.get_visible())
        y_z_plane.set_visible(not y_z_plane.get_visible())
        x_line.set_visible(not x_line.get_visible())
        y_line.set_visible(not y_line.get_visible())
        z_line.set_visible(not z_line.get_visible())
        if tog.get_toggle_state():
            ax2.set_axis_off()
        else:
            ax2.set_axis_on()
        tog.toggle()
        fig.canvas.draw_idle()

    but.on_clicked(func)

    # noinspection PyTypeChecker
    axphase = plt.axes([0.1, 0.1, 0.75, 0.03])
    sphase = Slider(axphase, 'phase', 0, 1, valinit=0, valstep=0.01)

    def update(dummy):  # dummy argument is required by matplotlib
        phase = sphase.val
        newvrad, newtheta, newecc_anom = orbit.radial_velocity(phase, getAngles=True)
        newx = orbit.x(newtheta)
        newy = orbit.y(newtheta)
        newz = orbit.z(newtheta)
        marker_point2.set_data_3d(newx, newy, newz)
        x_line.set_data_3d([0, newx], [newy, newy], [newz, newz])
        y_line.set_data_3d([newx, newx], [0, newy], [newz, newz])
        z_line.set_data_3d([newx, newx], [newy, newy], [0, newz])
        fig.canvas.draw_idle()

    sphase.on_changed(update)


# Modify the arrow class to draw the line of sight
class Arrow3D(FancyArrowPatch):
    def __init__(self, local_xs, local_ys, local_zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0, 0), (0, 0), *args, **kwargs)
        self._verts3d = local_xs, local_ys, local_zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        local_xs, local_ys, local_zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((local_xs[0], local_ys[0]), (local_xs[1], local_ys[1]))
        FancyArrowPatch.draw(self, renderer)
