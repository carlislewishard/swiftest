import swiftestio as swio
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import subprocess

radscale = 20
RMars = 3389500.0
xmin = 1.0
xmax = 10.0
ymin = 1e-6
ymax = 1.0

class AnimatedScatter(object):
    """An animated scatter plot using matplotlib.animations.FuncAnimation."""
    def __init__(self, ds, config):
        #outf = 'aescatter.mp4'

        frame = 0
        nframes = ds['time'].size
        self.ds = ds
        self.config = config
        self.ds['mass'] = self.ds['mass'] / config['GU']
        self.ds['radmarker'] = self.ds['radius'].fillna(0)
        self.ds['radmarker'] = self.ds['radmarker'] / self.ds['radmarker'].max() * radscale


        self.stream = self.data_stream(frame)
        # Setup the figure and axes...
        self.fig, self.ax = plt.subplots(figsize=(8,4.5))
        # Then setup FuncAnimation.
        self.ani = animation.FuncAnimation(self.fig, self.update, interval=1, frames=nframes,
                                          init_func=self.setup_plot, blit=True)
        self.ani.save('aescatter.mp4', fps=60, dpi=300,
                      extra_args=['-vcodec', 'libx264'])

    def setup_plot(self):
        # First frame
        """Initial drawing of the scatter plot."""
        t, name, mass, radius, npl, pl, radmarker = next(self.data_stream(0))

        # set up the figure
        self.ax = plt.axes(xlim=(xmin, xmax), ylim=(ymin, ymax))
        self.ax.margins(x=10, y=1)
        self.ax.set_xlabel('Semi Major Axis ($R_{Mars}$)', fontsize='16', labelpad=1)
        self.ax.set_ylabel('Eccentricity', fontsize='16', labelpad=1)
        self.ax.set_yscale('log')

        self.title = self.ax.text(0.50, 1.05, "", bbox={'facecolor': 'w', 'alpha': 0.5, 'pad': 5}, transform=self.ax.transAxes,
                        ha="center")

        self.title.set_text(f'Time = ${t / 24 / 3600:4.1f}$ days with ${npl:f}$ particles')
        self.scat = self.ax.scatter(pl[:, 0], pl[:, 1], marker='o', s=radmarker, c='red', alpha=0.25)
        return self.scat, self.title

    def data_stream(self, frame=0):
        while True:
            d = self.ds.isel(time=frame)
            radius = d['radmarker'].values
            mass = d['mass'].values
            a = d['a'].values / RMars
            e = d['e'].values
            name = d['id'].values
            npl = d['npl'].values
            radmarker = d['radmarker']
            t = self.ds.coords['time'].values[frame]
            frame += 1
            yield t, name, mass, radius, npl, np.c_[a, e], radmarker

    def update(self,frame):
        """Update the scatter plot."""
        t, name, mass, radius, npl, pl, radmarker = next(self.data_stream(frame))

        self.title.set_text(f'Time = ${t / 24 / 3600:4.1f}$ days with ${npl:4.0f}$ particles')

        # We need to return the updated artist for FuncAnimation to draw..
        # Note that it expects a sequence of artists, thus the trailing comma.
        self.scat.set_sizes(radmarker)
        self.scat.set_offsets(pl)
        self.scat.set_facecolor('red')

        return self.scat, self.title,


config = swio.read_swiftest_config("param.in")
marsdisk = swio.swiftest2xr(config)
plt.show()


anim = AnimatedScatter(marsdisk,config)
