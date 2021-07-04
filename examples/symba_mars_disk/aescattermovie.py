import swiftest.io as swio
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import matplotlib.colors as mcolors

radscale = 20
RMars = 3389500.0
xmin = 1.0
xmax = 10.0
ymin = 1e-6
ymax = 1.0

class AnimatedScatter(object):
    """An animated scatter plot using matplotlib.animations.FuncAnimation."""
    def __init__(self, ds, param):
        #outf = 'aescatter.mp4'

        frame = 0
        nframes = ds['time'].size
        self.ds = ds
        self.param = param
        self.ds['Mass'] = self.ds['Mass'] / param['GU']
        self.ds['radmarker'] = self.ds['Radius'].fillna(0)
        self.ds['radmarker'] = self.ds['radmarker'] / self.ds['radmarker'].max() * radscale

        self.clist = {'Initial conditions' : 'xkcd:faded blue',
                      'Disruption' : 'xkcd:marigold',
                      'Supercatastrophic' : 'xkcd:shocking pink',
                      'Hit and run fragment' : 'xkcd:baby poop green'}

        self.stream = self.data_stream(frame)
        # Setup the figure and axes...
        self.fig, self.ax = plt.subplots(figsize=(8,4.5))
        # Then setup FuncAnimation.
        self.ani = animation.FuncAnimation(self.fig, self.update, interval=1, frames=nframes,
                                          init_func=self.setup_plot, blit=True)
        self.ani.save('aescatter.mp4', fps=60, dpi=300,
                      extra_args=['-vcodec', 'libx264'])
        print('Finished writing aescattter.mp4')

    def scatters(self, pl, radmarker, origin):
        scat = []
        for key, value in self.clist.items():
            idx = origin == value
            s = self.ax.scatter(pl[idx, 0], pl[idx, 1], marker='o', s=radmarker[idx], c=value, alpha=0.25, label=key)
            scat.append(s)
        return scat

    def setup_plot(self):
        # First frame
        """Initial drawing of the scatter plot."""
        t, name, Mass, Radius, npl, pl, radmarker, origin = next(self.data_stream(0))

        # set up the figure
        self.ax = plt.axes(xlim=(xmin, xmax), ylim=(ymin, ymax))
        self.ax.margins(x=10, y=1)
        self.ax.set_xlabel('Semi Major Axis ($R_{Mars}$)', fontsize='16', labelpad=1)
        self.ax.set_ylabel('Eccentricity', fontsize='16', labelpad=1)
        self.ax.set_yscale('log')

        self.title = self.ax.text(0.50, 1.05, "", bbox={'facecolor': 'w', 'alpha': 0.5, 'pad': 5}, transform=self.ax.transAxes,
                        ha="center")

        self.title.set_text(f'Time = ${t / 24 / 3600:4.1f}$ days with ${npl:f}$ particles')
        slist = self.scatters(pl, radmarker, origin)
        self.s0 = slist[0]
        self.s1 = slist[1]
        self.s2 = slist[2]
        self.s3 = slist[3]
        self.ax.legend(loc='upper right')
        return self.s0, self.s1, self.s2, self.s3, self.title

    def data_stream(self, frame=0):
        while True:
            d = self.ds.isel(time=frame)
            Radius = d['radmarker'].values
            Mass = d['Mass'].values
            a = d['a'].values / RMars
            e = d['e'].values
            name = d['id'].values
            npl = d['npl'].values
            radmarker = d['radmarker']
            origin = d['origin_type']

            t = self.ds.coords['time'].values[frame]

            frame += 1
            yield t, name, Mass, Radius, npl, np.c_[a, e], radmarker, origin

    def update(self,frame):
        """Update the scatter plot."""
        t, name, Mass, Radius, npl, pl, radmarker, origin = next(self.data_stream(frame))

        self.title.set_text(f'Time = ${t / 24 / 3600:4.1f}$ days with ${npl:4.0f}$ particles')

        # We need to return the updated artist for FuncAnimation to draw..
        # Note that it expects a sequence of artists, thus the trailing comma.
        s = [self.s0, self.s1, self.s2, self.s3]
        for i, (key, value) in enumerate(self.clist.items()):
            idx = origin == key
            s[i].set_sizes(radmarker[idx])
            s[i].set_offsets(pl[idx,:])
            s[i].set_facecolor(value)

        self.s0 = s[0]
        self.s1 = s[1]
        self.s2 = s[2]
        self.s3 = s[3]
        return self.s0, self.s1, self.s2, self.s3, self.title,


param = swio.read_swiftest_param("param.in")
marsdisk = swio.swiftest2xr(param)
print('Making animation')
anim = AnimatedScatter(marsdisk,param)
print('Animation finished')