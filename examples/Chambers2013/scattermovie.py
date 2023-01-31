#!/usr/bin/env python3
import swiftest 
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import matplotlib.colors as mcolors
from collections import namedtuple
plt.switch_backend('agg')

titletext = "Chambers (2013)"
valid_plot_styles = ["aescatter", "aiscatter"]
xlim={"aescatter" : (0.0, 2.5),
      "aiscatter" : (0.0, 2.5)}
ylim={"aescatter" : (0.0, 1.0),
      "aiscatter" : (0.0, 40.0)}
xlabel={"aescatter": "Semimajor axis (AU)",
        "aiscatter": "Semimajor axis (AU)"}
ylabel={"aescatter": "Eccentricity",
        "aiscatter": "Inclination (deg)"}


plot_style = valid_plot_styles[0]
framejump = 1
animation_file = f"Chambers2013-{plot_style}.mp4"

origin_types = ["Initial conditions", "Merger", "Disruption", "Supercatastrophic", "Hit and run fragmentation"]

class AnimatedScatter(object):
    """An animated scatter plot using matplotlib.animations.FuncAnimation."""
    def __init__(self, ds, param):

        self.radscale = 2000
        nframes = int(ds['time'].size / framejump)
        self.ds = ds
        self.param = param
        self.Rcb = self.ds['radius'].sel(name="Sun").isel(time=0).values[()]
        colors = ["k", "xkcd:faded blue", "xkcd:marigold", "xkcd:shocking pink", "xkcd:baby poop green"]
        self.clist = dict(zip(origin_types,colors))


        # Setup the figure and axes...
        fig = plt.figure(figsize=(8,4.5), dpi=300)
        plt.tight_layout(pad=0)
        # set up the figure
        self.ax = plt.Axes(fig, [0.1, 0.15, 0.8, 0.75])
        fig.add_axes(self.ax)

        self.make_artists()
        
        self.ani = animation.FuncAnimation(fig, func=self.update, interval=1, frames=nframes, init_func=self.init_plot, blit=True)
        self.ani.save(animation_file, fps=60, dpi=300, extra_args=['-vcodec', 'libx264'])
        print(f'Finished writing {animation_file}')
        
    def make_artists(self):
        scatter_names = [f"s{i}" for i,k in enumerate(origin_types)]
        self.scatter_artist_names = dict(zip(origin_types,scatter_names))
        
        animated_elements = [self.ax.text(0.50, 1.05, "", bbox={'facecolor': 'w', 'alpha': 0.5, 'pad': 5}, transform=self.ax.transAxes,ha="center", animated=True)]
        element_names = ["title"]
        for key, value in self.clist.items():
            animated_elements.append(self.ax.scatter([], [], marker='o', s=[], c=value, alpha=0.75, label=key, animated=True))
            element_names.append(self.scatter_artist_names[key])
        
        Artists = namedtuple("Artists",tuple(element_names))
        self.artists = Artists(*animated_elements)
        return 

    def init_plot(self):
        self.ax.set_xlim(xlim[plot_style])
        self.ax.set_ylim(ylim[plot_style])
        
        # set up the figure
        self.ax.margins(x=10, y=1)
        self.ax.set_xlabel(xlabel[plot_style], fontsize='16', labelpad=1)
        self.ax.set_ylabel(ylabel[plot_style], fontsize='16', labelpad=1) 
        
        leg = plt.legend(loc="upper left", scatterpoints=1, fontsize=10)
        for i,l in enumerate(leg.legendHandles):
           leg.legendHandles[i]._sizes = [20]      
        
        return self.artists

    def get_data(self, frame=0):
        d = self.ds.isel(time = frame)
        n=len(d['name'])
        d = d.isel(name=range(1,n))
        d['radmarker'] = (d['radius'] / self.Rcb) * self.radscale

        t = d['time'].values
        npl = d['npl'].values
        radmarker = d['radmarker'].values
        origin = d['origin_type'].values
        
        if plot_style == "aescatter":
            pl = np.c_[d['a'].values,d['e'].values]
        elif plot_style == "aiscatter":
            pl = np.c_[d['a'].values,d['inc'].values]

        return t, npl, pl, radmarker, origin

    def update(self,frame):
        """Update the scatter plot."""
        t,  npl, pl, radmarker, origin = self.get_data(framejump * frame)

        self.artists.title.set_text(f"{titletext} - Time = ${t*1e-6:6.3f}$ My with ${npl:4.0f}$ particles")

        for key,name in self.scatter_artist_names.items():
            idx = origin == key
            if any(idx) and any(~np.isnan(radmarker[idx])):
                scatter = self.artists._asdict()[name]
                scatter.set_sizes(radmarker[idx])
                scatter.set_offsets(pl[idx,:])
                scatter.set_facecolor(self.clist[key])

        return self.artists

sim = swiftest.Simulation(read_data=True)
print('Making animation')
anim = AnimatedScatter(sim.data,sim.param)
print('Animation finished')
