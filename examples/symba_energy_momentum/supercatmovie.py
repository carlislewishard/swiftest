import swiftestio as swio
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import matplotlib.collections as clt

xmin = -8.0
xmax = 8.0
ymin = -8.0
ymax = 8.0

outfile = 'supercat.mp4'

def scale_sim(ds, config):

    dsscale = ds

    dsscale['mass'] = ds['mass'] / config['GU']
    Mtot = dsscale['mass'].sum(skipna=True, dim="id").isel(time=0)
    rscale = sum(ds['radius'].sel(id=[2, 3], time=0)).item()
    ds['radius'] /= rscale

    dsscale['radmarker'] = dsscale['radius'].fillna(0)

    dsscale['px'] /= rscale
    dsscale['py'] /= rscale
    dsscale['pz'] /= rscale

    mpx = dsscale['mass'] * dsscale['px']
    mpy = dsscale['mass'] * dsscale['py']
    mpz = dsscale['mass'] * dsscale['pz']
    xbsys = mpx.sum(skipna=True, dim="id") / Mtot
    ybsys = mpy.sum(skipna=True, dim="id") / Mtot
    zbsys = mpz.sum(skipna=True, dim="id") / Mtot

    mvx = dsscale['mass'] * dsscale['vx']
    mvy = dsscale['mass'] * dsscale['vy']
    mvz = dsscale['mass'] * dsscale['vz']
    vxbsys = mvx.sum(skipna=True, dim="id") / Mtot
    vybsys = mvy.sum(skipna=True, dim="id") / Mtot
    vzbsys = mvz.sum(skipna=True, dim="id") / Mtot

    dsscale['pxb'] = dsscale['px'] - xbsys
    dsscale['pyb'] = dsscale['py'] - ybsys
    dsscale['pzb'] = dsscale['pz'] - zbsys

    dsscale['vxb'] = dsscale['vx'] - vxbsys
    dsscale['vyb'] = dsscale['vy'] - vybsys
    dsscale['vzb'] = dsscale['vz'] - vzbsys

    return dsscale

class UpdatablePatchCollection(clt.PatchCollection):
    def __init__(self, patches, *args, **kwargs):
        self.patches = patches
        clt.PatchCollection.__init__(self, patches, *args, **kwargs)

    def get_paths(self):
        self.set_paths(self.patches)
        return self._paths

class AnimatedScatter(object):
    """An animated scatter plot using matplotlib.animations.FuncAnimation."""

    def __init__(self, ds, config):

        frame = 0
        nframes = ds['time'].size
        self.ds = scale_sim(ds, config)
        self.config = config

        self.clist = {'Initial conditions' : 'xkcd:faded blue',
                      'Disruption' : 'xkcd:marigold',
                      'Supercatastrophic' : 'xkcd:shocking pink',
                      'Hit and run fragment' : 'xkcd:baby poop green',
                      'Central body'    : 'black'}

        self.stream = self.data_stream(frame)
        # Setup the figure and axes...
        self.fig, self.ax = plt.subplots(figsize=(8,8))
        plt.tight_layout()
        # Then setup FuncAnimation.
        self.ani = animation.FuncAnimation(self.fig, self.update, interval=1, frames=nframes,
                                          init_func=self.setup_plot, blit=False)
        self.ani.save(outfile, fps=60, dpi=300,
                      extra_args=['-vcodec', 'libx264'])

    def plot_pl_circles(self, pl, radmarker):
        patches = []
        for i in range(pl.shape[0]):
            s = plt.Circle((pl[i, 0], pl[i, 1]), radmarker[i])
            patches.append(s)
        return patches

    def vec_props(self, c):
        arrowprops = {
            'arrowstyle': '<|-',
            'mutation_scale': 20,
            'connectionstyle': 'arc3',
        }

        arrow_args = {
            'xycoords': 'data',
            'textcoords': 'data',
            'arrowprops': arrowprops,
            'annotation_clip': True,
            'zorder': 100,
            'animated' : True
        }
        aarg = arrow_args.copy()
        aprop = arrowprops.copy()
        aprop['color'] = c
        aarg['arrowprops'] = aprop
        aarg['color'] = c
        return aarg

    def plot_pl_vectors(self, pl, cval):

        arrows = []
        for i in range(pl.shape[0]):
            aarg = self.vec_props(cval[i])
            varrowend = (pl[i, 0], pl[i, 1])
            varrowtip = (pl[i, 0] + pl[i, 2] * self.v_length, pl[i, 1] + pl[i, 3] * self.v_length)
            a = self.ax.annotate("",xy=varrowend,xytext=varrowtip, **aarg)
            arrows.append(a)
        return arrows

    def origin_to_color(self, origin):
        cval = []
        for o in origin:
           c = self.clist[o]
           cval.append(c)

        return cval

    def setup_plot(self):
        # First frame
        """Initial drawing of the scatter plot."""
        t, name, mass, radius, npl, pl, radmarker, origin = next(self.data_stream(0))

        cval = self.origin_to_color(origin)
        # set up the figure
        self.ax = plt.axes(xlim=(xmin, xmax), ylim=(ymin, ymax))
        self.ax.set_aspect(1)

        # Scale markers to the size of the system
        self.v_length = 0.50  # Length of arrow as fraction of velocity

        self.ax.margins(x=1, y=1)
        self.ax.set_xlabel('x distance / ($R_1 + R_2$)', fontsize='16', labelpad=1)
        self.ax.set_ylabel('y distance / ($R_1 + R_2$)', fontsize='16', labelpad=1)

        self.title = self.ax.text(0.50, 1.05, "", bbox={'facecolor': 'w', 'alpha': 0.5, 'pad': 5}, transform=self.ax.transAxes,
                        ha="center")

        self.title.set_text('Supercatastrophic')
        self.patches = self.plot_pl_circles(pl, radmarker)

        self.collection = UpdatablePatchCollection(self.patches, color=cval)
        self.ax.add_collection(self.collection)
        self.arrows = self.plot_pl_vectors(pl, cval)

        return self.collection, self.arrows

    def update(self,frame):
        """Update the scatter plot."""
        t, name, mass, radius, npl, pl, radmarker, origin = next(self.data_stream(frame))
        cval = self.origin_to_color(origin)
        for i, p in enumerate(self.patches):
            p.set_center((pl[i, 0], pl[i,1]))
            p.set_radius(radmarker[i])
            p.set_color(cval[i])
            varrowend = (pl[i, 0], pl[i, 1])
            varrowtip = (pl[i, 0] + pl[i, 2] * self.v_length, pl[i, 1] + pl[i, 3] * self.v_length)
            self.arrows[i].set_position(varrowtip)
            self.arrows[i].xy = varrowend

        self.collection.set_paths(self.patches)
        return self.collection, self.arrows

    def data_stream(self, frame=0):
        while True:
            d = self.ds.isel(time=frame)
            radius = d['radmarker'].values
            mass = d['mass'].values
            x = d['pxb'].values
            y = d['pyb'].values
            vx = d['vxb'].values
            vy = d['vyb'].values
            name = d['id'].values
            npl = d['npl'].values
            radmarker = d['radmarker'].values
            origin = d['origin_type'].values

            t = self.ds.coords['time'].values[frame]

            frame += 1
            yield t, name, mass, radius, npl, np.c_[x, y, vx, vy], radmarker, origin

config = swio.read_swiftest_config("param.supercatastrophic.in")
ds = swio.swiftest2xr(config)
print('Making animation')
anim = AnimatedScatter(ds,config)
print('Animation finished')
