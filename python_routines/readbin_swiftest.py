import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
import os
import sys
from scipy.io import FortranFile
#CALL python3 readbin_swiftest.py [output_path] [tstop] [dt] [mars/sun]

figure = plt.figure(1, figsize=(8,6))
filepath = sys.argv[1]
tstop = np.float(sys.argv[2])
dt = np.float(sys.argv[3])
runtype = sys.argv[4]
RMars = 3.394e6 #Rmars in cm

if runtype == 'mars':
    xmin = 2.0
    xmax = 8.0
    ymin = 0.0
    ymax = 0.02
else:
    xmin = 0.0
    xmax = 5.0
    ymin = 0.0
    ymax = 1.0

class AnimatedScatter(object):
    """An animated scatter plot using matplotlib.animations.FuncAnimation."""
    def __init__(self):
        self.stream = self.data_stream()
        self.binfilename = filepath+"/"+'bin.dat'

        # Setup the figure and axes...
        self.fig, self.ax = plt.subplots()
        # Then setup FuncAnimation.
        self.frame_skip = 1000
        nframes = int(tstop/dt/self.frame_skip)
        print(f'Num frames: {nframes}')
        self.ani = animation.FuncAnimation(self.fig, self.update, interval=1, frames=nframes,
                                          init_func=self.setup_plot, blit=True)

        #self.ani.save('frames/charnoz2010-saturn-ringmoons.png', writer = "imagemagick")
        if runtype == "mars":
            self.ani.save(filepath+"/"+'martiansystem-ae.mp4', fps=60, dpi=600, extra_args=['-vcodec', 'libx264'])
        else:
            self.ani.save(filepath + "/" + 'solarsystem-ae.mp4', fps=60, dpi=600, extra_args=['-vcodec', 'libx264'])

    def setup_plot(self):
        """Initial drawing of the scatter plot."""
        t, name, mass, radius, npl, pl = next(self.stream)


        self.ax = plt.axes(xlim=(xmin, xmax), ylim=(ymin, ymax))
        if runtype == "mars":
            self.ax.set_xlabel('Semi Major Axis ($Rmars$)', fontsize='12')
        else:
            self.ax.set_xlabel('Semi Major Axis (AU)', fontsize='12')

        self.ax.set_ylabel('Eccentricity', fontsize='12')

        self.title = self.ax.text(0.80, 0.1, "", bbox={'facecolor': 'w', 'alpha': 0.5, 'pad': 5},
                        transform=self.ax.transAxes, ha="center")
        if runtype == "mars":
            self.title.set_text(f'Time = ${t[0]/24/3600:4.2f}$ days with ${npl[0]:f}$ particles')
            self.line = self.ax.scatter(pl[:, 0], pl[:, 1], marker='o', color="red", s=mass * 1e-4, zorder=50)
        else:
            self.title.set_text(f'Time = ${t[0]:4.2f}$ y with ${npl[0]:4.0f}$ particles')
            self.line = self.ax.scatter(pl[:, 0], pl[:, 1], marker='o', color="red", s=mass * 1e5, zorder=50)

        # For FuncAnimation's sake, we need to return the artist we'll be using
        # Note that it expects a sequence of artists, thus the trailing comma.
        return self.line, self.title,

    def data_stream(self):
        with FortranFile(self.binfilename, 'r') as f:
            while True:  # Loop until you read the end of file
                try:
                    t = f.read_reals(np.float64)  # Try first part of the header
                    print(t[0])
                except:
                    break
                if t[0] > tstop: break
                npl = f.read_ints()
                ntp = f.read_ints()
                iout_form = f.read_ints()
                name = f.read_ints()
                px = f.read_reals(np.float64)
                py = f.read_reals(np.float64)
                pz = f.read_reals(np.float64)
                vx = f.read_reals(np.float64)
                vy = f.read_reals(np.float64)
                vz = f.read_reals(np.float64)
                mass = f.read_reals(np.float64)
                radius = f.read_reals(np.float64)
                a = f.read_reals(np.float64)
                e = f.read_reals(np.float64)
                inc = f.read_reals(np.float64)
                yield t, name, mass, radius, npl, np.c_[a,e]


    def update(self, i):
        """Update the scatter plot."""
        t,name, mass, radius, npl, pl = next(self.stream)
        # Set x and y data...
        self.line.set_offsets(pl)

        # Set sizes...

        if runtype == "mars":
            self.title.set_text(f'Time = ${t[0] / 24 / 3600:4.2f}$ days with ${npl[0]:4.0f}$ particles')
            pl[:,0] = np.divide(pl[:,0],RMars)
        else:
            self.title.set_text(f'Time = ${t[0]:4.0f}$ y with ${npl[0]:4.0f}$ particles')
        # We need to return the updated artist for FuncAnimation to draw..
        # Note that it expects a sequence of artists, thus the trailing comma.
        return self.line, self.title,


if __name__ == '__main__':
    anim = AnimatedScatter()

    plt.show()

