#!/usr/bin/env python3
"""
 Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
 This file is part of Swiftest.
 Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
 as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
 Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
 of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 You should have received a copy of the GNU General Public License along with Swiftest. 
 If not, see: https://www.gnu.org/licenses. 
"""

"""
Generates a movie of a fragmentation event from set of Swiftest output files.

Inputs
_______
param.in : ASCII text file
    Swiftest parameter input file.
out.nc   : NetCDF file
    Swiftest output file.

Returns
-------
fragmentation.mp4 : mp4 movie file
    Movie of a fragmentation event.
"""

import swiftest
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# ----------------------------------------------------------------------------------------------------------------------
# Define the names and initial conditions of the various fragmentation simulation types
# ----------------------------------------------------------------------------------------------------------------------
available_movie_styles = ["disruption_headon", "disruption_off_axis", "supercatastrophic_headon", "supercatastrophic_off_axis","hitandrun_disrupt", "hitandrun_pure", "merge"]
movie_title_list = ["Head-on Disruption", "Off-axis Disruption", "Head-on Supercatastrophic", "Off-axis Supercatastrophic", "Hit and Run w/ Runner Disruption", "Pure Hit and Run", "Merge"]
movie_titles = dict(zip(available_movie_styles, movie_title_list))
num_movie_frames = 1200

# These initial conditions were generated by trial and error
names = ["Target","Projectile"]
pos_vectors = {"disruption_headon"         : [np.array([1.0, -5.0e-05, 0.0]),
                                              np.array([1.0,  5.0e-05 ,0.0])],
              "disruption_off_axis"        : [np.array([1.0, -5.0e-05, 0.0]),
                                              np.array([1.0,  5.0e-05 ,0.0])], 
               "supercatastrophic_headon":   [np.array([1.0, -5.0e-05, 0.0]),
                                              np.array([1.0,  5.0e-05, 0.0])],
               "supercatastrophic_off_axis": [np.array([1.0, -5.0e-05, 0.0]),
                                              np.array([1.0,  5.0e-05, 0.0])],
               "hitandrun_disrupt"         : [np.array([1.0, -4.2e-05, 0.0]),
                                              np.array([1.0,  4.2e-05, 0.0])],
               "hitandrun_pure"            : [np.array([1.0, -4.2e-05, 0.0]),
                                              np.array([1.0,  4.2e-05, 0.0])],
               "merge"                      : [np.array([1.0, -5.0e-05, 0.0]),
                                              np.array([1.0,  5.0e-05 ,0.0])]                
               }

vel_vectors = {"disruption_headon"         : [np.array([ 0.00,  6.280005, 0.0]),
                                              np.array([ 0.00, -6.280005, 0.0])],
               "disruption_off_axis"       : [np.array([ 0.00,  6.280005, 0.0]),
                                              np.array([ 0.50, -6.280005, 0.0])], 
               "supercatastrophic_headon":   [np.array([ 0.00,  6.28,     0.0]),
                                              np.array([ 0.00, -6.28,     0.0])],
               "supercatastrophic_off_axis": [np.array([ 0.00,  6.28,     0.0]),
                                              np.array([ 0.50, -6.28,     0.0])],
               "hitandrun_disrupt"         : [np.array([ 0.00,  6.28,     0.0]),
                                              np.array([-1.45, -6.28,     0.0])],
               "hitandrun_pure"            : [np.array([ 0.00,  6.28,     0.0]),
                                              np.array([-1.52, -6.28,     0.0])],
               "merge"                     : [np.array([ 0.00,  0.0, 0.0]),
                                              np.array([ 0.01, -0.100005, 0.0])] 
               }

rot_vectors = {"disruption_headon"         : [np.array([0.0, 0.0, 0.0]),
                                              np.array([0.0, 0.0, 0.0])],
               "disruption_off_axis":        [np.array([0.0, 0.0, -6.0e3]),
                                              np.array([0.0, 0.0, 1.0e4])],
               "supercatastrophic_headon":   [np.array([0.0, 0.0, 0.0]),
                                              np.array([0.0, 0.0, 0.0])],
               "supercatastrophic_off_axis": [np.array([0.0, 0.0, -6.0e3]),
                                              np.array([0.0, 0.0, 1.0e4])],
               "hitandrun_disrupt"         : [np.array([0.0, 0.0, 6.0e3]),
                                              np.array([0.0, 0.0, 1.0e4])],
               "hitandrun_pure"            : [np.array([0.0, 0.0, 6.0e3]),
                                              np.array([0.0, 0.0, 1.0e4])],
               "merge"                     : [np.array([0.0, 0.0, -6.0e3]),
                                              np.array([0.0, 0.0, 1.0e4])] 
               }

body_Gmass = {"disruption_headon"        : [1e-7, 1e-10],
             "disruption_off_axis"       : [1e-7, 1e-10],
             "supercatastrophic_headon"  : [1e-7, 1e-8],
             "supercatastrophic_off_axis": [1e-7, 1e-8],
             "hitandrun_disrupt"         : [1e-7, 7e-10],
             "hitandrun_pure"            : [1e-7, 7e-10],
             "merge"                     : [1e-7, 1e-8] 
               }

tstop = {"disruption_headon"         : 5.0e-4,
         "disruption_off_axis"       : 5.0e-4,
         "supercatastrophic_headon"  : 5.0e-4,
         "supercatastrophic_off_axis": 5.0e-4,
         "hitandrun_disrupt"         : 2.0e-4,
         "hitandrun_pure"            : 2.0e-4,
         "merge"                     : 2.0e-3,
         }

density = 3000 * swiftest.AU2M**3 / swiftest.MSun
GU = swiftest.GMSun * swiftest.YR2S**2 / swiftest.AU2M**3
body_radius = body_Gmass.copy()
for k,v in body_Gmass.items():
    body_radius[k] = [((Gmass/GU)/(4./3.*np.pi*density))**(1./3.) for Gmass in v]

body_radius["hitandrun_disrupt"] = [7e-6, 3.25e-6] 
body_radius["hitandrun_pure"] = [7e-6, 3.25e-6] 

# ----------------------------------------------------------------------------------------------------------------------
# Define the animation class that will generate the movies of the fragmentation outcomes
# ----------------------------------------------------------------------------------------------------------------------


def encounter_combiner(sim):
    """
    Combines simulation data with encounter data to produce a dataset that contains the position,
    mass, radius, etc. of both. It will interpolate over empty time values to fill in gaps.
    """

    # Only keep a minimal subset of necessary data from the simulation and encounter datasets
    keep_vars = ['rh','vh','Gmass','radius']
    data = sim.data[keep_vars]
    enc = sim.encounters[keep_vars].load()

    # Remove any encounter data at the same time steps that appear in the data to prevent duplicates
    t_not_duplicate = ~enc['time'].isin(data['time'])
    enc = enc.where(t_not_duplicate,drop=True)
    tgood=enc.time.where(~np.isnan(enc.time),drop=True)
    enc = enc.sel(time=tgood)

    # The following will combine the two datasets along the time dimension, sort the time dimension, and then fill in any time gaps with interpolation
    ds = xr.combine_nested([data,enc],concat_dim='time').sortby("time").interpolate_na(dim="time")
   
    # Interpolate in time to make a smooth, constant time step dataset 
    smooth_time = np.linspace(start=tgood.isel(time=0), stop=ds.time[-1], num=num_movie_frames)
    ds = ds.interp(time=smooth_time)

    return ds

class AnimatedScatter(object):
    """An animated scatter plot using matplotlib.animations.FuncAnimation."""

    def __init__(self, sim, animfile, title, style, nskip=1):

        self.ds = encounter_combiner(sim)
        nframes = int(self.ds['time'].size)
        self.sim = sim
        self.title = title
        self.body_color_list = {'Initial conditions': 'xkcd:windows blue',
                      'Disruption': 'xkcd:baby poop',
                      'Supercatastrophic': 'xkcd:shocking pink',
                      'Hit and run fragmention': 'xkcd:blue with a hint of purple',
                      'Central body': 'xkcd:almost black'}

        # Set up the figure and axes...
        self.figsize = (4,4)
        self.fig, self.ax = self.setup_plot()

        # Then setup FuncAnimation.
        self.ani = animation.FuncAnimation(self.fig, self.update_plot, interval=1, frames=range(0,nframes,nskip), blit=True)
        self.ani.save(animfile, fps=60, dpi=300, extra_args=['-vcodec', 'libx264'])
        print(f"Finished writing {animfile}")

    def setup_plot(self):
        fig = plt.figure(figsize=self.figsize, dpi=300)
        plt.tight_layout(pad=0)

        # Calculate the distance along the y-axis between the colliding bodies at the start of the simulation.
        # This will be used to scale the axis limits on the movie.
        rhy1 = self.ds['rh'].sel(name="Target",space='y').isel(time=0).values[()]
        rhy2 = self.ds['rh'].sel(name="Projectile",space='y').isel(time=0).values[()]

        scale_frame =   abs(rhy1) + abs(rhy2)
        if "hitandrun" in style:
           scale_frame *= 2
           
        ax = plt.Axes(fig, [0.1, 0.1, 0.8, 0.8])
        self.ax_pt_size = self.figsize[0] *  72 / scale_frame
        ax.set_xlim(-scale_frame, scale_frame)
        ax.set_ylim(-scale_frame, scale_frame)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_title(self.title)
        fig.add_axes(ax)

        self.scatter_artist = ax.scatter([], [], animated=True, c='k', edgecolors='face')
        return fig, ax

    def update_plot(self, frame):
        
        # Define a function to calculate a reference frame for the animation
        # This will be based on the initial velocity of the Target body
        def reference_frame(r_ref, v_ref, t):
            coord_pos = r_ref + v_ref * t 
            return coord_pos.values[0], coord_pos.values[1]

        t, Gmass, rh, point_rad = next(self.data_stream(frame))
        x_ref, y_ref = reference_frame(self.r_ref,self.v_ref, t)
        self.scatter_artist.set_offsets(np.c_[rh[:,0] - x_ref, rh[:,1] - y_ref])
        self.scatter_artist.set_sizes(point_rad**2)
        return self.scatter_artist,

    def data_stream(self, frame=0):
        while True:
            ds = self.ds.isel(time=frame)
            ds = ds.where(ds['name'] != "Sun", drop=True)
            t = ds['time'].values[()]
            radius = ds['radius'].values
            Gmass = ds['Gmass'].values
            rh = ds['rh'].values
            point_rad = radius * self.ax_pt_size
            
            # Save the initial velocity of body 1 to use as a reference
            if frame == 0:
                self.r_ref = ds.sel(name="Target")['rh'] 
                self.v_ref = ds.sel(name="Target")['vh'] 
            yield t, Gmass, rh, point_rad

if __name__ == "__main__":

    print("Select a fragmentation movie to generate.")
    print("1. Head-on disruption")
    print("2. Off-axis disruption")
    print("3. Head-on supercatastrophic")
    print("4. Off-axis supercatastrophic")
    print("5. Hit and run with disruption of the runner")
    print("6. Pure hit and run")
    print("7. Merge")
    print("8. All of the above")
    user_selection = int(input("? "))

    if user_selection > 0 and user_selection < 8:
        movie_styles = [available_movie_styles[user_selection-1]]
    else:
        print("Generating all movie styles")
        movie_styles = available_movie_styles.copy()

    for style in movie_styles:
        print(f"Generating {movie_titles[style]}")
        movie_filename = f"{style}.mp4"
        # Pull in the Swiftest output data from the parameter file and store it as a Xarray dataset.
        sim = swiftest.Simulation(simdir=style, rotation=True, init_cond_format = "XV", compute_conservation_values=True)
        sim.add_solar_system_body("Sun")
        sim.add_body(name=names, Gmass=body_Gmass[style], radius=body_radius[style], rh=pos_vectors[style], vh=vel_vectors[style], rot=rot_vectors[style])

        # Set fragmentation parameters
        minimum_fragment_gmass = 0.05 * body_Gmass[style][1] 
        gmtiny = 0.10 * body_Gmass[style][1] 
        sim.set_parameter(collision_model="fraggle", encounter_save="both", gmtiny=gmtiny, minimum_fragment_gmass=minimum_fragment_gmass, verbose=False)
        sim.run(dt=5e-4, tstop=tstop[style], istep_out=1, dump_cadence=0)

        print("Generating animation")
        anim = AnimatedScatter(sim,movie_filename,movie_titles[style],style,nskip=1)
