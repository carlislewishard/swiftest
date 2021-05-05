def makemovie(ds):
    # set up the figure
    fig = plt.figure(figsize=(16, 9))
    canvas_width, canvas_height = fig.canvas.get_width_height()

    def data_stream(frame):
        while True:
            d = ds.isel(time=frame)
            ds['Mass'] = ds['Mass'] / config['GU']
            ds['radmarker'] = ds['Radius'].fillna(0)
            ds['radmarker'] = ds['radmarker'] / ds['radmarker'].max()
            ds['density'] = 3 * ds['Mass'] / (4 * np.pi * (ds['Radius']) ** 3)

            mass = d['Mass'].values
            radius = d['Radius'].values / radscale
            density = d['density'].values
            a = d['a'].values / RMars
            e = d['e'].values
            name = d['id'].values
            npl = d['npl'].values
            radmarker = d['radmarker'] * radscale
            t = ds.coords['time'].values[frame]
            yield t, name, mass, radius, npl, np.c_[a, e], density, radmarker

    def update(frame):
        """Update the scatter plot."""
        t, name, mass, radius, npl, pl, density, radmarker = next(data_stream(frame))

        title.set_text(f'Time = ${t / 24 / 3600:4.1f}$ days with ${npl:4.0f}$ particles')

        # We need to return the updated artist for FuncAnimation to draw..
        # Note that it expects a sequence of artists, thus the trailing comma.
        scat.set_sizes(radmarker)
        scat.set_offsets(pl)
        scat.set_facecolor('red')

        return scat, title,

    # First frame
    """Initial drawing of the scatter plot."""
    t, name, mass, radius, npl, pl, density, radmarker = next(data_stream(0))

    ax = plt.axes(xlim=(xmin, xmax), ylim=(ymin, ymax))
    ax.margins(x=10, y=1)
    ax.set_xlabel('Semi Major Axis ($R_{Mars}$)', fontsize='16', labelpad=1)
    ax.set_ylabel('Eccentricity', fontsize='16', labelpad=1)
    ax.set_yscale('log')

    # ax.ticklabel_format(axis='both', style='', scilimits=None, useOffset=1, useLocale=None, useMathText=True)

    title = ax.text(0.50, 1.05, "", bbox={'facecolor': 'w', 'alpha': 0.5, 'pad': 5}, transform=ax.transAxes,
                    ha="center")

    # reference lines
    ax.scatter([a_Deimos, a_Phobos], [e_Deimos, e_Phobos], marker='o', s=np.power([r_Deimos, r_Phobos], 2) * 0.0000005,
               c=satcolor, alpha=0.5)
    #RRLline = ax.plot([RRL / RP, RRL / RP], [ymin, ymax], '--', color=RRLcolor, linewidth=1.5, zorder=50)
    #RRLlab = ax.text(RRL / RP - 0.00, 0.9 * ymax, "RRL", color=RRLcolor, rotation=0, fontsize=tsize, ha='center')
    #FRLline = ax.plot([FRL / RP, FRL / RP], [ymin, ymax], ':', color=FRLcolor, linewidth=1.5, zorder=50)
    #FRLlab = ax.text(FRL / RP - 0.00, 0.9 * ymax, "FRL", color=FRLcolor, rotation=0, fontsize=tsize, ha='center')
    #Rsyncline = ax.plot([Rsync / RP, Rsync / RP], [ymin, ymax], '-.', color=asynccolor, linewidth=1.5, zorder=50)
    #Rsynclab = ax.text(Rsync / RP + 0.10, 0.9 * ymax, "$a_{sync}$", color=asynccolor, rotation=0, fontsize=tsize,
    #                   ha='center')
    #Rlindline = ax.plot([alind / RP, alind / RP], [ymin, ymax], '-.', color=alindcolor, linewidth=1.5, zorder=50)
    #Rlindlab = ax.text(alind / RP - 0.02, 0.9 * ymax, "$a_{lind}$", color=alindcolor, rotation=0, fontsize=tsize,
    #                   ha='center')

    title.set_text(f'Time = ${t / 24 / 3600:4.1f}$ days with ${npl:f}$ particles')
    scat = ax.scatter(pl[:, 0], pl[:, 1], marker='o', s=radmarker, c='red', alpha=0.25)  # , vmin = 1000, vmax = 3000)

    # Open an ffmpeg process
    outf = 'aescatter.mp4'
    cmdstring = ('ffmpeg',
                 '-y', '-r', '1',  # overwrite, 1fps
                 '-s', '%dx%d' % (canvas_width, canvas_height),  # size of image string
                 '-pix_fmt', 'argb',  # format
                 '-f', 'rawvideo', '-i', '-',  # tell ffmpeg to expect raw video from the pipe
                 '-vcodec', 'mpeg4', outf)  # output encoding
    p = subprocess.Popen(cmdstring, stdin=subprocess.PIPE)

    # Draw frames and write to the pipe
    for frame in range(nframes):
        print(f'Generating frame {frame} of {nframes}')
        # draw the frame
        scat, title, = update(frame)
        fig.canvas.draw()

        # extract the image as an ARGB string
        string = fig.canvas.tostring_argb()

        # write to pipe
        p.stdin.write(string)

    # Finish up
    p.communicate()
