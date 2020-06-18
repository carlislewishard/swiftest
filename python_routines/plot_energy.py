###### Run file to read and plot energy.out file
###### Created 05/21/2020 by Jennifer Pouplin
# CALL python3 plot_energy.py [mars/sun] [output_path] [filename]
################################################  IMPORTS ################################################
import matplotlib.pyplot as plt
import numpy as np
import sys
import os
import pandas as pd
#################################################  MAIN ###################################################



def main():
    filename = sys.argv[3]
    raw =[]
    filepath = sys.argv[2] + '/' + filename
    if not os.path.isfile(filepath):
        print("File path {} does not exist. Exiting...".format(filepath))
        sys.exit()

    df = pd.read_csv(filepath, delim_whitespace=True, names=('t', 'ke', 'pe', 'te', 'htotx','htoty','htotz' ),
            dtype={'t': np.float64, 'ke': np.float64, 'pe': np.float64,'te': np.float64,'htotx': np.float64,'htoty': np.float64,'htotz': np.float64})
    if sys.argv[1] == 'mars' :
        df.t = df.t /3600 /24 /365
    else:
        df.t = df.t
    df.pe =  (100*(df.pe - df.pe[0])/df.pe[0])
    df.ke = (100*(df.ke - df.ke[0]) / df.ke[0])
    df.te = abs(100*(df.te - df.te[0]) / df.te[0])

    new_array = df.to_numpy(dtype=float)
    time = new_array[[],1]
    time = df.t
    ahtotx0 = df.htotx[0]
    print(ahtotx0)
    ahtotx0 = ahtotx0**2

    ahtoty0 = df.htoty[0]
    ahtoty0 = ahtoty0**2
    ahtotz0 = df.htotz[0]
    ahtotz0 = ahtotz0**2
    ahtot0= np.sqrt(ahtotx0+ahtoty0+ahtotz0)

    ahtotx = np.power(df.htotx,2)
    ahtoty = np.power(df.htoty,2)
    ahtotz = np.power(df.htotz,2)
    ahtot = np.power(ahtotx+ahtoty+ahtotz, 0.5)
    Lorig = ahtot0
    Lnow = ahtot
    Ltot = abs(100*(Lnow - Lorig) / Lorig)


    df.htotx = abs(100*(df.htotx - df.htotx[0]) / df.htotx[0])
    df.htoty = abs(100*(df.htoty - df.htoty[0]) / df.htoty[0])
    df.htotz = abs(100*(df.htotz - df.htotz[0]) / df.htotz[0])

    fig, ax = plt.subplots()

    df.plot(x='t', y='te', color="slategrey", ax=ax)
    plt.title('Relative Change in Total Energy over time')
    ax.set_ylabel("% change")
    ax.set_xlabel("Time in years")
    #plt.show()
    plt.savefig(sys.argv[2] + '/'+ "plot_energyte.png")

    fig, ax = plt.subplots()


    df.plot(x='t', y='htotx', color="sandybrown", ax=ax)
    df.plot(x='t', y='htoty', color="slateblue", ax=ax)
    df.plot(x='t', y='htotz', color="#FF3030", ax=ax)
    plt.title('Relative Change in Angular momentum over time')
    # ax.set_ylim(0,3)
    ax.set_ylabel("% change")
    ax.set_xlabel("Time in years")
    #plt.show()
    plt.savefig(sys.argv[2] + '/' + "plot_energyhtot.png")

    fig, ax = plt.subplots()

    df.plot(x='t', y='ke', color="cadetblue", ax=ax)
    df.plot(x='t', y='pe', color="maroon", ax=ax)
    plt.title('Relative Change in Kinetic and Potential Energy over time')
    ax.set_ylabel("% change")
    ax.set_xlabel("Time in years")
    #plt.show()
    plt.savefig(sys.argv[2] + '/' + "plot_energykpe.png")

    fig, ax = plt.subplots()

    #print(Ltot)
    plt.plot(time, Ltot)
    plt.title('Relative Change in Angular momentum over time')
    # ax.set_ylim(0,3)
    ax.set_ylabel("% change")
    ax.set_xlabel("Time in years")
    #plt.show()
    plt.savefig(sys.argv[2] + '/' + "plot_energydl.png")


if __name__ == '__main__':
    main()