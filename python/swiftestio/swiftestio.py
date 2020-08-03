import numpy as np
import pandas as pd
from scipy.io import FortranFile
import xarray as xr
import dask
import os
import tempfile



#I/O Routines for reading in Swifter and Swiftest parameter and binary data files


GC = 6.6743E-11
einstinC = 299792458.0

def read_swifter_param(inparfile):
    """
    Reads in a Swifter param.in file and saves it as a dictionary

    Parameters
    ----------
    inparfile : string
        File name of the input parameter file

    Returns
    -------
    param
        A dictionary containing the entries in the user parameter file
    """
    param = {
    'INPARFILE'      : inparfile,
    'NPLMAX'         : -1,
    'NTPMAX'         : -1,
    'T0'             : 0.0,
    'TSTOP'          : 0.0,
    'DT'             : 0.0,
    'PL_IN'          : "",
    'TP_IN'          : "",
    'IN_TYPE'        : "ASCII",
    'ISTEP_OUT'      : -1,
    'BIN_OUT'        : "",
    'OUT_TYPE'       : 'REAL8',
    'OUT_FORM'       : "XV",
    'OUT_STAT'       : "NEW",
    'ISTEP_DUMP'     : -1,
    'J2'             : 0.0,
    'J4'             : 0.0,
    'CHK_CLOSE'      : 'NO',
    'CHK_RMIN'       : -1.0,
    'CHK_RMAX'       : -1.0,
    'CHK_EJECT'      : -1.0,
    'CHK_QMIN'       : -1.0,
    'CHK_QMIN_COORD' : "HELIO",
    'CHK_QMIN_RANGE' : "",
    'QMIN_ALO'       : -1.0,
    'QMIN_AHI'       : -1.0,
    'ENC_OUT'        : "",
    'EXTRA_FORCE'    : 'NO',
    'BIG_DISCARD'    : 'NO',
    'RHILL_PRESENT'  : 'NO',
    'GR'             : 'NO',
    'C2'             : -1.0,
             }

    # Read param.in file
    print(f'Reading Swifter file {inparfile}')
    f = open(inparfile, 'r')
    swifterlines = f.readlines()
    f.close()
    for line in swifterlines:
        fields = line.split()
        if len(fields) > 0:
            for key in param:
                if (key == fields[0].upper()): param[key] = fields[1]
            #Special case of CHK_QMIN_RANGE requires a second input
            if (param['CHK_QMIN_RANGE'] == fields[0].upper()):
                param['QMIN_ALO'] = fields[1]
                param['QMIN_AHI'] = fields[2]

    param['NPLMAX']     = int(param['NPLMAX'])
    param['NTPMAX']     = int(param['NTPMAX'])
    param['ISTEP_OUT']  = int(param['ISTEP_OUT'])
    param['ISTEP_DUMP'] = int(param['ISTEP_DUMP'])
    param['T0']         = float(param['T0'])
    param['TSTOP']      = float(param['TSTOP'])
    param['DT']         = float(param['DT'])
    param['J2']         = float(param['J2'])
    param['J4']         = float(param['J4'])
    param['CHK_RMIN']   = float(param['CHK_RMIN'])
    param['CHK_RMAX']   = float(param['CHK_RMAX'])
    param['CHK_EJECT']  = float(param['CHK_EJECT'])
    param['CHK_QMIN']   = float(param['CHK_QMIN'])
    param['QMIN_ALO']   = float(param['QMIN_ALO'])
    param['QMIN_AHI']   = float(param['QMIN_AHI'])
    param['EXTRA_FORCE'] = param['EXTRA_FORCE'].upper()
    param['BIG_DISCARD'] = param['BIG_DISCARD'].upper()
    param['CHK_CLOSE']  = param['CHK_CLOSE'].upper()
    param['RHILL_PRESENT'] = param['RHILL_PRESENT'].upper()

    return param

def read_swiftest_config(config_file_name):
    """
    Reads in a Swiftest config.in file and saves it as a dictionary

    Parameters
    ----------
    config_file_name : string
        File name of the input parameter file

    Returns
    -------
    config : dict
        A dictionary containing the entries in the user parameter file
    """
    config = {
    'CONFIG_FILE_NAME' : config_file_name,
    'NPLMAX'         : -1,
    'NTPMAX'         : -1,
    'T0'             : 0.0,
    'TSTOP'          : 0.0,
    'DT'             : 0.0,
    'PL_IN'          : "",
    'TP_IN'          : "",
    'IN_TYPE'        : "ASCII",
    'ISTEP_OUT'      : -1,
    'BIN_OUT'        : "",
    'OUT_TYPE'       : 'REAL8',
    'OUT_FORM'       : "XV",
    'OUT_STAT'       : "NEW",
    'ISTEP_DUMP'     : -1,
    'J2'             : 0.0,
    'J4'             : 0.0,
    'CHK_RMIN'       : -1.0,
    'CHK_RMAX'       : -1.0,
    'CHK_EJECT'      : -1.0,
    'CHK_QMIN'       : -1.0,
    'CHK_QMIN_COORD' : "HELIO",
    'CHK_QMIN_RANGE' : "",
    'QMIN_ALO'       : -1.0,
    'QMIN_AHI'       : -1.0,
    'ENC_OUT'        : "",
    'MTINY'          : -1.0,
    'MU2KG'          : -1.0,
    'TU2S'           : -1.0,
    'DU2M'           : -1.0,
    'GU'             : -1.0,
    'INV_C2'         : -1.0,
    'EXTRA_FORCE'    : 'NO',
    'BIG_DISCARD'    : 'NO',
    'CHK_CLOSE'      : 'NO',
    'FRAGMENTATION'  : 'NO',
    'MTINY_SET'      : 'NO',
    'ROTATION'       : 'NO',
    'TIDES'          : 'NO',
    'ENERGY'         : 'NO',
    'GR'             : 'NO',
    'YARKOVSKY'      : 'NO',
    'YORP'           : 'NO',
    'EUCL_THRESHOLD' : 1000000,
    }

    # Read config.in file
    print(f'Reading Swiftest file {config_file_name}' )
    f = open(config_file_name, 'r')
    swiftestlines = f.readlines()
    f.close()
    for line in swiftestlines:
        fields = line.split()
        if len(fields) > 0:
            for key in config:
                if (key == fields[0].upper()): config[key] = fields[1]
            #Special case of CHK_QMIN_RANGE requires a second input
            if (config['CHK_QMIN_RANGE'] == fields[0].upper()):
                config['QMIN_ALO'] = fields[1]
                config['QMIN_AHI'] = fields[2]

    config['NPLMAX']     = int(config['NPLMAX'])
    config['NTPMAX']     = int(config['NTPMAX'])
    config['ISTEP_OUT']  = int(config['ISTEP_OUT'])
    config['ISTEP_DUMP'] = int(config['ISTEP_DUMP'])
    config['EUCL_THRESHOLD'] = int(config['EUCL_THRESHOLD'])
    config['T0']         = float(config['T0'])
    config['TSTOP']      = float(config['TSTOP'])
    config['DT']         = float(config['DT'])
    config['J2']         = float(config['J2'])
    config['J4']         = float(config['J4'])
    config['CHK_RMIN']   = float(config['CHK_RMIN'])
    config['CHK_RMAX']   = float(config['CHK_RMAX'])
    config['CHK_EJECT']  = float(config['CHK_EJECT'])
    config['CHK_QMIN']   = float(config['CHK_QMIN'])
    config['QMIN_ALO']   = float(config['QMIN_ALO'])
    config['QMIN_AHI']   = float(config['QMIN_AHI'])
    config['MTINY']      = float(config['MTINY'])
    config['DU2M']       = float(config['DU2M'])
    config['MU2KG']      = float(config['MU2KG'])
    config['TU2S']       = float(config['TU2S'])
    config['EXTRA_FORCE'] = config['EXTRA_FORCE'].upper()
    config['BIG_DISCARD'] = config['BIG_DISCARD'].upper()
    config['CHK_CLOSE']   = config['CHK_CLOSE'].upper()
    config['FRAGMENTATION'] = config['FRAGMENTATION'].upper()

    config['GU']         = GC / (config['DU2M']**3 / (config['MU2KG'] * config['TU2S']**2))
    return config

def swifter_stream(f, param):
    """
    Reads in a Swifter bin.dat file and returns a single frame of data as a datastream

    Parameters
    ----------
    f : file object
    param : dict

    Yields
    -------
    t    : float
        Time of this frame
    npl  : int
        Number of massive bodies
    plid : int array
        IDs of massive bodies
    pvec : float array
        (npl,N) - vector of N quantities or each particle (6 of XV/EL + Mass, Radius, etc)
    plab : string list
        Labels for the pvec data
    ntp  : int
        Number of test particles
    tpid : int array
        Ids of test particles
    tvec : float array
        (ntp,N) - vector of N quantities for each particle (6 of XV/EL, etc.)
    tlab : string list
        Labels for the tvec data
    """
    while True:  # Loop until you read the end of file
        try:
            # Read single-line header
            record =  f.read_record('<f8', '<i4', '<i4', '<i4')
        except:
            break
        t = record[0]
        npl = record[1][0] - 1
        ntp = record[2][0]

        pvec = np.empty((6, npl))
        plid = np.empty(npl, dtype='int')
        tvec = np.empty((6, ntp))
        tpid = np.empty(ntp, dtype='int')
        if npl > 0:
            Mpl = np.empty(npl)
            Rpl = np.empty(npl)
            for i in range(npl):
                #Read single-line pl frame for
                record =  f.read_record('<i4', '<f8', '<f8', '(6,)<f8')
                plid[i] = record[0]
                Mpl[i] = record[1]
                Rpl[i] = record[2]
                pvec[:,i] = record[3]
        if ntp > 0:
            for i in range(ntp):
                record =  f.read_record('<i4', '(6,)<f8')
                tpid[i] = record[0]
                tvec[:,i] = record[1]

        tlab = []
        if param['OUT_FORM'] == 'XV':
            tlab.append('px')
            tlab.append('py')
            tlab.append('pz')
            tlab.append('vx')
            tlab.append('vy')
            tlab.append('vz')
        elif param['OUT_FORM'] == 'EL':
            tlab.append('a')
            tlab.append('e')
            tlab.append('inc')
            tlab.append('capom')
            tlab.append('omega')
            tlab.append('capm')
        plab = tlab.copy()
        plab.append('Mass')
        plab.append('Radius')
        pvec = np.vstack([pvec,Mpl,Rpl])

        yield t, npl, plid, pvec.T, plab, \
              ntp, tpid, tvec.T, tlab

def swiftest_stream(f, config):
    """
    Reads in a Swifter bin.dat file and returns a single frame of data as a datastream

    Parameters
    ----------
    f : file object
    param : dict

    Yields
    -------
    t    : float
        Time of this frame
    npl  : int
        Number of massive bodies
    plid : int array
        IDs of massive bodies
    pvec : float array
        (npl,N) - vector of N quantities or each particle (6 of XV/EL + Mass, Radius, etc)
    plab : string list
        Labels for the pvec data
    ntp  : int
        Number of test particles
    tpid : int array
        Ids of test particles
    tvec : float array
        (ntp,N) - vector of N quantities for each particle (6 of XV/EL, etc.)
    tlab : string list
        Labels for the tvec data
    """
    while True:  # Loop until you read the end of file
        try:
            # Read multi-line header
            t = f.read_reals(np.float64)  # Try first part of the header
        except:
            break
        npl = f.read_ints()
        ntp = f.read_ints()
        iout_form = f.read_reals('c')
        if npl[0] > 0:
            plid = f.read_ints()
            p1 = f.read_reals(np.float64)
            p2 = f.read_reals(np.float64)
            p3 = f.read_reals(np.float64)
            p4 = f.read_reals(np.float64)
            p5 = f.read_reals(np.float64)
            p6 = f.read_reals(np.float64)
            Mpl = f.read_reals(np.float64)
            Rpl = f.read_reals(np.float64)
        if ntp[0] > 0:
            tpid = f.read_ints()
            t1 = f.read_reals(np.float64)
            t2 = f.read_reals(np.float64)
            t3 = f.read_reals(np.float64)
            t4 = f.read_reals(np.float64)
            t5 = f.read_reals(np.float64)
            t6 = f.read_reals(np.float64)

        tlab = []
        if config['OUT_FORM'] == 'XV':
            tlab.append('px')
            tlab.append('py')
            tlab.append('pz')
            tlab.append('vx')
            tlab.append('vy')
            tlab.append('vz')
        elif config['OUT_FORM'] == 'EL':
            tlab.append('a')
            tlab.append('e')
            tlab.append('inc')
            tlab.append('capom')
            tlab.append('omega')
            tlab.append('capm')
        plab = tlab.copy()
        plab.append('Mass')
        plab.append('Radius')

        if npl > 0:
            pvec = np.vstack([p1,p2,p3,p4,p5,p6,Mpl,Rpl])
        else:
            pvec = np.empty((8,0))
            plid = np.empty(0)
        if ntp > 0:
            tvec = np.vstack([t1,t2,t3,t4,t5,t6])
        else:
            tvec = np.empty((6,0))
            tpid = np.empty(0)
        yield t, npl, plid, pvec.T, plab, \
              ntp, tpid, tvec.T, tlab

def swifter2xr(param):
    dims  = ['time','id', 'vec']
    pl = []
    tp = []
    nplarr = []
    ntparr = []
    with FortranFile(param['BIN_OUT'], 'r') as f:
        for t, npl, plid, pvec, plab, \
              ntp, tpid, tvec, tlab in swifter_stream(f, param):

            #Prepare frames by adding an extra axis for the time coordinate
            plframe = np.expand_dims(pvec, axis=0)
            tpframe = np.expand_dims(tvec, axis=0)

            #Create xarray DataArrays out of each body type
            plxr = xr.DataArray(plframe, dims = dims, coords = {'time' : t, 'id' : plid, 'vec' : plab})
            tpxr = xr.DataArray(tpframe, dims = dims, coords = {'time' : t, 'id' : tpid, 'vec' : tlab})
            plda = xr.concat(pl, dim='time')
            tpda = xr.concat(tp, dim='time')
            plds = plda.to_dataset(dim='vec')
            tpds = tpda.to_dataset(dim='vec')
            dsi = xr.combine_by_coords([plds, tpds])
            #nplarr.append(npl)
            #ntparr.append(ntp)

            pl.append(plxr)
            tp.append(tpxr)




        ds = xr.combine_by_coords([plds, tpds])

    return ds

def swiftest2xr(config, subsize):
    """Reads in the """

    dims  = ['time','id', 'vec']
    dsframes = []
    numsubs = 0

    subdat = 0
    filehead = 'bin'

    tmpdir = tempfile.mkdtemp()

    # Ensure the file is read/write by the creator only

    #tmppath = os.path.join(tmpdir, predictable_filename)

    with FortranFile(config['BIN_OUT'], 'r') as f:
        for t, npl, plid, pvec, plab, \
            ntp, tpid, tvec, tlab in swiftest_stream(f, config):
            print(f'Time = {t[0]}')
            plframe = np.expand_dims(pvec, axis=0)
            tpframe = np.expand_dims(tvec, axis=0)
            bd = []

            # Create xarray DataArrays out of each body type
            if npl[0] > 0:
                plxr = xr.DataArray(plframe, dims=dims, coords={'time': t, 'id': plid, 'vec': plab})
                bd.append(plxr)
            if ntp[0] > 0:
                tpxr = xr.DataArray(tpframe, dims=dims, coords={'time': t, 'id': tpid, 'vec': tlab})
                bd.append(tpxr)
            bdxr = xr.concat(bd, dim='time')
            bdxr = bdxr.to_dataset(dim='vec')
            bdxr = bdxr.assign(npl=npl[0])
            bdxr = bdxr.assign(ntp=ntp[0])
            bdxr = bdxr.chunk()
            dsframes.append(bdxr)
            subdat += 1
            if subdat == subsize:
                numsubs += 1
                filename = f'{filehead}{numsubs:03d}.tmp.nc'
                path = os.path.join(tmpdir, filename)
                print(f'Concatenating sub Dataset to file {path}')
                xr.concat(dsframes, dim='time').to_netcdf(path)
                dsframes = []
                chunk = 0
    if (chunk > 0):
        numsubs += 1
        filename = f'{filehead}{numsubs:03d}.tmp.nc'
        path = os.path.join(tmpdir, filename)
        print(f'Concatenating sub Dataset to file {path}')
        xr.concat(dsframes, dim='time').to_netcdf(path)

    print('Combining sub Datasets into a single Dataset')
    filename = f'{filehead}'
    path = os.path.join(tmpdir, filename)
    ds = xr.open_mfdataset(f'{path}*.tmp.nc',combine='by_coords',concat_dim='time',
                           data_vars='minimal', coords='minimal', compat='override')

    os.rmdir(tmpdir)
    filename = f'{filehead}.nc'
    path = os.path.join(config['WORKINGDIR'], filename)
    print(f'Saving complete dataset to file {path}')
    ds.to_netcdf(path=path)

    return ds

if __name__ == '__main__':

    workingdir = '/Users/daminton/work/Projects/Swiftest/Elliott_Performance/high_high_1500_1/'

    config_file_name = workingdir + 'param.in'
    config = read_swiftest_config(config_file_name)
    config['BIN_OUT'] = workingdir + config['BIN_OUT']
    config['WORKINGDIR'] = workingdir

    swiftestdat = swiftest2xr(config, subsize=10)

