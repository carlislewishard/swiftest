import numpy as np
import pandas as pd
from scipy.io import FortranFile
import xarray as xr
from astroquery.jplhorizons import Horizons
import astropy.constants as const
import datetime
import sys
from swiftestpy import orbel

# Constants in SI units
AU2M = np.longdouble(const.au.value)
GMSunSI = np.longdouble(const.GM_sun.value)
RSun = np.longdouble(const.R_sun.value)
GC = np.longdouble(const.G.value)
JD2S = 86400
year = np.longdouble(365.25 * JD2S)
einsteinC = np.longdouble(299792458.0)
# Solar oblatenes values: From Mecheri et al. (2004), using Corbard (b) 2002 values (Table II)
J2Sun = np.longdouble(2.198e-7)
J4Sun = np.longdouble(-4.805e-9)

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
   'INPARFILE'     : inparfile,
   'NPLMAX'       : -1,
   'NTPMAX'       : -1,
   'T0'          : 0.0,
   'TSTOP'        : 0.0,
   'DT'          : 0.0,
   'PL_IN'        : "",
   'TP_IN'        : "",
   'IN_TYPE'      : "ASCII",
   'ISTEP_OUT'     : -1,
   'BIN_OUT'      : "",
   'OUT_TYPE'      : 'REAL8',
   'OUT_FORM'      : "XV",
   'OUT_STAT'      : "NEW",
   'ISTEP_DUMP'    : -1,
   'J2'          : 0.0,
   'J4'          : 0.0,
   'CHK_CLOSE'     : 'NO',
   'CHK_RMIN'      : -1.0,
   'CHK_RMAX'      : -1.0,
   'CHK_EJECT'     : -1.0,
   'CHK_QMIN'      : -1.0,
   'CHK_QMIN_COORD' : "HELIO",
   'CHK_QMIN_RANGE' : "",
   'QMIN_ALO'      : -1.0,
   'QMIN_AHI'      : -1.0,
   'ENC_OUT'      : "",
   'EXTRA_FORCE'   : 'NO',
   'BIG_DISCARD'   : 'NO',
   'RHILL_PRESENT'  : 'NO',
   'GR'          : 'NO',
   'C2'          : -1.0,
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

   param['NPLMAX']    = int(param['NPLMAX'])
   param['NTPMAX']    = int(param['NTPMAX'])
   param['ISTEP_OUT']  = int(param['ISTEP_OUT'])
   param['ISTEP_DUMP'] = int(param['ISTEP_DUMP'])
   param['T0']       = float(param['T0'])
   param['TSTOP']     = float(param['TSTOP'])
   param['DT']       = float(param['DT'])
   param['J2']       = float(param['J2'])
   param['J4']       = float(param['J4'])
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
   'NPLMAX'       : -1,
   'NTPMAX'       : -1,
   'T0'          : 0.0,
   'TSTOP'        : 0.0,
   'DT'          : 0.0,
   'PL_IN'        : "",
   'TP_IN'        : "",
   'IN_TYPE'      : "ASCII",
   'ISTEP_OUT'    : -1,
   'BIN_OUT'      : "",
   'PARTICLE_FILE' : "",
   'OUT_TYPE'      : 'REAL8',
   'OUT_FORM'      : "XV",
   'OUT_STAT'      : "NEW",
   'ISTEP_DUMP'    : -1,
   'J2'          : 0.0,
   'J4'          : 0.0,
   'CHK_RMIN'      : -1.0,
   'CHK_RMAX'      : -1.0,
   'CHK_EJECT'     : -1.0,
   'CHK_QMIN'      : -1.0,
   'CHK_QMIN_COORD' : "HELIO",
   'CHK_QMIN_RANGE' : "",
   'QMIN_ALO'      : -1.0,
   'QMIN_AHI'      : -1.0,
   'ENC_OUT'      : "",
   'MTINY'        : -1.0,
   'MU2KG'        : -1.0,
   'TU2S'         : -1.0,
   'DU2M'         : -1.0,
   'GU'          : -1.0,
   'INV_C2'       : -1.0,
   'EXTRA_FORCE'   : 'NO',
   'BIG_DISCARD'   : 'NO',
   'CHK_CLOSE'     : 'NO',
   'FRAGMENTATION'  : 'NO',
   'MTINY_SET'     : 'NO',
   'ROTATION'      : 'NO',
   'TIDES'        : 'NO',
   'ENERGY'       : 'NO',
   'GR'          : 'NO',
   'YARKOVSKY'     : 'NO',
   'YORP'         : 'NO',
   'EUCL_THRESHOLD' : 1000000,
   'SKIP'         : 0,
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

   config['NPLMAX']    = int(config['NPLMAX'])
   config['NTPMAX']    = int(config['NTPMAX'])
   config['ISTEP_OUT']  = int(config['ISTEP_OUT'])
   config['ISTEP_DUMP'] = int(config['ISTEP_DUMP'])
   config['EUCL_THRESHOLD'] = int(config['EUCL_THRESHOLD'])
   config['T0']       = float(config['T0'])
   config['TSTOP']     = float(config['TSTOP'])
   config['DT']       = float(config['DT'])
   config['J2']       = float(config['J2'])
   config['J4']       = float(config['J4'])
   config['CHK_RMIN']   = float(config['CHK_RMIN'])
   config['CHK_RMAX']   = float(config['CHK_RMAX'])
   config['CHK_EJECT']  = float(config['CHK_EJECT'])
   config['CHK_QMIN']   = float(config['CHK_QMIN'])
   config['QMIN_ALO']   = float(config['QMIN_ALO'])
   config['QMIN_AHI']   = float(config['QMIN_AHI'])
   config['MTINY']     = float(config['MTINY'])
   config['DU2M']      = float(config['DU2M'])
   config['MU2KG']     = float(config['MU2KG'])
   config['TU2S']      = float(config['TU2S'])
   config['EXTRA_FORCE'] = config['EXTRA_FORCE'].upper()
   config['BIG_DISCARD'] = config['BIG_DISCARD'].upper()
   config['CHK_CLOSE']   = config['CHK_CLOSE'].upper()
   config['FRAGMENTATION'] = config['FRAGMENTATION'].upper()
   config['ROTATION'] = config['ROTATION'].upper()

   config['GU']       = GC / (config['DU2M']**3 / (config['MU2KG'] * config['TU2S']**2))
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
   t   : float
      Time of this frame
   npl  : int
      Number of massive bodies
   plid : int array
      IDs of massive bodies
   pvec : float array
      (npl,N) - vector of N quantities or each particle (6 of XV/EL + mass, radius, etc)
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
      plab.append('mass')
      plab.append('radius')
      pvec = np.vstack([pvec,Mpl,Rpl])

      yield t, npl, plid, pvec.T, plab, \
           ntp, tpid, tvec.T, tlab

def make_swiftest_labels(config):
   tlab = []
   if config['OUT_FORM'] == 'XV':
      tlab.append('px')
      tlab.append('py')
      tlab.append('pz')
      tlab.append('vx')
      tlab.append('vy')
      tlab.append('vz')
      tlab.append('a')
      tlab.append('e')
      tlab.append('inc')
      tlab.append('capom')
      tlab.append('omega')
      tlab.append('capm')
   elif config['OUT_FORM'] == 'EL':
      tlab.append('a')
      tlab.append('e')
      tlab.append('inc')
      tlab.append('capom')
      tlab.append('omega')
      tlab.append('capm')
   plab = tlab.copy()
   plab.append('mass')
   plab.append('radius')
   if config['ROTATION'] == 'YES':
      plab.append('rot_x')
      plab.append('rot_y')
      plab.append('rot_z')
      plab.append('Ip_x')
      plab.append('Ip_y')
      plab.append('Ip_z')
   if config['TIDES'] == 'YES':
      plab.append('k2')
      plab.append('Q')
   return plab, tlab


def swiftest_stream(f, config):
   """
   Reads in a Swiftest bin.dat file and returns a single frame of data as a datastream

   Parameters
   ----------
   f : file object
   config : dict

   Yields
   -------
   t   : float
      Time of this frame
   npl  : int
      Number of massive bodies
   plid : int array
      IDs of massive bodies
   pvec : float array
      (npl,N) - vector of N quantities or each particle (6 of XV/EL + mass, radius, etc)
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
         if config['ROTATION'] == 'YES':
            rot_x = f.read_reals(np.float64)
            rot_y = f.read_reals(np.float64)
            rot_z = f.read_reals(np.float64)
            Ip_1  = f.read_reals(np.float64)
            Ip_2  = f.read_reals(np.float64)
            Ip_3  = f.read_reals(np.float64)
      if ntp[0] > 0:
         tpid = f.read_ints()
         t1 = f.read_reals(np.float64)
         t2 = f.read_reals(np.float64)
         t3 = f.read_reals(np.float64)
         t4 = f.read_reals(np.float64)
         t5 = f.read_reals(np.float64)
         t6 = f.read_reals(np.float64)

      plab, tlab = make_swiftest_labels(config)

      if config['OUT_FORM'] == 'XV':
         mu = np.empty_like(p1)
         mu[0] = Mpl[0]
         mu[1:] = Mpl[0] + Mpl[1:]
         p7 = []
         p8 = []
         p9 = []
         p10 = []
         p11 = []
         p12 = []
         for i in range(mu.size):
            elem = orbel.xv2el(mu[i], p1[i], p2[i], p3[i], p4[i], p5[i], p6[i])
            p7.append(elem[0])
            p8.append(elem[1])
            p9.append(elem[2])
            p10.append(elem[3])
            p11.append(elem[4])
            p12.append(elem[5])
         p7 = np.array(p7)
         p8 = np.array(p8)
         p9 = np.array(p9)
         p10 = np.array(p10)
         p11 = np.array(p11)
         p12 = np.array(p12)
         if ntp[0] > 0:
            mu = np.full_like(t1,Mpl[0])
            t7 = []
            t8 = []
            t9 = []
            t10 = []
            t11 = []
            t12 = []
            for i in range(mu.size):
               elem = orbel.xv2el(mu[i], t1[i], t2[i], t3[i], t4[i], t5[i], t6[i])
               t7.append(elem[0])
               t8.append(elem[1])
               t9.append(elem[2])
               t10.append(elem[3])
               t11.append(elem[4])
               t12.append(elem[5])
            t7 = np.array(t7)
            t8 = np.array(t8)
            t9 = np.array(t9)
            t10 = np.array(t10)
            t11 = np.array(t11)
            t12 = np.array(t12)

      if npl > 0:
         if config['ROTATION'] == 'YES':
            if config['OUT_FORM'] == 'EL':
               pvec = np.vstack([p1,p2,p3,p4,p5,p6,Mpl,Rpl,rot_x,rot_y,rot_z,Ip_1,Ip_2,Ip_3])
            else:
               pvec = np.vstack([p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, Mpl, Rpl, rot_x, rot_y, rot_z, Ip_1, Ip_2, Ip_3])
         else:
            if config['OUT_FORM'] == 'EL':
               pvec = np.vstack([p1,p2,p3,p4,p5,p6,Mpl,Rpl])
            else:
               pvec = np.vstack([p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, Mpl, Rpl])
      else:
         pvec = np.empty((8,0))
         plid = np.empty(0)
      if ntp > 0:
         if config['OUT_FORM'] == 'EL':
            tvec = np.vstack([t1,t2,t3,t4,t5,t6])
         else:
            tvec = np.vstack([t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12])
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

def swiftest2xr(config):
   """Reads in the Swiftest BIN_OUT file and converts it to an xarray Dataset"""
   dims  = ['time','id', 'vec']
   dsframes = []

   with FortranFile(config['BIN_OUT'], 'r') as f:
      for t, npl, plid, pvec, plab, \
         ntp, tpid, tvec, tlab in swiftest_stream(f, config):
         plframe = np.expand_dims(pvec, axis=0)
         tpframe = np.expand_dims(tvec, axis=0)
         bd = []

         # Create xarray DataArrays out of each body type
         if npl[0] > 0 :
            plxr = xr.DataArray(plframe, dims=dims, coords={'time': t, 'id': plid, 'vec': plab})
            bd.append(plxr)
         if ntp[0] > 0 :
            tpxr = xr.DataArray(tpframe, dims=dims, coords={'time': t, 'id': tpid, 'vec': tlab})
            bd.append(tpxr)
         bdxr = xr.concat(bd, dim='time')
         bdxr = bdxr.to_dataset(dim='vec')
         bdxr = bdxr.assign(npl=npl[0])
         bdxr = bdxr.assign(ntp=ntp[0])
         dsframes.append(bdxr)
         sys.stdout.write('\r'+f"Reading in time {t[0]:.3e}")
         sys.stdout.flush()
      print('\nCreating Dataset')
      ds = xr.concat(dsframes, dim='time')
   if not config['PARTICLE_FILE'] == '':
      ds = swiftest_particle_2xr(ds, config)
   print(f"Successfully converted {ds.sizes['time']} output frames.")
   return ds

def swiftest_particle_stream(f):
   """
   Reads in a Swiftest particle.dat file and returns a single frame of particle data as a datastream

   Parameters
   ----------
   f : file object
   config : dict

   Yields
   -------
   plid : int
      ID of massive bodie
   origin_type : string
      The origin type for the body (Initial conditions, disruption, supercatastrophic, hit and run, etc)
   origin_xh : float array
      The origin heliocentric position vector
   origin_vh : float array
      The origin heliocentric velocity vector
   """
   while True:  # Loop until you read the end of file
      try:
         # Read multi-line header
         plid = f.read_ints()  # Try first part of the header
      except:
         break
      origin_rec = f.read_record(np.dtype('a32'), np.dtype(('<f8', (7))))
      origin_type = np.char.strip(str(origin_rec[0], encoding='utf-8'))
      origin_vec = origin_rec[1]
      yield plid, origin_type, origin_vec

def swiftest_particle_2xr(ds, config):
   """Reads in the Swiftest PARTICLE_FILE  and converts it to an xarray Dataset"""
   veclab = ['time_origin', 'px_origin', 'py_origin', 'pz_origin', 'vx_origin', 'vy_origin', 'vz_origin']
   id_list = []
   origin_type_list = []
   origin_vec_list = []

   with FortranFile(config['PARTICLE_FILE'], 'r') as f:
      for plid, origin_type, origin_vec in swiftest_particle_stream(f):
         id_list.append(plid)

         origin_type_list.append(origin_type)
         origin_vec_list.append(origin_vec)

   id_list =  np.asarray(id_list)[:,0]
   origin_type_list = np.asarray(origin_type_list)
   origin_vec_list = np.vstack(origin_vec_list)

   typeda = xr.DataArray(origin_type_list, dims=['id'], coords={'id' : id_list})
   vecda = xr.DataArray(origin_vec_list, dims=['id', 'vec'], coords={'id' : id_list, 'vec' : veclab})

   infoxr = vecda.to_dataset(dim='vec')
   infoxr['origin_type'] = typeda

   print('\nAdding particle info to Dataset')
   ds = xr.merge([ds, infoxr])
   return ds


def solar_system_pl(config, ephemerides_start_date):
   """
   Initializes a Swiftest dataset containing the major planets of the Solar System at a particular data from JPL/Horizons

   Parameters
   ----------
   config : dict
       Swiftest Configuration parameters. This method uses the unit conversion factors to convert from JPL's AU-day system into the system specified in the config file
   ephemerides_start_date : string
       Date to use when obtaining the ephemerides in the format YYYY-MM-DD

   Returns
   -------
   xarray dataset
   """
   # Planet ids
   planetid = {
      'mercury': '1',
      'venus': '2',
      'earthmoon': '3',
      'mars': '4',
      'jupiter': '5',
      'saturn': '6',
      'uranus': '7',
      'neptune': '8',
      'plutocharon': '9'
   }

   # Planet Msun/M ratio
   MSun_over_Mpl = {
      'mercury': np.longdouble(6023600.0),
      'venus': np.longdouble(408523.71),
      'earthmoon': np.longdouble(328900.56),
      'mars': np.longdouble(3098708.),
      'jupiter': np.longdouble(1047.3486),
      'saturn': np.longdouble(3497.898),
      'uranus': np.longdouble(22902.98),
      'neptune': np.longdouble(19412.24),
      'plutocharon': np.longdouble(1.35e8)
   }

   # Planet radii in meters
   planetradius = {
      'mercury': np.longdouble(2439.4e3),
      'venus': np.longdouble(6051.8e3),
      'earthmoon': np.longdouble(6371.0084e3),  # Earth only for radius
      'mars': np.longdouble(3389.50e3),
      'jupiter': np.longdouble(69911e3),
      'saturn': np.longdouble(58232.0e3),
      'uranus': np.longdouble(25362.e3),
      'neptune': np.longdouble(24622.e3),
      'plutocharon': np.longdouble(1188.3e3)
   }

   # Unit conversion factors
   DCONV = AU2M / config['DU2M']
   VCONV = (AU2M / JD2S) / (config['DU2M'] / config['TU2S'])
   THIRDLONG = np.longdouble(1.0) / np.longdouble(3.0)

   # Central body value vectors
   GMcb = np.array([GMSunSI * config['TU2S'] ** 2 / config['DU2M'] ** 3])
   Rcb = np.array([RSun / config['DU2M']])
   J2RP2 = np.array([J2Sun * (RSun / config['DU2M']) ** 2])
   J4RP4 = np.array([J4Sun * (RSun / config['DU2M']) ** 4])
   cbid = np.array([0])
   cvec = np.vstack([GMcb, Rcb, J2RP2, J4RP4])

   # Horizons date time internal variables
   tstart = datetime.date.fromisoformat(ephemerides_start_date)
   tstep = datetime.timedelta(days=1)
   tend = tstart + tstep
   ephemerides_end_date = tend.isoformat()
   ephemerides_step = '1d'

   pldata = {}
   p1 = []
   p2 = []
   p3 = []
   p4 = []
   p5 = []
   p6 = []
   Rhill = []
   Rpl = []
   GMpl = []

   # Fetch solar system ephemerides from Horizons
   for key, val in planetid.items():
      pldata[key] = Horizons(id=val, id_type='majorbody', location='@sun',
                             epochs={'start': ephemerides_start_date, 'stop': ephemerides_end_date,
                                     'step': ephemerides_step})
      if config['OUT_FORM'] == 'XV':
         p1.append(pldata[key].vectors()['x'][0] * DCONV)
         p2.append(pldata[key].vectors()['y'][0] * DCONV)
         p3.append(pldata[key].vectors()['z'][0] * DCONV)
         p4.append(pldata[key].vectors()['vx'][0] * VCONV)
         p5.append(pldata[key].vectors()['vy'][0] * VCONV)
         p6.append(pldata[key].vectors()['vz'][0] * VCONV)
      elif config['OUT_FORM'] == 'EL':
         p1.append(pldata[key].elements()['a'][0] * DCONV)
         p2.append(pldata[key].elements()['e'][0])
         p3.append(pldata[key].elements()['inc'][0] * np.pi / 180.0)
         p4.append(pldata[key].elements()['Omega'][0] * np.pi / 180.0)
         p5.append(pldata[key].elements()['w'][0] * np.pi / 180.0)
         p6.append(pldata[key].elements()['M'][0] * np.pi / 180.0)
      Rhill.append(pldata[key].elements()['a'][0] * (3 * MSun_over_Mpl[key]) ** (-THIRDLONG))
      Rpl.append(planetradius[key] * DCONV)
      GMpl.append(GMcb[0] / MSun_over_Mpl[key])
   # Generate planet value vectors
   plid = np.fromiter(planetid.values(), dtype=int)
   pvec = np.vstack([p1, p2, p3, p4, p5, p6, GMpl, Rpl])

   dims = ['time', 'id', 'vec']
   cb = []
   pl = []
   tp = []
   t = np.array([0.0])

   clab, plab, tlab = make_swiftest_labels(config)

   # Prepare frames by adding an extra axis for the time coordinate
   cbframe = np.expand_dims(cvec.T, axis=0)
   plframe = np.expand_dims(pvec.T, axis=0)

   # Create xarray DataArrays out of each body type
   cbxr = xr.DataArray(cbframe, dims=dims, coords={'time': t, 'id': cbid, 'vec': clab})
   plxr = xr.DataArray(plframe, dims=dims, coords={'time': t, 'id': plid, 'vec': plab})

   cb.append(cbxr)
   pl.append(plxr)

   cbda = xr.concat(cb, dim='time')
   plda = xr.concat(pl, dim='time')

   cbds = cbda.to_dataset(dim='vec')
   plds = plda.to_dataset(dim='vec')
   ds = xr.combine_by_coords([cbds, plds])
   return ds


def swiftest_xr2_infile(ds, config, framenum=-1):
   """
   Writes a set of Swiftest input files from a single frame of a Swiftest xarray dataset

   Parameters
   ----------
   ds : xarray dataset
       Dataset containing Swiftest n-body data in XV format
   framenum : int
       Time frame to use to generate the initial conditions. If this argument is not passed, the default is to use the last frame in the dataset.
   config : dict
       Swiftest Configuration parameters. This method uses the names of the cb, pl, and tp files from the configuration

   Returns
   -------
   A set of three input files for a Swiftest run
   """
   frame = ds.isel(time=framenum)
   pl = frame.where(np.invert(np.isnan(frame['Mass'])), drop=True)
   tp = frame.where(np.isnan(frame['Mass']), drop=True).drop_vars(['Mass', 'Radius'])
   tp = tp.where(np.invert(np.isnan(tp['px'])), drop=True)

   GMSun = np.double(pl['Mass'].isel(id=0))
   RSun = np.double(pl['Radius'].isel(id=0))

   if config['IN_TYPE'] == 'ASCII':
      # Swiftest PL file
      plfile = open(config['PL_IN'], 'w')
      print(pl.id.count().values, file=plfile)
      for i in pl.id:
         pli = pl.sel(id=i)
         print(i.values, pli['Mass'].values, file=plfile)
         if i != 0:
            print(pli['Radius'].values, file=plfile)
         print(pli['px'].values, pli['py'].values, pli['pz'].values, file=plfile)
         print(pli['vx'].values, pli['vy'].values, pli['vz'].values, file=plfile)
         if config['ROTATION'] == 'YES':
            print(pli['Ip_x'].values, pli['Ip_y'].values, pli['Ip_z'].values, file=plfile)
            print(pli['rot_x'].values, pli['rot_y'].values, pli['rot_z'].values, file=plfile)
      plfile.close()

      # TP file
      tpfile = open(config['TP_IN'], 'w')
      print(tp.id.count().values, file=tpfile)
      for i in tp.id:
         tpi = tp.sel(id=i)
         print(i.values, file=tpfile)
         print(tpi['px'].values, tpi['py'].values, tpi['pz'].values, file=tpfile)
         print(tpi['vx'].values, tpi['vy'].values, tpi['vz'].values, file=tpfile)
      tpfile.close()
   elif config['IN_TYPE'] == 'REAL8':
      # Now make Swiftest files
      cbfile = FortranFile(swiftest_cb, 'w')
      Msun = np.double(1.0)
      cbfile.write_record(np.double(GMSun))
      cbfile.write_record(np.double(rmin))
      cbfile.write_record(np.double(J2))
      cbfile.write_record(np.double(J4))
      cbfile.close()

      plfile = FortranFile(swiftest_pl, 'w')
      plfile.write_record(npl)

      plfile.write_record(plid)
      plfile.write_record(p_pl[0])
      plfile.write_record(p_pl[1])
      plfile.write_record(p_pl[2])
      plfile.write_record(v_pl[0])
      plfile.write_record(v_pl[1])
      plfile.write_record(v_pl[2])
      plfile.write_record(mass)
      plfile.write_record(radius)
      plfile.close()
      tpfile = FortranFile(swiftest_tp, 'w')
      ntp = 1
      tpfile.write_record(ntp)
      tpfile.write_record(tpid)
      tpfile.write_record(p_tp[0])
      tpfile.write_record(p_tp[1])
      tpfile.write_record(p_tp[2])
      tpfile.write_record(v_tp[0])
      tpfile.write_record(v_tp[1])
      tpfile.write_record(v_tp[2])
   else:
      print(f"{config['IN_TYPE']} is an unknown file type")

if __name__ == '__main__':

   workingdir = '/Users/daminton/git/swiftest/examples/symba_mars_disk/'
   config_file_name = workingdir + 'param.in'
   config = read_swiftest_config(config_file_name)
   config['BIN_OUT'] = workingdir + config['BIN_OUT']
   config['PARTICLE_FILE'] = workingdir + config['PARTICLE_FILE']
   config['WORKINGDIR'] = workingdir
   swiftestdat = swiftest2xr(config)



