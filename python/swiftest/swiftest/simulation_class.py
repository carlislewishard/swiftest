from swiftest import io
from swiftest import init_cond
from swiftest import tool
from swiftest import constants
from datetime import date
import xarray as xr
import numpy as np
import os
import shutil

class Simulation:
    """
    This is a class that defines the basic Swift/Swifter/Swiftest simulation object
    """
    def __init__(self, codename="Swiftest", param_file="param.in", readbin=True, verbose=True):
        self.ds = xr.Dataset()
        self.param = {
            '! VERSION': f"Swiftest parameter input",
            'T0': "0.0",
            'TSTOP': "0.0",
            'DT': "0.0",
            'PL_IN': "pl.in",
            'TP_IN': "tp.in",
            'CB_IN': "cb.in",
            'IN_TYPE': "ASCII",
            'IN_FORM': "EL",
            'ISTEP_OUT': "1",
            'ISTEP_DUMP': "1",
            'BIN_OUT': "bin.nc",
            'OUT_TYPE': 'NETCDF_DOUBLE',
            'OUT_FORM': "XVEL",
            'OUT_STAT': "REPLACE",
            'CHK_RMAX': "-1.0",
            'CHK_EJECT': "-1.0",
            'CHK_RMIN': "-1.0",
            'CHK_QMIN': "-1.0",
            'CHK_QMIN_COORD': "HELIO",
            'CHK_QMIN_RANGE': "-1.0 -1.0",
            'ENC_OUT': "",
            'MU2KG': constants.MSun,
            'TU2S': constants.JD2S,
            'DU2M': constants.AU2M,
            'EXTRA_FORCE': "NO",
            'DISCARD_OUT': "",
            'PARTICLE_OUT' : "",
            'BIG_DISCARD': "NO",
            'CHK_CLOSE': "YES",
            'RHILL_PRESENT': "YES",
            'FRAGMENTATION': "NO",
            'ROTATION': "NO",
            'TIDES': "NO",
            'ENERGY': "NO",
            'GR': "YES",
            'INTERACTION_LOOPS': "ADAPTIVE",
            'ENCOUNTER_CHECK': "ADAPTIVE"
        }
        self.codename = codename
        self.verbose = verbose
        if param_file != "" :
            dir_path = os.path.dirname(os.path.realpath(param_file))
            self.read_param(param_file, codename=codename, verbose=self.verbose)
            if readbin:
                binpath = os.path.join(dir_path,self.param['BIN_OUT'])
                if os.path.exists(binpath):
                    self.param['BIN_OUT'] = binpath
                    self.bin2xr()
                else:
                    print(f"BIN_OUT file {self.param['BIN_OUT']} not found.")
        return
    
    def add(self, plname, date=date.today().isoformat(), idval=None):
        """
        Adds a solar system body to an existing simulation DataSet.
        
        Parameters
        ----------
           plname : string
                Name of planet to add (e.g. "Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune"
           date : string
                 Date to use when obtaining the ephemerides in the format YYYY-MM-DD. Defaults to "today"
        Returns
        -------
            self.ds : xarray dataset
        """
        self.ds = init_cond.solar_system_horizons(plname, idval, self.param, date, self.ds)
        return
    
    
    def addp(self, idvals, namevals, t1, t2, t3, t4, t5, t6, GMpl=None, Rpl=None, rhill=None, Ip1=None, Ip2=None, Ip3=None, rotx=None, roty=None, rotz=None, t=None):
        """
        Adds a body (test particle or massive body) to the internal DataSet given a set up 6 vectors (orbital elements
        or cartesian state vectors, depending on the value of self.param). Input all angles in degress

        Parameters
        ----------
            idvals : integer 
                Array of body index values.
            t1     : float
                xh for param['IN_FORM'] == "XV"; a for param['IN_FORM'] == "EL"
            t2     : float
                yh for param['IN_FORM'] == "XV"; e for param['IN_FORM'] == "EL"
            t3     : float
                zh for param['IN_FORM'] == "XV"; inc for param['IN_FORM'] == "EL"
            t4     : float
                vhxh for param['IN_FORM'] == "XV"; capom for param['IN_FORM'] == "EL"
            t5     : float
                vhyh for param['IN_FORM'] == "XV"; omega for param['IN_FORM'] == "EL"
            t6     : float
                vhzh for param['IN_FORM'] == "XV"; capm for param['IN_FORM'] == "EL"
            Gmass  : float
                Optional: Array of G*mass values if these are massive bodies
            radius : float
                Optional: Array radius values if these are massive bodies
            rhill  : float
                Optional: Array rhill values if these are massive bodies
            Ip1,y,z : float
                Optional: Principal axes moments of inertia
            rotx,y,z: float
                Optional: Rotation rate vector components
            t      :  float
                Optional: Time at start of simulation
        Returns
        -------
            self.ds : xarray dataset
        """
        if t is None:
            t = self.param['T0']

        dsnew = init_cond.vec2xr(self.param, idvals, namevals, t1, t2, t3, t4, t5, t6, GMpl, Rpl, rhill, Ip1, Ip2, Ip3, rotx, roty, rotz, t)
        if dsnew is not None:
            self.ds = xr.combine_by_coords([self.ds, dsnew])
        self.ds['ntp'] = self.ds['id'].where(np.isnan(self.ds['Gmass'])).count(dim="id")
        self.ds['npl'] = self.ds['id'].where(np.invert(np.isnan(self.ds['Gmass']))).count(dim="id") - 1

        return
    
    
    def read_param(self, param_file, codename="Swiftest", verbose=True):
        """
        Reads in a param.in file and determines whether it is a Swift/Swifter/Swiftest parameter file.
        
        Parameters
        ----------
           param_file : string
                File name of the input parameter file
           codename : string
                 Type of parameter file, either "Swift", "Swifter", or "Swiftest"
        Returns
        -------
            self.ds : xarray dataset
        """
        if codename == "Swiftest":
            self.param = io.read_swiftest_param(param_file, self.param, verbose=verbose)
            self.codename = "Swiftest"
        elif codename == "Swifter":
            self.param = io.read_swifter_param(param_file, verbose=verbose)
            self.codename = "Swifter"
        elif codename == "Swift":
            self.param = io.read_swift_param(param_file, verbose=verbose)
            self.codename = "Swift"
        else:
            print(f'{codename} is not a recognized code name. Valid options are "Swiftest", "Swifter", or "Swift".')
            self.codename = "Unknown"
        return
    
    
    def write_param(self, param_file, param=None):
        """
        Writes to a param.in file and determines whether the output format needs to be converted between Swift/Swifter/Swiftest.
        
        Parameters
        ----------
           param_file : string
                File name of the input parameter file
        Returns
        -------
            self.ds : xarray dataset
        """
        if param is None:
            param = self.param
        # Check to see if the parameter type matches the output type. If not, we need to convert
        codename = param['! VERSION'].split()[0]
        if codename == "Swifter" or codename == "Swiftest":
            io.write_labeled_param(param, param_file)
        elif codename == "Swift":
            io.write_swift_param(param, param_file)
        else:
            print('Cannot process unknown code type. Call the read_param method with a valid code name. Valid options are "Swiftest", "Swifter", or "Swift".')
        return
    
    
    def convert(self, param_file, newcodename="Swiftest", plname="pl.swiftest.in", tpname="tp.swiftest.in", cbname="cb.swiftest.in", conversion_questions={}):
        """
        Converts simulation input files from one format to another (Swift, Swifter, or Swiftest). 

        Parameters
        ----------
           param_file : string
                File name of the input parameter file
            newcodename : string
                Name of the desired format (Swift/Swifter/Swiftest)
            plname : string
                File name of the massive body input file
            tpname : string
                File name of the test particle input file
            cbname : string
                File name of the central body input file
            conversion_questions : dictronary
                Dictionary of additional parameters required to convert between formats

        Returns
        -------
            oldparam : xarray dataset
                The old parameter configuration.
        """
        oldparam = self.param
        if self.codename == newcodename:
            print(f"This parameter configuration is already in {newcodename} format")
            return oldparam
        if newcodename != "Swift" and newcodename != "Swifter" and newcodename != "Swiftest":
            print(f'{newcodename} is an invalid code type. Valid options are "Swiftest", "Swifter", or "Swift".')
            return oldparam
        goodconversion = True
        if self.codename == "Swifter":
            if newcodename == "Swiftest":
                self.param = io.swifter2swiftest(self.param, plname, tpname, cbname, conversion_questions)
            else:
                goodconversion = False
        elif self.codename == "Swift":
            if newcodename == "Swifter":
                self.param = io.swift2swifter(self.param, plname, tpname, conversion_questions)
            elif newcodename == "Swiftest":
                self.param = io.swift2swiftest(self.param, plname, tpname, cbname, conversion_questions)
            else:
                goodconversion = False
        else:
            goodconversion = False
            
        if goodconversion:
            self.write_param(param_file)
        else:
            print(f"Conversion from {self.codename} to {newcodename} is not supported.")
        return oldparam
    
    
    def bin2xr(self):
        """
        Converts simulation output files from a flat binary file to a xarray dataset. 

        Parameters
        ----------

        Returns
        -------
            self.ds : xarray dataset
        """
        if self.codename == "Swiftest":
            self.ds = io.swiftest2xr(self.param, verbose=self.verbose)
            if self.verbose: print('Swiftest simulation data stored as xarray DataSet .ds')
        elif self.codename == "Swifter":
            self.ds = io.swifter2xr(self.param, verbose=self.verbose)
            if self.verbose: print('Swifter simulation data stored as xarray DataSet .ds')
        elif self.codename == "Swift":
            print("Reading Swift simulation data is not implemented yet")
        else:
            print('Cannot process unknown code type. Call the read_param method with a valid code name. Valid options are "Swiftest", "Swifter", or "Swift".')
        return
    
    
    def follow(self, codestyle="Swifter"):
        """
        An implementation of the Swift tool_follow algorithm. Under development. Currently only for Swift simulations. 

        Parameters
        ----------
            codestyle : string
                Name of the desired format (Swift/Swifter/Swiftest)

        Returns
        -------
            fol : xarray dataset
        """
        if self.ds is None:
            self.bin2xr()
        if codestyle == "Swift":
            try:
                with open('follow.in', 'r') as f:
                    line = f.readline() # Parameter file (ignored because bin2xr already takes care of it
                    line = f.readline() # PL file (ignored)
                    line = f.readline() # TP file (ignored)
                    line = f.readline() # ifol
                    i_list = [i for i in line.split(" ") if i.strip()]
                    ifol = int(i_list[0])
                    line = f.readline()  # nskp
                    i_list = [i for i in line.split(" ") if i.strip()]
                    nskp = int(i_list[0])
            except IOError:
                print('No follow.in file found')
                ifol = None
                nskp = None
            fol = tool.follow_swift(self.ds, ifol=ifol, nskp=nskp)
        else:
            fol = None
        
        if self.verbose: print('follow.out written')
        return fol
    
    
    def save(self, param_file, framenum=-1, codename="Swiftest"):
        """
        Saves an xarray dataset to a set of input files.

        Parameters
        ----------
            param_file : string
                Name of the parameter input file
            framenum : integer (default=-1)
                Time frame to use to generate the initial conditions. If this argument is not passed, the default is to use the last frame in the dataset.
            codename : string
                Name of the desired format (Swift/Swifter/Swiftest)

        Returns
        -------
            self.ds : xarray dataset
        """

        if codename == "Swiftest":
            io.swiftest_xr2infile(self.ds, self.param, framenum)
            self.write_param(param_file)
        elif codename == "Swifter":
            if self.codename == "Swiftest":
                swifter_param = io.swiftest2swifter_param(self.param)
            else:
                swifter_param = self.param
            io.swifter_xr2infile(self.ds, swifter_param, framenum)
            self.write_param(param_file, param=swifter_param)
        else:
            print(f'Saving to {codename} not supported')

        return

    def initial_conditions_from_bin(self, framenum=-1, new_param=None, new_param_file="param.new.in", new_initial_conditions_file="bin_in.nc",  restart=False, codename="Swiftest"):
        """
        Generates a set of input files from a old output file.

        Parameters
        ----------
            framenum : integer (default=-1)
                Time frame to use to generate the initial conditions. If this argument is not passed, the default is to use the last frame in the dataset.
            new_param : string
                File to copy parameters from. Default is the old parameter file.
            new_param_file : string
                Name of the new parameter file.
            new_initial_conditions_file : string
                Name of the new NetCDF file containing the new initial conditions.
            restart : True or False
                If True, overwrite the old output file. If False, generate a new output file.
            codename : string
                Name of the desired format (Swift/Swifter/Swiftest)

        Returns
        -------
            frame : NetCDF dataset 
        """


        if codename != "Swiftest":
            self.save(new_param_file, framenum, codename)
            return
        if new_param is None:
            new_param = self.param.copy()

        if codename == "Swiftest":
            if restart:
                new_param['T0'] = self.ds.time.values[framenum]
            if self.param['OUT_TYPE'] == 'NETCDF_DOUBLE' or self.param['OUT_TYPE'] == 'REAL8':
                new_param['IN_TYPE'] = 'NETCDF_DOUBLE'
            elif self.param['OUT_TYPE'] == 'NETCDF_FLOAT' or self.param['OUT_TYPE'] == 'REAL4':
                new_param['IN_TYPE'] = 'NETCDF_FLOAT'
            else:
                print(f"{self.param['OUT_TYPE']} is an invalid OUT_TYPE file")
                return
            if self.param['BIN_OUT'] != new_param['BIN_OUT'] and restart:
               print(f"Restart run with new output file. Copying {self.param['BIN_OUT']} to {new_param['BIN_OUT']}")
               shutil.copy2(self.param['BIN_OUT'],new_param['BIN_OUT'])
            new_param['IN_FORM'] = 'XV'
            if restart:
                new_param['OUT_STAT'] = 'APPEND'
            new_param['FIRSTKICK'] = 'T'
            new_param['NC_IN'] = new_initial_conditions_file
            new_param.pop('PL_IN', None)
            new_param.pop('TP_IN', None)
            new_param.pop('CB_IN', None)
            print(f"Extracting data from dataset at time frame number {framenum} and saving it to {new_param['NC_IN']}")
            frame = io.swiftest_xr2infile(self.ds, self.param, infile_name=new_param['NC_IN'],framenum=framenum)
            print(f"Saving parameter configuration file to {new_param_file}")
            self.write_param(new_param_file, param=new_param)

        return frame
