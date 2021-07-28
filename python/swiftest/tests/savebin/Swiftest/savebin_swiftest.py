import swiftest
import swiftest.io as swio

param_file = "param.swiftest.in"
sim = swiftest.Simulation(param_file="param.swiftest.in")
ds = swio.swiftest2xr(sim.param)
swio.savebin(sim.ds, outname="testfile.nc")