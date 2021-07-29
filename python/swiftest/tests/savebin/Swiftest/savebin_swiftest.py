import swiftest

sim = swiftest.Simulation(param_file="param.in")
ds = sim.bin2xr()
sim.savebin(sim.ds, outname="testfile.nc")
