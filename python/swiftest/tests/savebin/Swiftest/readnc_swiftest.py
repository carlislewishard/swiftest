import swiftest

sim = swiftest.Simulation(param_file="param.in")
ds = sim.readnc(inname="data.nc")
