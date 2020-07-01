#CALL python3 create_disk_mars.py disk_type ecc_type n_particles seed tracking output_path
# disk type low, medium, high
# ecc_type low, high
# n_particles number
# seed number between 0-9
# tracking density "yes" or "no"
# output_path such as: /Users/jpouplin/Documents/canup/
###########################################################Imports######################################################
import rebound
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import timeit
import random
import numpy as np
import openpyxl
import os
import sys
import pandas as pd
import seaborn as sns
from mpl_toolkits import mplot3d
sns.set(color_codes=True)
###############################################Argument import from shell ##############################################
start = timeit.default_timer()
disk_type = sys.argv[1]
ecc_type = sys.argv[2]
n_particles = int(sys.argv[3])
seed = sys.argv[4]
np.random.seed(int(seed))
tracking = sys.argv[5]
output_path = sys.argv[6]
dir = os.getcwd()
create = 1               #create =1 for excel spreadhseet of particles diameter in the debris disk
create_canup=1           #create = 1 for  canup input file for swifter canup_input.in
###########################################################Constants ###################################################
R_Mars=3.394e8 #Rmars in cm
R_Planet=R_Mars/100. #R_planet in m
M_Mars=6.4185e26 #MMars in grams
G_aud = 1.4878e-34 #G in AU days (kg, sec, m)
G = 6.67428e-11 #G in SI units
M_P = 10.6e15 #kg Phobos mass
R_P = 11.1e3 #radius phobos m
M_D = 1.476188e15 #kg deimos mass
R_D = 6.2e3 #radius deimos in m
e_D = 0.001 #current eccentricity deimos
a_D = 6.92*R_Planet #current semi major axis deimos in m
AU=1.496e+11 #1 AU in meters
day = 60*60*24
M_Planet=M_Mars*1e-3 #MPlanet in kg
mars_density=3200 #density mars mantle in kg/m3
impactor_density=1000  #comet density in kg/m3
fimpactor = 0.4 #fraction of impactor material in the disk
fmars = 0.6     #fraction of mars material in the disk
particle_density=3000 #julien salmon 0.66*mars_density+0.33*impactor_density



#Mdisk in kg #from Salmon personal communication
if disk_type == 'low':
    ### low mass disk
    M_disk=(5e-7)*M_Planet

elif disk_type == 'medium':
    M_disk=(1e-5)*M_Planet

## High Mass Disk
elif disk_type == 'high':
    M_disk=(0.25e-4)*M_Planet



######################################################Rebound setup ####################################################
sim=rebound.Simulation()
sim.units = ('s','m','kg')
sim.G=6.67428e-11
sim.add(m=M_Planet,r=R_Planet)
########################################################################################################################
#Make folder to put all files created
NewDirName = output_path + "/"+disk_type+"_"+ ecc_type+"_"+str(n_particles)+"_"+seed  #store in output_path/disk_type/eccentricity_type/#particles/#seed
print("Making directory :", NewDirName)
if not os.path.exists(NewDirName):
    os.makedirs(NewDirName)

########################################################################################################################
# functions needed to create disk
def powerlaw(slope, min_v, max_v):
    y = np.random.uniform()
    pow_max = pow(max_v, slope+1.)
    pow_min = pow(min_v, slope+1.)
    return pow((pow_max-pow_min)*y + pow_min, 1./(slope+1.))

def rndm(a,b,g,size=1):
    #power law generation for pdf (x)\ x^{g-1} for a <=x<=b
    r =np.random.random(size=size)
    ag, bg = a**g, b**g
    return (ag + (bg -ag)*r)**(1./g)

def pickdensity(a,b,tracking,c):
    if tracking == 'no':
        return c
    if tracking == 'yes':
        test = np.random.random() # random float 0.0 <= x < 1.0
        if test < fimpactor:
            return a
        else:
            return b
def density(a,b):
    rho = a / ((4 / 3) * np.pi * np.power(np.divide(b, 2), 3))
    return rho
####################################################   CREATE DISK  ################################################
d=np.zeros((n_particles, 1))
m=d
total=0

### low mass disk
if disk_type == 'low':
    dmax=((3*3.99e-9*M_Planet/4/np.pi/particle_density)**(1/3))*2

    if n_particles == 1500:
        dmin = ((3*1.79e-10*M_Planet/4/np.pi/particle_density)**(1/3))*2
        print(dmax,dmin)
        while abs((total - M_disk) / M_disk) > 0.3:
            iexp = random.randrange(-50, -30, 5)/10
            d = rndm(dmin, dmax, iexp + 1,
                     n_particles)  # normally -3 + 1 due to power law d ~ -3 , but had to change to -5 + 1
            darray= [pickdensity(impactor_density,mars_density,tracking,particle_density) for i in range(0,n_particles)]
            m = np.multiply((4 / 3) * np.pi * np.power(np.divide(d, 2), 3),darray)
            total_particles = len(d)
            total = sum(m)
        print('disk created')


    elif n_particles == 5000:
        dmin = ((3*1.79e-10*M_Planet/4/np.pi/particle_density)**(1/3))*2 /3
        while abs((total - M_disk) / M_disk) > 0.3:
            d = rndm(dmin, dmax, -3 + 1,
                     n_particles)  # normally -3 + 1 due to power law d ~ -3 , but had to change to -5 + 1
            darray= [pickdensity(impactor_density,mars_density,tracking,particle_density) for i in range(0,n_particles)]
            m = np.multiply((4 / 3) * np.pi * np.power(np.divide(d, 2), 3),darray)
            total_particles = len(d)
            total = sum(m)
            # print(total/M_disk)
        print('disk created')


    elif n_particles == 15000:
        dmin = ((3 * 1.79e-10 * M_Planet / 4 / np.pi / particle_density) ** (1 / 3)) * 2 / 5

    while abs((total - M_disk) / M_disk) > 0.3:
        d = rndm(dmin, dmax, -3 + 1,
                 n_particles)  # normally -3 + 1 due to power law d ~ -3 , but had to change to -5 + 1
        darray = [pickdensity(impactor_density, mars_density, tracking, particle_density) for i in
                  range(0, n_particles)]
        m = np.multiply((4 / 3) * np.pi * np.power(np.divide(d, 2), 3), darray)
        total_particles = len(d)
        total = sum(m)
        #print(total/M_disk)
        print('disk created')


### Medium Mass disk
elif disk_type == 'medium':

    dmax=((6*1.0e-8*M_Planet/4/np.pi/particle_density)**(1/3))*2
    if n_particles ==1500:
        dmin = ((5*5.0e-10*M_Planet/4/np.pi/particle_density)**(1/3))*2
        #print(dmin, dmax, "dmin, dmax")
        while abs((total - M_disk) / M_disk) > 0.3:
            iexp = random.randrange(11, 31, 5)/10.0
            d = rndm(dmin, dmax, 1 - iexp, n_particles)  # normally -3 + 1 due to power law d ~ -3 , but had to change to 1 + 1
            darray= [pickdensity(impactor_density,mars_density,tracking,particle_density) for i in range(0,n_particles)]
            m = np.multiply((4 / 3) * np.pi * np.power(np.divide(d, 2), 3),darray)
            total_particles = len(d)
            total = sum(m)
            print((total-M_disk)/M_disk)
        print('disk created')


    elif n_particles == 5000:
        dmin = ((3*5.0e-10*M_Planet/4/np.pi/particle_density)**(1/3))*2 /3

        while abs((total - M_disk) / M_disk) > 0.3:
            d = rndm(dmin, dmax, -0.5 + 1,
                     n_particles)  # normally -3 + 1 due to power law d ~ -3 , but had to change to -0.5 + 1
            darray= [pickdensity(impactor_density,mars_density,tracking,particle_density) for i in range(0,n_particles)]
            m = np.multiply((4 / 3) * np.pi * np.power(np.divide(d, 2), 3),darray)
            total_particles = len(d)
            total = sum(m)

        print('disk created')


    elif n_particles == 15000:
        dmin = ((3*5.0e-10*M_Planet/4/np.pi/particle_density)**(1/3))*2 / 5

        while abs((total - M_disk) / M_disk) > 0.3:
            d = rndm(dmin, dmax, -1.5 + 1,
                     n_particles)  # normally -3 + 1 due to power law d ~ -3 , but had to change to -1.5 + 1
            darray= [pickdensity(impactor_density,mars_density,tracking,particle_density) for i in range(0,n_particles)]
            m = np.multiply((4 / 3) * np.pi * np.power(np.divide(d, 2), 3),darray)
            total_particles = len(d)
            total = sum(m)
            #print((total) / M_disk)
        print('disk created')

## High Mass Disk
elif disk_type == 'high':

    dmax=((3*4.4e-8*M_Planet/4/np.pi/particle_density)**(1/3))*2
    if n_particles ==1500:
        dmin = ((3*1.0e-9*M_Planet/4/np.pi/particle_density)**(1/3))*2

        while abs((total - M_disk) / M_disk) > 0.3:
            iexp = random.randrange(1, 99, 1)/100
            d = rndm(dmin, dmax, -iexp + 1, n_particles)  # normally -3 + 1 due to power law d ~ -3 , but had to change to -1.1 + 1
            darray= [pickdensity(impactor_density,mars_density,tracking,particle_density) for i in range(0,n_particles)]
            m = np.multiply((4 / 3) * np.pi * np.power(np.divide(d, 2), 3),darray)
            total_particles = len(d)
            total = sum(m)
            #print((total-M_disk) / M_disk)
        print('disk created')

    elif n_particles == 5000:
        dmin = ((3*1.0e-9*M_Planet/4/np.pi/particle_density)**(1/3))*2 /3
        while abs((total - M_disk) / M_disk) > 0.3:
            d = rndm(dmin, dmax, -1.5 + 1,
                     n_particles)  # normally -3 + 1 due to power law d ~ -3 , but had to change to -1.5 + 1
            darray= [pickdensity(impactor_density,mars_density,tracking,particle_density) for i in range(0,n_particles)]
            m = np.multiply((4 / 3) * np.pi * np.power(np.divide(d, 2), 3),darray)
            total_particles = len(d)
            total = sum(m)
            #print((total) / M_disk)
        print('disk created')

    elif n_particles == 15000:
        dmin = ((3*1.0e-9*M_Planet/4/np.pi/particle_density)**(1/3))*2 / 5
        while abs((total - M_disk) / M_disk) > 0.3:
            d = rndm(dmin, dmax, -2 + 1,
                     n_particles)  # normally -3 + 1 due to power law d ~ -3 , but had to change to -2 + 1
            darray= [pickdensity(impactor_density,mars_density,tracking,particle_density) for i in range(0,n_particles)]
            m = np.multiply((4 / 3) * np.pi * np.power(np.divide(d, 2), 3),darray)
            total_particles = len(d)
            total = sum(m)
        print('disk created')
textfile = open(NewDirName+"/disk_characteristics.txt",'w')
textfile.write("Debris disk Mass in kg = %f \n" % total)
textfile.write("Debris disk Mass / Mass Phobos = %f \n" % float(float(total)/float(M_P)) )
print(total, "Debris disk Mass in kg", total/M_P, "debris disk mass / Mass of Phobos")
print(min(m), "min mass in Debris disk Mass in kg")
## Create random eccentricities, inclinations, and semi-major axis
if ecc_type == "high":
    sigma_e = 1e-2
elif ecc_type == "low":
    sigma_e = 1e-3
sigma_inc = sigma_e /2
r = rndm(2.6 * R_Planet, 7 * R_Planet, -3, n_particles)
e = sigma_e*(-2*np.log(np.random.random(n_particles)))**(1/2)
print("Minimum and maximum eccentricity", min(e), max(e))
textfile.write("minimum eccentricity = %f \n" % min(e))
textfile.write("maximum eccentricity = %f \n" % max(e))
print("Minimum and maximum position in R_Planet", min(r)/R_Planet, float(max(r)/float(R_Planet)))
textfile.write("minimum semimajor axis in R_Planet = %f \n" % float(min(r)/float(R_Planet)))
textfile.write("maximum semimajor axis in R_Planet = %f \n" % float(max(r)/float(R_Planet)))
inc = sigma_inc*(-2*np.log(np.random.random(n_particles)))**(1/2)
print("Minimum and maximum inclination in degrees", min(inc*180/np.pi), max(inc*180/np.pi))
textfile.write("minimum inclination in degrees = %f \n" % min(inc*180/np.pi))
textfile.write("maximum inclination in degrees = %f \n" % max(inc*180/np.pi))
## Create workbook object of particles diameter in debris disk
if create==1:
    wb = openpyxl.Workbook()
    sheet = wb.active
    sheet.title = 'Planet_diameters_list'
    for i in range(0, len(d)):
        sheet.cell(row=i + 2, column=1).value = d[i]
## Finally, save the file and give it a name
    wb.save(NewDirName+'/canup.xlsx')
m_impactor =[]
a_impactor=[]
m_mars=[]
a_mars=[]

#################################### SIMULATION START #######################################
for i in range(0, n_particles):
    sim.add(m=m[i], r=d[i]/2, a=r[i], e=e[i], inc=inc[i], Omega=np.random.uniform(-np.pi, np.pi),
            omega=np.random.uniform(-np.pi, np.pi), M=np.random.uniform(-np.pi, np.pi))
    #print(density(m[i],d[i]/2))
    if abs(density(m[i],d[i])-impactor_density)/impactor_density <= 0.1:
        m_impactor.append(m[i])
        a_impactor.append(r[i])

    elif abs(density(m[i],d[i])-mars_density)/mars_density <= 0.1:
        m_mars.append(m[i])
        a_mars.append(r[i])
#print(a_impactor)
sim.move_to_com()
E0 = sim.calculate_energy()


hist = [m/M_Planet, r/R_Planet]
np.save(NewDirName +'/hist.npy',hist)
mint = [((sim.particles[j].a**3)*4*(np.pi**2)/(G*M_Planet))**(0.5)/60/60 for j in range(1, sim.N)] #%period orbits in hours
am = [sim.particles[j].a for j in range(1, sim.N)]
mm = [sim.particles[j].m * G for j in range(1, sim.N)]


print("E0", E0)
print("minimum semi major axis in m", min(am))
print(min(mm), "MTINY minimum G*Mass body in disk")
print(max(mm), "MMAX minimum G*Mass body in disk")
print(min(mint), "minimum orbital period of particles in hours")
print("dt =", int(min(mint)*3600/30), "DT in seconds")
print(max(am)/R_Planet,"maximum semi major axis")
textfile.write("E0 = %f \n" % E0)
textfile.write("minimum semi major axis in m = %f \n" % min(am))
textfile.write("MTINY minimum G*Mass body in disk = %f \n" % min(mm))
textfile.write("MMAX minimum G*Mass body in disk = %f \n" % max(mm))
textfile.write("minimum orbital period P of particles in hours = %f \n" % min(mint))
textfile.write("DT in seconds = P/30  = %f \n" % int(min(mint)*3600/30))
textfile.write("maximum semimajor axis in R_Planet = %f \n" % float(max(am)/float(R_Planet)))
textfile.close()

e = [sim.particles[j].e for j in range(1, sim.N)]
x = [sim.particles[j].x for j in range(1, sim.N)]
y = [sim.particles[j].y for j in range(1, sim.N)]
z = [sim.particles[j].z for j in range(1, sim.N)]
vx = [sim.particles[j].vx for j in range(1, sim.N)]
vy = [sim.particles[j].vy for j in range(1, sim.N)]
vz = [sim.particles[j].vz for j in range(1, sim.N)]
m = [sim.particles[j].m * G for j in range(1, sim.N)]
r = [sim.particles[j].r for j in range(1, sim.N)]
Rhill = [sim.particles[j].a * (sim.particles[j].m / (3 * M_Planet))**(1/3) for j in
         range(1, sim.N)]
density = np.multiply(3/4/np.pi/G,np.divide(m,np.power(r,3)))

if create_canup ==1:
    print("create input file")
    print("minimum period in hours", min(mint))
    with open(NewDirName+ '/mars.in', 'w') as output:
        output.write("%s        ! Mars System in SI units\n" % (sim.N))
        output.write("1 %s \n" % ("{:10.8e}".format(G * M_Planet)))
        output.write(".0 .0 .0        ! x y z\n")
        output.write(".0 .0 .0        !vx vy vz\n")
        for i in range(0, (sim.N - 1)):
            output.write("%s %s %s       ! particle number mass Rhill\n" % (
            i + 2, "{:10.8e}".format(m[i]), "{:10.8e}".format(Rhill[i])))
            output.write("%s                    !particle radius in AU\n" % ("{:10.8e}".format(r[i])))
            output.write(
                "%s %s %s      ! x y z \n" % ("{:10.8e}".format(x[i]), "{:10.8e}".format(y[i]), "{:10.8e}".format(z[i])))
            output.write("%s %s %s       ! vx vy vz \n" % (
            "{:10.8e}".format(vx[i]), "{:10.8e}".format(vy[i]), "{:10.8e}".format(vz[i])))



## Start initial 2D edge Figure
fig= plt.figure()
#cm = sns.coolwarm_palette(dark=.3, light=.8, as_cmap=True)
plt.scatter(np.divide(am,R_Planet), e, s=np.multiply(1e-7,np.multiply(r,r)),c=density, alpha=0.9, cmap = 'coolwarm') #s=np.power(np.multiply(np.pi / 2e7, r), 2)
plt.xlabel('Semimajor axis (R_Planet)')
cbar = plt.colorbar()
plt.ylabel('Eccentricity')
plt.xlim(2, 8)
#plt.yscale('log')
plt.ylim(1.0e-6, 5e-3)
cbar.set_label('Particles density', rotation =270)
plt.grid(b=True, which='major', color='#666666', linestyle='-',alpha=0.5)
plt.minorticks_on()
plt.grid(b=True, which='minor', color='#999999', linestyle='--', alpha = 0.2)
plt.title("Eccentricity of particles in the simulation at t=%6.1f: %d" % (int(0/365/24/60/60), (sim.N)-1), y=1.08)
plt.savefig(NewDirName + "/"+'disk_ecc_distance%s.png' % (i + 1))

plt.close()

## Start initial 2D
fig= plt.figure()
#cm = sns.coolwarm_palette(dark=.3, light=.8, as_cmap=True)
plt.scatter(np.divide(x,R_Planet), np.divide(y,R_Planet), s=np.multiply(1e-7,np.multiply(r,r)),c=density, alpha=0.9, cmap = 'coolwarm') #s=np.power(np.multiply(np.pi / 2e7, r), 2)
plt.xlabel('X position (R_Planet)')
cbar = plt.colorbar()
plt.ylabel('Y position (R_Planet')
plt.xlim(-8, 8)
#plt.yscale('log')
plt.ylim(-8, 8)
cbar.set_label('Particles density', rotation =270)
plt.grid(b=True, which='major', color='#666666', linestyle='-',alpha=0.5)
plt.minorticks_on()
plt.grid(b=True, which='minor', color='#999999', linestyle='--', alpha = 0.2)
plt.title("Top-Down view of martian disk at t=%6.1f: %d" % (int(0/365/24/60/60), (sim.N)-1), y=1.08)
plt.savefig(NewDirName + "/"+'2Ddisk%s.png' % (i + 1))

plt.close()

## Start initial 2D
fig= plt.figure()
#cm = sns.coolwarm_palette(dark=.3, light=.8, as_cmap=True)
plt.scatter(np.divide(x,R_Planet), np.divide(z,R_Planet), s=np.multiply(1e-7,np.multiply(r,r)),c=density, alpha=0.9, cmap = 'coolwarm') #s=np.power(np.multiply(np.pi / 2e7, r), 2)
plt.xlabel('X position (R_Planet)')
cbar = plt.colorbar()
plt.ylabel('Y position (R_Planet')
plt.xlim(-8, 8)
#plt.yscale('log')
plt.ylim(-0.05, 0.05)
cbar.set_label('Particles density', rotation =270)
plt.grid(b=True, which='major', color='#666666', linestyle='-',alpha=0.5)
plt.minorticks_on()
plt.grid(b=True, which='minor', color='#999999', linestyle='--', alpha = 0.2)
plt.title("Side view of martian disk at t=%6.1f: %d" % (int(0/365/24/60/60), (sim.N)-1), y=1.08)
plt.savefig(NewDirName + "/"+'2Dedgedisk%s.png' % (i + 1))

plt.close()



## Start initial 3D
fig= plt.figure()
ax = plt.axes(projection ="3d")
sctt = ax.scatter3D(np.divide(x,R_Planet), np.divide(y,R_Planet),np.divide(z,R_Planet), s=np.multiply(1e-7,np.multiply(r,r)),c=density, alpha=0.9, cmap = 'coolwarm') #s=np.power(np.multiply(np.pi / 2e7, r), 2)
plt.xlabel('X position (R_Planet)')
plt.ylabel('Y position (R_Planet')
plt.xlim(-8, 8)
ax.view_init(20,-120)
cbar = plt.colorbar(sctt, ax=ax)
#plt.yscale('log')
plt.ylim(-8, 8)
ax.set_zlim3d(-1, 1)
cbar.set_label('Particles density', rotation =270)
plt.grid(b=True, which='major', color='#666666', linestyle='-',alpha=0.5)
plt.minorticks_on()
plt.grid(b=True, which='minor', color='#999999', linestyle='--', alpha = 0.2)
plt.title("Top-Down view of martian disk at t=%6.1f: %d" % (int(0/365/24/60/60), (sim.N)-1), y=1.08)
plt.savefig(NewDirName + "/"+'3Ddisk%s.png' % (i + 1))

plt.close()


bins = np.linspace(2, 7, 100)
plt.hist(np.divide(a_mars,R_Planet), bins, alpha=0.5, label='Mars mantle')
plt.hist(np.divide(a_impactor,R_Planet), bins, alpha=0.5, label='Impactor material')
plt.legend(loc='upper right')
plt.title("Particle density profile of the debris disk")
plt.savefig(NewDirName+"/"+"disk.png")
plt.close()


plot=sns.distplot(np.divide(a_mars,R_Planet), label="Martian material")
plot=sns.distplot(np.divide(a_impactor,R_Planet),label="Impactor material")
fig = plot.get_figure()
plt.title("Particle density profile of the debris disk")
plt.legend(loc='upper right')
plt.xlabel('Semimajor axis (R_Planet)')
plt.ylabel('Frequency')
fig.savefig(NewDirName+"/snshist.png")
plt.close()


data = np.column_stack((np.divide(am,R_Planet),np.divide(m,M_P*G)))
data = pd.DataFrame(data, columns=['Semi-major axis','Mass'])
fig= plt.figure()
#df = pd.DataFrame(data, columns=["Mass", "Semi-major axis"])
plt.scatter(np.divide(am,R_Planet), np.divide(m,M_P*G),c=density,s=np.multiply(1.0e-7,np.multiply(r,r)), cmap="coolwarm")
cbar=plt.colorbar()
#sns_plot.ax_joint.set_xticks(np.arange(2, 7.5+1, 0.5))
plt.xlim(2, 8)
plt.yscale('log')
plt.ylim(1.0e-3, 1e2)
plt.xlabel('Semi-major axis of particles in Mars radii')
plt.ylabel('Mass of particles as a factor of Mass of Phobos')
plt.title("Semi-major axis of particles in the disk in function of mass at t=%6.1f: %d" % (int(0/365/24/60/60), (sim.N)-1),y =1.08 )
cbar.set_label('Particles density [$kg/m^{3}$]', rotation =270, labelpad = +20)
plt.grid(b=True, which='major', color='#666666', linestyle='-',alpha=0.4)
plt.minorticks_on()
plt.grid(b=True, which='minor', color='#999999', linestyle='--', alpha = 0.2)
#sns_plot.ax_joint.set_yscale('log')
#sns_plot.ax_joint.set_yticks([1e-5, 1.0e1, 1, 1.1])
plt.savefig(NewDirName+"/"+"disk_mass_distance.png")

MU2KG=1
DU2M =1
TU2S = 1
lfrag = 'no'
rmin = 10
rmax = 33890000.
J2 = 0.0
J4 = 0.0
t_0=0.0
end_sim = 365*3600*24*100
deltaT = 800
iout = 100
lenergy ="yes"
mtiny = 10
sys.stdout = open(NewDirName+"/"+"param.in", "w")
print(f'!Parameter file for the SyMBA-RINGMOONS test')
print(f'!NPLMAX         {n_particles}')
print(f'!NTPMAX         {0}')
print(f'T0             {t_0} ')
print(f'TSTOP          {end_sim}')
print(f'DT             {deltaT}')
print(f'PL_IN          mars.in')
print(f'TP_IN          tp.in')
print(f'IN_TYPE        ASCII')
print(f'ISTEP_OUT      {iout:d}')
print(f'BIN_OUT        bin.dat')
print(f'OUT_TYPE       REAL8')
print(f'OUT_FORM       EL')
print(f'OUT_STAT       NEW')
print(f'ISTEP_DUMP     {iout:d}')
print(f'J2             {J2}')
print(f'J4             {J4}')
print(f'CHK_CLOSE      yes')
print(f'CHK_RMIN       {rmin}')
print(f'CHK_RMAX       {rmax}')
print(f'CHK_EJECT      {rmax}')
print(f'CHK_QMIN       {rmin}')
print(f'CHK_QMIN_COORD HELIO')
print(f'CHK_QMIN_RANGE {rmin} {rmax}')
print(f'ENC_OUT        enc.dat')
print(f'EXTRA_FORCE    no')
print(f'BIG_DISCARD    no')
print(f'RHILL_PRESENT  yes')
print(f'MTINY          {mtiny}')
print(f'ENERGY        {lenergy}')
print(f'FRAGMENTATION  {lfrag}')
print(f'MU2KG          {MU2KG}')
print(f'DU2M          {DU2M}')
print(f'TU2S           {TU2S}')







