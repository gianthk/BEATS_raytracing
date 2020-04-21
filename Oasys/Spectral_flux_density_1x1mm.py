import numpy as np
import scipy
from srxraylib.plot.gol import plot
import matplotlib.pyplot as plt

# ##############################################################
# print(dir(in_object_1))

# print(in_object_1.get_contents("xoppy_data").shape)

#Xoppy data format (in this case):

#column 0: energy in eV
#column 1: flux in photons/sec/0.1%bw
#column 2: spectral power (W/eV)
#column 3: cumulated power(W)

eV = in_object_1.get_contents("xoppy_data")[:,0]
# EEc=in_object_1.get_contents("xoppy_data")[:,2]
flux = in_object_1.get_contents("xoppy_data")[:,1]
spectral_power=in_object_1.get_contents("xoppy_data")[:,2]
cumulated_power = in_object_1.get_contents("xoppy_data")[:,3]

# ##############################################################
# Plot flux and power
print("Max flux: %g at: %d [eV]"%(flux.max(),eV[flux.argmax()]))
#print("Integrated power [W]: %4.2f  " % (cumulated_power[-1]))

power=eV*spectral_power
# plot(eV,flux,xlog=True,ylog=True,xtitle="Photon energy [eV]", ytitle="Flux [ph/s/0.1%bw]", title="Photon flux vs Energy")
# plot(eV,power,xlog=True,ylog=True,xtitle="Photon energy [eV]", ytitle="Power[W]", title="Power emitted vs Energy")


# ##############################################################
# read SHADOW output (n_rays, intensity..)
# print(dir(in_object_3))
rays_tot = in_object_3.get_number_of_rays(0)         # total rays
rays_good = in_object_3.get_number_of_rays(1)        # good rays
#print(dir(in_object_3._beam))
intensity = in_object_3._beam.intensity(1)            # beam intensity

# ##############################################################
# photon flux density @ specific eV through 1 x 1 mm2
# from Abdellatief et al. 2017}
# F_sample(E) = F_source_0.1%BW(E) * (dE_source/(0.001*E)) * (intensity/rays_tot)

E = 50e3                      # scanning energy [10 20]
E_source_0 = 49000            # 9850 19700 
E_source_1 = 51000            # 10150 20300

# print(eV[np.argmin(np.abs(eV-E))])
F_source = flux[np.argmin(np.abs(eV-E))]
print("Flux_source @ %d [eV]: %g [Ph/s/0.1%%BW]"%(E, F_source))

# source bandwidth
dE_source = E_source_1 - E_source_0
print("dE_source: %d [eV] (%g%% BW)"%(dE_source, 100*(dE_source/E)))

# F_sample = F_source * (dE_source/(0.001*E)) * (intensity/rays_tot) * (1/(10*2))        # the last term is due to the screen @ sample [10x2] larger than 1x1 mm2
F_sample = F_source * (dE_source/(0.001*E)) * (intensity/rays_tot)
print("Flux_sample @ %d [eV]: %g [Ph/s/DMM BW]"%(E, F_sample))
