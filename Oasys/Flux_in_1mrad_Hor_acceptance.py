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
# print("Max flux: %g at: %d [eV]"%(flux.max(),eV[flux.argmax()]))
#print("Integrated power [W]: %4.2f  " % (cumulated_power[-1]))

power=eV*spectral_power
# plot(eV,flux,xlog=True,ylog=True,xtitle="Photon energy [eV]", ytitle="Flux [ph/s/0.1%bw]", title="Photon flux vs Energy")
# plot(eV,power,xlog=True,ylog=True,xtitle="Photon energy [eV]", ytitle="Power[W]", title="Power emitted vs Energy")


# ##############################################################
# photon flux density through 1 x 1 mm2)
# print(dir(in_object_3))
#print(in_object_3.get_contents
rays_tot = in_object_2.get_number_of_rays(0)        # total rays
rays_good = in_object_2.get_number_of_rays(1)        # good rays

# IF delta_E @ source and sample are equal:
flux_at_sample = flux*(rays_good/rays_tot)
# plot(eV,flux_at_sample,xlog=True,ylog=True,xtitle="Photon energy [eV]", ytitle="Spectral flux density [ph/s/0.1%bw]", title="Vacuum")

# PLOT FLUX DENSITY @ SAMPLE in 1 mrad Hor acceptance angle
# FWHM_X = 6360.7716e-3    # HOW THE HELL DO I GET FWHM_X and FWHM_Z FOM SHADOW?
# FWHM_Z = 3455.5716e-3
# window_sample = FWHM_X*FWHM_Z
plot(eV,flux_at_sample,xlog=True,ylog=True, xtitle="Photon energy [eV]", ytitle="Spectral flux density [ph/s/mrad/0.1%bw]")

# plt.plot(eV,flux_1mm, marker='', color='black', linewidth=2, linestyle='', label="Vacuum")
plt.grid(True, which="both")