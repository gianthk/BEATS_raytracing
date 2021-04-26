# DMM_plots
#   Plot acceptance of DMM for different d-spacings and energies
#
#   References
#   -------
#   ______________________________________________________
#
#   Author:         Gianluca Iori (gianthk.iori@gmail.com)
#   SESAME - BEATS
#   Created on:   13/04/2020
#   Last update:  07/06/2020
#   ______________________________________________________

# import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import scipy as sc
font = {'weight' : 'normal',
        'size'   : 14}

matplotlib.rc('font', **font)

def lin_interp(x, y, i, half):
    return x[i] + (x[i + 1] - x[i]) * ((half - y[i]) / (y[i + 1] - y[i]))


def fwhm(x, y):
    #   FULL WIDTH HALF MAXIMUM
    x = x.tolist()
    y = y.tolist()
    half = max(y) / 2.0
    signs = np.sign(np.add(y, -half))
    zero_crossings = (signs[0:-2] != signs[1:-1])
    zero_crossings_i = np.where(zero_crossings)[0]
    return [lin_interp(x, y, zero_crossings_i[1], half) - lin_interp(x, y, zero_crossings_i[0], half)]


def fwpm(x, y, p):
    #   FULL WIDTH AT GIVEN PERCENT OF MAXIMUM
    x = x.tolist()
    y = y.tolist()
    per = max(y) * p
    signs = np.sign(np.add(y, -per))
    zero_crossings = (signs[0:-2] != signs[1:-1])
    zero_crossings_i = np.where(zero_crossings)[0]
    return [lin_interp(x, y, zero_crossings_i[1], per) - lin_interp(x, y, zero_crossings_i[0], per)]


def deg2rad(angle_deg):
    return angle_deg*(np.pi)/180

def rad2deg(angle_rad):
    return angle_rad*180/np.pi

# variables definition      ##################################################################################
d_MM1 = 15.165 # [m] 16.155 distance of the first mirror from the source
d_sample = 43 # [m] distance of the sample from the source
d_DMM1_sample = d_sample-d_MM1 # [m] distance of the sample from the DMM

l = 0.5 # [m] mirror length
E_1 = 8e3 # [eV] lowest foreseen energy
Y_1 = 6e-3 # [m] beam height @ first DMM mirror @ E_1
E = np.arange(5e3, 65e3, 5e2) # [eV] array of the energies
lambda_1 = 1239.842/E_1 # nm
lambda_E = 1239.842/E # nm
d_spacing = np.array([3.0, 4.0]) # [nm]
d_spacing_0 = 4.0 # [nm] d-spacing of the stripe with thicker layers

theta_1 = np.arcsin(lambda_1/(2*d_spacing_0)) # [rad] grazing angle at lowest foreseen energy
theta = np.arcsin(np.transpose(np.array([lambda_E, lambda_E]))/(2*d_spacing)) # [rad] working grazing angle

# OFFSET                    ##################################################################################
# the OFFSET is defined by grazing @ the lowest possible energy
o_0 = l*np.sin(theta_1) # [m] OFFSET for zero mirror overlap
# the OFFSET can be reduced in this way if the beam height is less than the acceptance of the first mirror
if l*np.sin(theta_1) > Y_1:
    o_0_p = 0.5*(o_0 + Y_1) # [m] reduced OFFSET for zero shadow
else:
    # print("Offset remains O_0")
    o_0_p = o_0

print(f"Min. Energy: {(E_1):.3} [eV]")
print(f"Max. beam height @ entrance: {(1e3*Y_1):.3} [mm]")

# PLOT grazing angle and acceptance @ different energies and for different d-spacings          ###############
E2 = np.array([5, 8, 10, 15, 20, 25, 30, 35, 40, 45, 50])*1e3 # [eV]
# sigmap_Y = np.array([128.3, 117.5, 92.9, 74.1, 65.7, 59.9])*1e-6 # [rad]
FWHM_Y_MM1 = 1e-6*np.array([7798.1898, 6334.7404, 5572.6568, 4427.8019, 3884.0863, 3410.4061, 3076.2361, 2823.0118, 2647.6735, 2539.5069, 2358.9008]) # [m]

f, ax = plt.subplots()
plt.plot(E/1e3, (rad2deg(theta[:,0])), 'k')
plt.plot(E/1e3, (rad2deg(theta[:,1])), 'r')
# plt.plot(E/1e3, (rad2deg(theta[:,2])), 'b')
leg = [f"d = {d_spacing[0]} nm", f"d = {d_spacing[1]} nm"]
# plt.legend(['d-spacing = 3nm', 'd-spacing = 4nm'], loc='upper right') # 'd-spacing = 2.5nm',
plt.legend(leg, loc='upper right') # 'd-spacing = 2.5nm',
plt.xticks(np.arange(5, 55, 5))
plt.yticks(np.arange(0.15, 2.7, 0.1))
ax.set(xlim=(8, 50), ylim=(0.12, 1.55))
plt.xlabel('Photon energy [keV]')
plt.ylabel('grazing angle [deg]')
plt.grid(True, which="both")
plt.show()
# f.savefig(f"figures/DMM_theta_{d_spacing[0]}nm_{d_spacing[1]}nm.png", bbox_inches='tight', dpi=600)

# f, ax = plt.subplots()
# plt.plot(E/1e3, (1e3*l*np.sin(theta[:,0])), 'k')
# plt.plot(E/1e3, (1e3*l*np.sin(theta[:,1])), 'r')
# # plt.plot(E/1e3, (1e3*l*np.sin(theta[:,2])), 'b')
# plt.plot(E2/1e3, 1e3*FWHM_Y_MM1, '--m', linewidth=2)
# # plt.legend(['d-spacing = 2.5nm', 'd-spacing = 3nm', 'd-spacing = 4nm', 'beam height (FWHM)'], loc='upper right')
# plt.xticks(np.arange(8, 62, 2))
# plt.yticks(np.arange(1, 15, 1))
# ax.set(xlim=(8, 50), ylim=(0.7, 9))
# plt.xlabel('Photon energy [keV]')
# plt.ylabel('DMM mirror1 height [mm] \n (mirror length = 300mm)')
# plt.grid(True, which="both")
# plt.show()
# f.savefig(f"figures/DMM_mirror1_height_L{l*1e3}mm.png", bbox_inches='tight', dpi=600)

#
# # zoom plot
# f, ax = plt.subplots()
# plt.plot(E/1e3, (1e3*l*np.sin(theta[:,0])), 'k')
# plt.plot(E/1e3, (1e3*l*np.sin(theta[:,1])), '--k')
# plt.plot(E/1e3, (1e3*l*np.sin(theta[:,2])), '-.k')
# plt.xticks(np.arange(24, 60, 2))
# plt.yticks(np.arange(0.7, 4, 0.1))
# ax.set(xlim=(25, 60), ylim=(0.7, 4))
# plt.xlabel('Photon energy [keV]')
# plt.ylabel('mirror height [mm]')
# plt.grid(True, which="both")
# plt.legend(['d-spacing = 2nm', 'd-spacing = 3nm', 'd-spacing = 4nm'], loc='upper right')
# plt.show()

# PLOT beam size @ sample for different energies and d-spacings          ###############
f, ax = plt.subplots()
plt.plot(E/1e3, 1e3*(((l*np.sin(theta[:,0]))/d_MM1)*d_sample), 'k')
plt.plot(E/1e3, 1e3*(((l*np.sin(theta[:,1]))/d_MM1)*d_sample), 'r')
# plt.plot(E/1e3, 1e3*(((l*np.sin(theta[:,2]))/d_MM1)*d_sample), 'b')
plt.plot(E2/1e3, (1e3*FWHM_Y_MM1/d_MM1)*d_sample, '--m', linewidth=2)
plt.xlabel('Photon energy [keV]')
plt.ylabel(f'beam height @ {d_sample}m [mm] \n (mirror length = {l*1e3}mm)')
plt.grid(True, which="both")
plt.yticks(np.arange(0, 24, 2))
# beam size and divergences from Shadow at different energies          #################
plt.plot(E2/1e3, 1e3*(FWHM_Y_MM1/d_MM1)*d_sample, '--m', linewidth=2)
leg = [f"d = {d_spacing[0]} nm", f"d = {d_spacing[1]} nm", "WB available (FWHM)"]
plt.legend(leg, loc='upper right')
plt.xticks(np.arange(5, 55, 5))
plt.yticks(np.arange(1, 19, 1))
ax.set(xlim=(8, 50), ylim=(4, 19))
plt.show()
f.savefig(f"figures/DMM_beam_height_43m_L{l*1e3}mm.png", bbox_inches='tight', dpi=600)

# f, ax = plt.subplots()
# plt.plot(E2/1e3, 2.355e3*sigmap_Y*d_MM1, 'k')
# plt.xlabel('Photon energy [keV]')
# plt.ylabel("$beam height @ DMM1 [mm]$ ")
# plt.grid(True, which="both")
# plt.show()

