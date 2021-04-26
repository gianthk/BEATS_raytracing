# DMM_plots_coatings
#   Plot reflectivity and dE/E of different DMM coatinga from XOP Oasys simulations
#
#   References
#   -------
#   ______________________________________________________
#
#   Author:         Gianluca Iori (gianthk.iori@gmail.com)
#   SESAME - BEATS
#   Created on:   26/04/2020
#   Last update:  03/09/2020
#   ______________________________________________________

import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('TkAgg')
# import pandas as pd
import numpy as np
import statistics as st
import scipy as sc
from scipy.constants import c, h
font = {'weight' : 'normal',
        'size'   : 14}

matplotlib.rc('font', **font)


def lin_interp(x, y, i, half):
    return x[i] + (x[i+1] - x[i]) * ((half - y[i]) / (y[i+1] - y[i]))

def moving_av(list, N):
    cumsum, moving_aves = [0], []
    for i, x in enumerate(list, 1):
        cumsum.append(cumsum[i-1] + x)
        if i>=N:
            moving_ave = (cumsum[i] - cumsum[i-N])/N
            moving_aves.append(moving_ave)
    return moving_aves

def fwhm(x, y):
    x = x.tolist()
    y = y.tolist()
    half = max(y)/2.0
    signs = np.sign(np.add(y, -half))
    zero_crossings = (signs[0:-2] != signs[1:-1])
    zero_crossings_i = np.where(zero_crossings)[0]
    if zero_crossings_i.size < 2:
        return np.array([0])
    else:
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
Oasys_folder = '/home/gianthk/Oasys/BEATS/BEATS_BL_Oasys/DMM_BEATS/CfT/'
d_spacing = np.array([2.5, 3.0])  # [nm] [2.5, 3, 4]
d_spacing_Mo = np.array([4.0, 4.0]) # [nm]
n_bilayers = np.array([80, 100]) # [100, 60, 30]

En_lim = np.array([20e3, 50e3]) # eV
theta_lim = np.array([[0.22, 0.75],
                      [0.22, 0.75]]) # deg

En_lim_Mo = np.array([7e3, 23e3]) # eV
theta_lim_Mo = np.array([[0.37, 1.3],
                        [0.37, 1.3]]) # deg
n_bilayers_Mo = np.array([60, 80, 100]) # [100, 60, 30]

# energy array for the X axis
array_length = np.array([201, 501])
en_step = np.round((En_lim[1]-En_lim[0])/array_length[0], decimals=1)
en_step_Mo = np.round((En_lim_Mo[1]-En_lim_Mo[0])/array_length[0], decimals=8)
en = np.arange(En_lim[0],En_lim[1],en_step)
en_Mo = np.arange(En_lim_Mo[0],En_lim_Mo[1],en_step_Mo)
en2 = np.array([19e3, 29e3]) # just 2 energy points to show the possibility of total reflection
# thetA2_array = np.array([[0.728, 0.485],
#                          [0.728, 0.485],
#                          [0.613, 0.409],
#                          [0.470, 0.313]])
thetA2_array = np.array([[0.37, 1.3],
                         [0.37, 1.3]])

l = (h*c)/(en*1.6e-19) # wavelength [m]
l2 = (h*c)/(en2*1.6e-19) # wavelength [m]
l_Mo = (h*c)/(en_Mo*1.6e-19) # wavelength [m]

# initialize arrays for reflectivity and BW
Refl_peak = np.zeros([d_spacing.shape[0], n_bilayers.shape[0], array_length[0]])
theta_peak_ind = np.zeros([d_spacing.shape[0], n_bilayers.shape[0], array_length[1]], dtype='int')
theta_peak = np.zeros([d_spacing.shape[0], n_bilayers.shape[0], array_length[1]])
BW_peak = np.zeros([d_spacing.shape[0], n_bilayers.shape[0], array_length[0]])
Refl_int = np.zeros([d_spacing.shape[0], n_bilayers.shape[0], array_length[0]])

Refl_peak_Mo = np.zeros([d_spacing_Mo.shape[0], n_bilayers_Mo.shape[0], array_length[0]])
theta_peak_ind_Mo = np.zeros([d_spacing_Mo.shape[0], n_bilayers_Mo.shape[0], array_length[1]], dtype='int')
theta_peak_Mo = np.zeros([d_spacing_Mo.shape[0], n_bilayers_Mo.shape[0], array_length[1]])
BW_peak_Mo = np.zeros([d_spacing_Mo.shape[0], n_bilayers_Mo.shape[0], array_length[0]])
Refl_int_Mo = np.zeros([d_spacing_Mo.shape[0], n_bilayers_Mo.shape[0], array_length[0]])

# initialize arrays for plots of total reflection
thetA2_ind = np.zeros([d_spacing.shape[0], en2.shape[0]], dtype='int')
Refl_tot = np.zeros([d_spacing.shape[0], array_length[0], en2.shape[0]])
thetA2_Mo_ind = np.zeros([d_spacing_Mo.shape[0], en2.shape[0]], dtype='int')
thetA2_Mo = np.zeros([d_spacing_Mo.shape[0], en2.shape[0]])
Refl_Mo_tot = np.zeros([d_spacing_Mo.shape[0], array_length[0], en2.shape[0]])

# STRIPE_1 #############################################################################################################
# for each ML and d_spacing investigated
count_b = 0
for n in  n_bilayers:
    count = 0
    for d in d_spacing:
        # read 2D reflectivity data
        Refl_data = np.load(Oasys_folder + f'DMM_SiWB4C_{n_bilayers[count_b]}_d{d}_Refl.npy')
        # f, ax = plt.subplots()
        # plt.imshow(Refl_data, extent=[En_lim[0], En_lim[1], theta_lim[count,1], theta_lim[count,0]], aspect='auto')

        # angle array for the X axis
        angle_step = np.round((theta_lim[count, 1] - theta_lim[count, 0]) / array_length[1], decimals=5)
        angle = np.arange(theta_lim[count, 0], theta_lim[count, 1] + (angle_step / 2), angle_step)

        # theta from Bragg's law
        d_sp = float(d*1e-9)
        thetA = rad2deg(np.arcsin(l/(2*d_sp))) # [deg] grazing angle
        thetA2 = thetA2_array[count,:]  # [deg] grazing angle

        # find index for thetA2 angles and extract plots of total reflection
        count2 = 0
        for angle2 in thetA2:
            thetA2_ind[count, count2] = np.argmin(np.abs(angle - angle2))
            Refl_tot[count, :, count2] = Refl_data[thetA2_ind[count, count2], :]
            count2 = count2 + 1

        # remove area of total reflection
        Refl_data2 = Refl_data
        for col in range(0, Refl_data.shape[1]):
            for row in range(0, Refl_data.shape[0]):
                if thetA[col] > angle[row] + 0.1:
                    Refl_data2[row, col] = 0

        # get peak reflectivity and grazing angle for each photon energy
        for col in range(0, Refl_data2.shape[1]):
            Refl_peak[count, count_b, col] = max(Refl_data2[:, col])**2
            theta_peak_ind[count, count_b, col] = np.argmax(Refl_data2[:, col], axis=0)
            theta_peak[count, count_b, col] = angle[theta_peak_ind[count, count_b, col]]

        # # find index for thetA2 angles and extract plots of total reflection
        # count2 = 0
        # for angle2 in thetA2:
        #     thetA2_ind[count, count2] = np.argmin(np.abs(theta_peak[count, :] - angle2))
        #     print(thetA2_ind)
        #     Refl_tot[count, :] = Refl_data[thetA2_ind[count, count2], :]
        #     count2 = count2 + 1

        # get BW for each photon energy
        for col in range(0, Refl_data2.shape[1]):
            if theta_peak_ind[count, count_b, col] > 0:
                BW_peak[count, count_b, col] = 100 * fwhm(en, Refl_data2[theta_peak_ind[count, count_b, col], :])[0] / en[col]
                # Refl_int[count, count_b, col] = Refl_peak[count, count_b, col] * BW_peak[count, count_b, col] / 100
        BW_peak[count, count_b, 10:-10] = moving_av(BW_peak[count, count_b, :], 21)
        Refl_int[count, count_b, :] = Refl_peak[count, count_b, :] * (BW_peak[count, count_b, :]) / 100

        count = count + 1
    count_b = count_b + 1

# STRIPE_2 #############################################################################################################
# for each ML and n_layers investigated
count_b = 0
for n in n_bilayers_Mo:
    count = 0
    for d in d_spacing_Mo:
        # read 2D reflectivity data
        if count == 0:
            Refl_data = np.load(Oasys_folder + f'DMM_SiMoB4C_{n_bilayers_Mo[count_b]}_d{d}_Refl.npy')
        else:
            Refl_data = np.load(Oasys_folder + f'DMM_SiRuB4C_{n_bilayers_Mo[count_b]}_d{d}_Refl.npy')

        # energy array for the X axis
        angle_step = np.round((theta_lim_Mo[count, 1] - theta_lim_Mo[count, 0]) / array_length[1], decimals=5)
        angle = np.arange(theta_lim_Mo[count, 0], theta_lim_Mo[count, 1]+0.001*angle_step, angle_step)

        # theta from Bragg's law
        d_sp = float(d*1e-9)
        thetA = rad2deg(np.arcsin(l_Mo/(2*d_sp))) # [deg] grazing angle
        thetA2 = rad2deg(np.arcsin(l2 / (2 * d_sp)))  # [deg] grazing angle

        # find index for thetA2 angles and extract plots of total reflection
        count2 = 0
        for angle2 in thetA2:
            thetA2_Mo_ind[count, count2] = np.argmin(np.abs(angle - angle2))
            thetA2_Mo[count, count2] = angle[thetA2_Mo_ind[count, count2]]
            Refl_Mo_tot[count, :, count2] = Refl_data[thetA2_Mo_ind[count, count2], :]
            count2 = count2 + 1

        # remove area of total reflection
        Refl_data2 = Refl_data
        for col in range(0, Refl_data.shape[1]):
            for row in range(0, Refl_data.shape[0]):
                if thetA[col] > angle[row] + 0.1:
                    Refl_data2[row, col] = 0

        # get peak reflectivity and grazing angle for each photon energy
        for col in range(0, Refl_data2.shape[1]):
            Refl_peak_Mo[count, count_b, col] = max(Refl_data2[:, col])**2
            theta_peak_ind_Mo[count, count_b, col] = np.argmax(Refl_data2[:, col], axis=0)
            theta_peak_Mo[count, count_b, col] = angle[theta_peak_ind_Mo[count, count_b, col]]

        # get BW for each photon energy
        for col in range(0, Refl_data2.shape[1]):
            if theta_peak_ind_Mo[count, count_b, col] > 0:
                BW_peak_Mo[count, count_b, col] = 100 * fwhm(en, Refl_data2[theta_peak_ind_Mo[count, count_b, col], :])[0] / en[col]
                # Refl_int_Mo[count, count_b, col] = Refl_peak_Mo[count, count_b, col] * BW_peak_Mo[count, col] / 100
        BW_peak_Mo[count, count_b, 10:-10] = moving_av(BW_peak_Mo[count, count_b, :], 21)
        Refl_int_Mo[count, count_b, :] = Refl_peak_Mo[count, count_b, :] * (BW_peak_Mo[count, count_b, :]) / 100

        count = count + 1
    count_b = count_b + 1

# # plots of the Peak Reflectivity
# f, ax = plt.subplots()
# plt.plot(en*1e-3, Refl_peak[0,0,:], 'k')
# plt.plot(en*1e-3, Refl_peak[0,1,:], '--k')
# plt.plot(en*1e-3, Refl_peak[1,0,:], 'r')
# plt.plot(en*1e-3, Refl_peak[1,1,:], '--r')
# ax.set(xlim=(19, 51), ylim=(0.4, 0.9))
# # plt.plot(en, Refl_peak_Mo[0,:]**2, 'r')
# plt.grid(True, which="both")
# plt.xlabel('Photon Energy [keV]')
# plt.ylabel('$Refl.^2 (peak)$')
# plt.legend(['$[W/B_4C]_{80}$ d2.5nm', '$[W/B_4C]_{100}$ d2.5nm', '$[W/B_4C]_{80}$ d3.0nm', '$[W/B_4C]_{100}$ d3.0nm'], loc='lower right')
# plt.show()
# f.savefig(f"CfT_BEATS_WB4C_DMM_ReflPEAK.png", bbox_inches='tight', dpi=600)
#
# # plots of the dE/E
# f, ax = plt.subplots()
# plt.plot(en*1e-3, BW_peak[0,0,:], 'k')
# plt.plot(en*1e-3, BW_peak[0,1,:], '--k')
# plt.plot(en*1e-3, BW_peak[1,0,:], 'r')
# plt.plot(en*1e-3, BW_peak[1,1,:], '--r')
# ax.set(xlim=(19, 51), ylim=(1.5, 3.5))
# # plt.plot(en, Refl_peak_Mo[0,:]**2, 'r')
# plt.grid(True, which="both")
# plt.xlabel('Photon Energy [keV]')
# plt.ylabel('$dE/E$')
# # plt.legend(['$[W/B_4C]_{80}$ d2.5nm', '$[W/B_4C]_{100}$ d2.5nm', '$[W/B_4C]_{80}$ d3.0nm', '$[W/B_4C]_{100}$ d3.0nm'], loc='lower right')
# plt.show()
# f.savefig(f"CfT_BEATS_WB4C_DMM_BW.png", bbox_inches='tight', dpi=600)
#
# # plots of the integrated Reflectivity
# f, ax = plt.subplots()
# plt.plot(en*1e-3, Refl_int[0,0,:], 'k')
# plt.plot(en*1e-3, Refl_int[0,1,:], '--k')
# plt.plot(en*1e-3, Refl_int[1,0,:], 'r')
# plt.plot(en*1e-3, Refl_int[1,1,:], '--r')
# ax.set(xlim=(19, 51), ylim=(0.01, 0.03))
# # plt.plot(en, Refl_peak_Mo[0,:]**2, 'r')
# plt.grid(True, which="both")
# plt.xlabel('Photon Energy [keV]')
# plt.ylabel('$Refl (int)$')
# # plt.legend(['$[W/B_4C]_{80}$ d2.5nm', '$[W/B_4C]_{100}$ d2.5nm', '$[W/B_4C]_{80}$ d3.0nm', '$[W/B_4C]_{100}$ d3.0nm'], loc='lower right')
# plt.show()
# f.savefig(f"CfT_BEATS_WB4C_DMM_ReflINT.png", bbox_inches='tight', dpi=600)
#
# # Mo or Ru - B4C PLOTS  ################################################################################################
# # plots of the Peak Reflectivity
# f, ax = plt.subplots()
# plt.plot(en_Mo*1e-3, Refl_peak_Mo[0,0,:], 'k')
# plt.plot(en_Mo*1e-3, Refl_peak_Mo[0,1,:], '--k')
# plt.plot(en_Mo*1e-3, Refl_peak_Mo[0,2,:], '-.k')
# plt.plot(en_Mo*1e-3, Refl_peak_Mo[1,0,:], 'r')
# plt.plot(en_Mo*1e-3, Refl_peak_Mo[1,1,:], '--r')
# plt.plot(en_Mo*1e-3, Refl_peak_Mo[1,2,:], '-.r')
# ax.set(xlim=(7, 24), ylim=(0.4, 0.9))
# # plt.plot(en, Refl_peak_Mo[0,:]**2, 'r')
# plt.grid(True, which="both")
# plt.xlabel('Photon Energy [keV]')
# plt.ylabel('$Refl.^2 (peak)$')
# # plt.legend(['$[Mo/B_4C]_{60}$ d4nm', '$[Mo/B_4C]_{80}$ d4nm', '$[Mo/B_4C]_{100}$ d4nm', '$[Ru/B_4C]_{60}$ d4nm', '$[Ru/B_4C]_{80}$ d4nm', '$[Ru/B_4C]_{100}$ d4nm'], loc='lower right')
# plt.show()
# f.savefig(f"CfT_BEATS_RuB4C_DMM_ReflPEAK.png", bbox_inches='tight', dpi=600)
#
# # plots of the dE/E
# f, ax = plt.subplots()
# plt.plot(en_Mo*1e-3, BW_peak_Mo[0,0,:], 'k')
# plt.plot(en_Mo*1e-3, BW_peak_Mo[0,1,:], '--k')
# plt.plot(en_Mo*1e-3, BW_peak_Mo[0,2,:], '-.k')
# plt.plot(en_Mo*1e-3, BW_peak_Mo[1,0,:], 'r')
# plt.plot(en_Mo*1e-3, BW_peak_Mo[1,1,:], '--r')
# plt.plot(en_Mo*1e-3, BW_peak_Mo[1,2,:], '-.r')
# ax.set(xlim=(7, 24), ylim=(0.0, 3.7))
# # plt.plot(en, Refl_peak_Mo[0,:]**2, 'r')
# plt.grid(True, which="both")
# plt.xlabel('Photon Energy [keV]')
# plt.ylabel('$dE/E$')
# plt.legend(['$[Mo/B_4C]_{60}$ d4nm', '$[Mo/B_4C]_{80}$ d4nm', '$[Mo/B_4C]_{100}$ d4nm', '$[Ru/B_4C]_{60}$ d4nm', '$[Ru/B_4C]_{80}$ d4nm', '$[Ru/B_4C]_{100}$ d4nm'], loc='lower right')
# plt.show()
# f.savefig(f"CfT_BEATS_RuB4C_DMM_BW.png", bbox_inches='tight', dpi=600)

# plots of the integrated Reflectivity
f, ax = plt.subplots()
plt.plot(en_Mo*1e-3, Refl_int_Mo[0,0,:], 'k')
plt.plot(en_Mo*1e-3, Refl_int_Mo[0,1,:], '--k')
plt.plot(en_Mo*1e-3, Refl_int_Mo[0,2,:], '-.k')
plt.plot(en_Mo*1e-3, Refl_int_Mo[1,0,:], 'r')
plt.plot(en_Mo*1e-3, Refl_int_Mo[1,1,:], '--r')
plt.plot(en_Mo*1e-3, Refl_int_Mo[1,2,:], '-.r')
ax.set(xlim=(7, 24), ylim=(0.008, 0.03))
# plt.plot(en, Refl_peak_Mo[0,:]**2, 'r')
plt.grid(True, which="both")
plt.xlabel('Photon Energy [keV]')
plt.ylabel('$Refl (int)$')
# plt.legend(['$[Mo/B_4C]_{60}$ d4nm', '$[Mo/B_4C]_{80}$ d4nm', '$[Mo/B_4C]_{100}$ d4nm', '$[Ru/B_4C]_{60}$ d4nm', '$[Ru/B_4C]_{80}$ d4nm', '$[Ru/B_4C]_{100}$ d4nm'], loc='lower right')
plt.show()
f.savefig(f"CfT_BEATS_RuB4C_DMM_ReflINT.png", bbox_inches='tight', dpi=600)
