# DMM_plot_mirrors_Erange
#   Plot position and coordinates of the mirrors of a BEATS DMM for a given energy RANGE
#
#   References
#   -------
#   ______________________________________________________
#
#   Author:         Gianluca Iori (gianthk.iori@gmail.com)
#   SESAME - BEATS
#   Created on:   10/04/2020
#   Last update:  07/06/2020
#   ______________________________________________________

# import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import scipy as sc
font = {'weight' : 'normal',
        'size'   : 24}

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

def min_offset(length, theta, beam_height):
    # the OFFSET is defined by grazing @ the lowest possible energy
    offset_0 = length*np.sin(theta) # [m] OFFSET for zero mirror overlap
    # the OFFSET can be reduced in this way if the beam height is less than the acceptance of the first mirror
    if length*np.sin(theta) > beam_height:
        return 0.5*(offset_0 + beam_height) # [m] reduced OFFSET for zero shadow
    else:
        # print("Offset remains O_0")
        return offset_0

def M1_coordinates(l, theta):
    # coordinates of the extremes and of the midpoint of DMM Mirror 1
    #
    # we are looking at the beam from the side: coordinates are expressed as [z, y]
    # each mirror matrix (e.g. M1) contains the [z, y] coordinates of the three points:
    #   a: left edge
    #   o: midpoint
    #   b: right edge
    return np.array([[-(l / 2) * np.cos(theta), -(l / 2) * np.sin(theta)],
                     [0, 0],
                     [(l / 2) * np.cos(theta), (l / 2) * np.sin(theta)]])

# variables definition      ##################################################################################
fixed_exit = False # flag for fixed exit solution
l = 0.38 # [m] mirror length
E_1 = 6e3 # [eV] lowest foreseen energy
Y_1 = 5.5e-3 # [m] beam height @ first DMM mirror @ E_1
E_in = float(input("Insert start Energy [eV]: ")) # [eV] working energy
if E_in < E_1:
    print(f"warning! This is less than the Min En. Changed to Min. En.: {E_1} [eV]")
    E_in = E_1

E_end = float(input("Insert final Energy [eV]: ")) # [eV] working energy
if E_end < E_in + 4001:
    E_end = E_in + 4001

E = np.arange(E_in, E_end, ((E_end-E_in)/8))
lambda_1 = 1239.842/E_1 # nm
lambda_E = 1239.842/E # nm
d_spacing = float(input("Insert d-spacing [nm]: ")) # [nm]
# d_spacing = 4.0 # [nm]
d_spacing_0 = 4.0 # [nm] d-spacing of the stripe with thicker layers

o_0 = float(input("Insert target OFFSET [mm] (leave empty for MIN. OFFSET): "))*1e-3 # [mm]
# if o_0.isnumeric():
#     o_0 = float(o_0)*1e-3
# else:
#     o_0 = 0.0

theta_1 = np.arcsin(lambda_1/(2*d_spacing_0)) # [rad] grazing angle at lowest foreseen energy
theta_array = np.arcsin(lambda_E/(2*d_spacing)) # [rad] working grazing angle

# OFFSET                    ##################################################################################
if fixed_exit:
    # check if the OFFSET given by user is < than Min. affordable one
    if min_offset(l, theta_1, Y_1) > o_0:
        o_0_p = min_offset(l, theta_1, Y_1)
    else:
        o_0_p = o_0
else:
    # only check if input OFFSET was left empty and if @ Min. En. the second mirror projects shadow on the first one!
    o_min = min_offset(l, theta_1, Y_1)

    if o_0 == 0.0:
        o_0_p = o_min
    else:
        o_0_p = o_0
    if (l / 2) * np.sin(theta_array[0]) > o_0_p - (l / 2) * np.sin(theta_array[0]):
        print(f"warning! Given OFFSET is too small. Changed to Min. OFFSET: {o_min*1e3} [mm]")
        o_0_p = o_min

print(f"Min. Energy: {(E_1):.3} [eV]")

f, ax = plt.subplots()
exit_ray_x = o_0_p/np.tan(2*theta_array[-1]) + 0.5
for theta in theta_array:
    # DMM mirrors coordinates   ##################################################################################
    # we are looking at the beam from the side: coordinates are expressed as [z, y]
    # each mirror matrix (e.g. M1) contains the [z, y] coordinates of the three points: [left edge, midpoint, right edge]

    # FIRST MIRROR
    M1 = M1_coordinates(l, theta)

    # HORIZONTAL DISTANCE BETWEEN THE MIRRORS
    L_1_p = o_0_p/np.tan(2*theta) # [m]

    # SECOND MIRROR
    # M2 = np.array([[L_1_p-(l/2)*np.cos(theta), o_0_p-(l/2)*np.sin(theta)],
    #                [L_1_p, o_0_p],
    #                [L_1_p+(l/2)*np.cos(theta), o_0_p+(l/2)*np.sin(theta)]])
    M2 = M1 + L_1_p*(np.array([[1,0],[1,0],[1,0]], dtype='float')) + o_0_p*(np.array([[0,1],[0,1],[0,1]], dtype='float'))

    # TRACE EXTREME RAYS        ##################################################################################
    footprint = Y_1/np.sin(theta) # [m] footprint height on the mirrors
    if footprint > l:
        footprint = l

    # Extreme rays on FIRST MIRROR
    M1_ray = np.array([[-(footprint/2)*np.cos(theta), -(footprint/2)*np.sin(theta)],
                       [(footprint/2)*np.cos(theta), (footprint/2)*np.sin(theta)]])
    # Extreme rays on SECOND MIRROR
    M2_ray = np.array([[L_1_p-(footprint/2)*np.cos(theta), o_0_p-(footprint/2)*np.sin(theta)],
                       [L_1_p+(footprint/2)*np.cos(theta), o_0_p+(footprint/2)*np.sin(theta)]])

    ray_low = np.array([[-0.5, M1_ray[0,0], M2_ray[0,0], exit_ray_x],
                        [M1_ray[0,1], M1_ray[0,1], M2_ray[0,1], M2_ray[0,1]]])
    ray_high = np.array([[-0.5, M1_ray[1,0], M2_ray[1,0], exit_ray_x],
                        [M1_ray[1,1], M1_ray[1,1], M2_ray[1,1], M2_ray[1,1]]])

    # PLOT DMM MIRRORS          ##################################################################################
    plt.plot(M1[:,0], M1[:,1], 'k', linewidth=2.5)
    plt.plot(M2[:,0], M2[:,1], 'k', linewidth=2.5)
    plt.plot(ray_low[0,:], ray_low[1,:], 'c', linewidth=0.5)
    plt.plot(ray_high[0,:], ray_high[1,:], 'c', linewidth=0.5)
    plt.plot(M1[1, 0], M1[1, 1], 'Xr',
             M2[1,0], M2[1,1], 'Xr')

plt.xlabel('Z [m]')
plt.ylabel('Y [m]')
plt.grid(True, which="both")
# ax.set(xlim=(-0.28, 1.28), ylim=(-0.003, 0.019))

# ax.set_aspect('equal', 'box')
plt.show()
# f.savefig(f"DMM_mirrors_di{int(d_spacing_0)}_Ei{int(E_1/1e3)}_d{int(d_spacing)}_E{int(E/1e3)}_OFFSET{(o_0_p*1e3):.3}mm_Gr{rad2deg(theta):.3}deg_D{L_1_p:.3}m.png", bbox_inches='tight', dpi=600)

# PRINT OUTPUT              ##################################################################################
print(f"Energy: {E} [eV]")
print(f"d-spacing: {d_spacing:.3} [nm]")
print(f"Grazing angle: {rad2deg(theta):.3} [deg] ({(1e3*theta):.3} [mrad])")
print(f"DMM OFFSET: {(o_0_p*1e3):.3} [mm]")
print(f"Fixed exit: {fixed_exit}")
print(f"Min. mirrors distance: {o_0_p/np.tan(2*theta_array[0]):.3} [m] (E: {E[0]} [eV]; theta: {rad2deg(theta_array[0]):.5} [deg])")
print(f"Max. mirrors distance: {o_0_p/np.tan(2*theta_array[-1]):.3} [m] (E: {E[-1]} [eV]; theta: {rad2deg(theta_array[-1]):.5} [deg])")
