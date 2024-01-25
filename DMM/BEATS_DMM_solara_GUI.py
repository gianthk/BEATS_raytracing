# BEATS DMM solara GUI
#   Small solara GUI to plot position and coordinates of the mirrors of a BEATS DMM for a given energy
#
#   ______________________________________________________
#
#   Author:         Gianluca Iori (gianthk.iori@gmail.com)
#   SESAME - BEATS
#   Created on:   05/01/2024
#   Last update:  05/01/2024
#   ______________________________________________________

import pandas as pd
import solara
import solara.website
from pathlib import Path
import numpy as np
# import plotly.express as px
# import solara.express as spx
from matplotlib.figure import Figure

stripes = ["STRIPE 1 - [Ru/B4C]65 - d-spacing: 4.0 nm", "STRIPE 2 - [W/B4C]100 - d-spacing: 2.5 nm"]
fixed_exit = False # flag for fixed exit solution
l = 0.5 # [m] mirror length
L_1_p = 510*1e-3 # [m] horizontal distance between the mirrors
E_1 = 8e3 # [eV] lowest foreseen energy
Y_1 = 11e-3 # [m] beam height @ first DMM mirror @ E_1

stripe = solara.reactive("STRIPE 2 - [W/B4C]100 - d-spacing: 2.5 nm")
# energy = solara.reactive(20000.)
energies = np.arange(7000., 60000., 10.)
d_spacing1 = solara.reactive(2.5)
d_spacing2 = solara.reactive(4.0)
n_bilayers_1 = solara.reactive(100)
n_bilayers_2 = solara.reactive(65)

M1 = solara.reactive(np.array([[0, 0, 0],[0, 0, 0]]))
M2 = solara.reactive(np.array([[0, 0, 0],[0, 0, 0]]))
ray_low = solara.reactive(np.array([[0, 0, 0],[0, 0, 0]]))
ray_high = solara.reactive(np.array([[0, 0, 0],[0, 0, 0]]))

# theta_E_plot_path = Path(solara.website.__file__).parent / "resources/SESAME_DMM_Refl_Th_E_exp.GIF"
theta_E_plot_path = "./resources/BEATS_DMM_angle_versus_energy_calibration/SESAME_DMM_Refl_Th_E_exp.GIF"

# path to Bragg look-up tables for BEATS DMM
resources_dir = 'resources/BEATS_DMM_angle_versus_energy_calibration'

# load Bragg look-up table data
file_path_stripe1 = resources_dir + '/BEATS_230619_RuB4C_E_th.txt'
file_path_stripe2 = resources_dir + '/BEATS_230626_WB4C_E_th.txt'

# Column names for the DataFrame
column_names = ['E', 'theta']

data_stripe1 = pd.read_csv(file_path_stripe1, names=column_names, sep=' ', dtype=float)
data_stripe2 = pd.read_csv(file_path_stripe2, names=column_names, sep=' ', dtype=float)

def deg2rad(angle_deg):
    return angle_deg*(np.pi)/180

def rad2deg(angle_rad):
    return angle_rad*180/np.pi

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

# @solara.component
# def plot_DMM_position():
#     # check energy VS stripe
#     if energy.value > 25000:
#         stripe.set("STRIPE 2 - [W/B4C]100 - d-spacing: 2.5 nm")
#
#     if energy.value < 18000:
#         stripe.set("STRIPE 1 - [Ru/B4C]65 - d-spacing: 4.0 nm")
#
#     if stripe.value == "STRIPE 2 - [W/B4C]100 - d-spacing: 2.5 nm":
#         data_frame = data_stripe2
#     elif stripe.value == "STRIPE 1 - [Ru/B4C]65 - d-spacing: 4.0 nm":
#         data_frame = data_stripe1
#
#     index = data_frame.E[data_frame.E == energy.value].index[0]
#     lambda_E = 1239.842 / energy.value # nm
#     theta = 1e-3*data_frame.theta[index] # rad
#
#     # FIRST MIRROR
#     M1.set(M1_coordinates(l, theta))
#
#     # OFFSET
#     o_0_p = L_1_p * np.tan(2 * theta)  # [m]
#
#     # SECOND MIRROR
#     # M2 = np.array([[L_1_p-(l/2)*np.cos(theta), o_0_p-(l/2)*np.sin(theta)],
#     #                [L_1_p, o_0_p],
#     #                [L_1_p+(l/2)*np.cos(theta), o_0_p+(l/2)*np.sin(theta)]])
#     # ML2 will have the same angle as ML1.
#     # Its center position is translated by the distance between ML towers along X and by the DMM offset along Z.
#     M2.set(M1.value + L_1_p * (np.array([[1, 0], [1, 0], [1, 0]], dtype='float')) + o_0_p * (np.array([[0, 1], [0, 1], [0, 1]], dtype='float')))
#
#     # TRACE EXTREME RAYS
#     footprint = Y_1 / np.sin(theta)  # [m] footprint height on the mirrors
#     if footprint > l:
#         footprint = l
#
#     # Extreme rays on FIRST MIRROR
#     M1_ray = np.array([[-(footprint / 2) * np.cos(theta), -(footprint / 2) * np.sin(theta)],
#                        [(footprint / 2) * np.cos(theta), (footprint / 2) * np.sin(theta)]])
#     # Extreme rays on SECOND MIRROR
#     M2_ray = np.array([[L_1_p - (footprint / 2) * np.cos(theta), o_0_p - (footprint / 2) * np.sin(theta)],
#                        [L_1_p + (footprint / 2) * np.cos(theta), o_0_p + (footprint / 2) * np.sin(theta)]])
#
#     ray_low.set(np.array([[-1, M1_ray[0, 0], M2_ray[0, 0], M2_ray[0, 0] + 1],
#                         [M1_ray[0, 1], M1_ray[0, 1], M2_ray[0, 1], M2_ray[0, 1]]]))
#     ray_high.set(np.array([[-1, M1_ray[1, 0], M2_ray[1, 0], M2_ray[0, 0] + 1],
#                          [M1_ray[1, 1], M1_ray[1, 1], M2_ray[1, 1], M2_ray[1, 1]]]))
#
#     # PRINT OUTPUT
#     print(f"Energy: {energy.value} [eV]")
#     # print(f"d-spacing: {d_spacing:.3} [nm]")
#     print(f"Grazing angle: {rad2deg(theta):.3} [deg] ({(1e3 * theta):.3} [mrad])")
#     print(f"DMM OFFSET: {(o_0_p * 1e3):.3} [mm]")
#     # print(f"Fixed exit: {fixed_exit}")
#     print(f"Mirrors distance: {L_1_p:.3} [m]")
#     print(f"Mirrors length: {l * 1e3:.3} [mm]")
#     print(ray_high.value)


@solara.component
def Page():
    energy, set_energy = solara.use_state(20000.0)

    # check energy VS stripe
    if energy > 25000:
        stripe.set("STRIPE 2 - [W/B4C]100 - d-spacing: 2.5 nm")

    if energy < 18000:
        stripe.set("STRIPE 1 - [Ru/B4C]65 - d-spacing: 4.0 nm")

    if stripe.value == "STRIPE 2 - [W/B4C]100 - d-spacing: 2.5 nm":
        data_frame = data_stripe2
    elif stripe.value == "STRIPE 1 - [Ru/B4C]65 - d-spacing: 4.0 nm":
        data_frame = data_stripe1

    index = data_frame.E[data_frame.E == energy].index[0]
    lambda_E = 1239.842 / energy # nm
    theta = 1e-3*data_frame.theta[index] # rad

    # FIRST MIRROR
    M1 = M1_coordinates(l, theta)

    # OFFSET
    o_0_p = L_1_p * np.tan(2 * theta)  # [m]

    # SECOND MIRROR
    # M2 = np.array([[L_1_p-(l/2)*np.cos(theta), o_0_p-(l/2)*np.sin(theta)],
    #                [L_1_p, o_0_p],
    #                [L_1_p+(l/2)*np.cos(theta), o_0_p+(l/2)*np.sin(theta)]])
    # ML2 will have the same angle as ML1.
    # Its center position is translated by the distance between ML towers along X and by the DMM offset along Z.
    M2 = M1 + L_1_p * (np.array([[1, 0], [1, 0], [1, 0]], dtype='float')) + o_0_p * (np.array([[0, 1], [0, 1], [0, 1]], dtype='float'))

    # TRACE EXTREME RAYS
    footprint = Y_1 / np.sin(theta)  # [m] footprint height on the mirrors
    if footprint > l:
        footprint = l

    # Extreme rays on FIRST MIRROR
    M1_ray = np.array([[-(footprint / 2) * np.cos(theta), -(footprint / 2) * np.sin(theta)],
                       [(footprint / 2) * np.cos(theta), (footprint / 2) * np.sin(theta)]])
    # Extreme rays on SECOND MIRROR
    M2_ray = np.array([[L_1_p - (footprint / 2) * np.cos(theta), o_0_p - (footprint / 2) * np.sin(theta)],
                       [L_1_p + (footprint / 2) * np.cos(theta), o_0_p + (footprint / 2) * np.sin(theta)]])

    ray_low = np.array([[-1, M1_ray[0, 0], M2_ray[0, 0], M2_ray[0, 0] + 1],
                        [M1_ray[0, 1], M1_ray[0, 1], M2_ray[0, 1], M2_ray[0, 1]]])
    ray_high = np.array([[-1, M1_ray[1, 0], M2_ray[1, 0], M2_ray[0, 0] + 1],
                         [M1_ray[1, 1], M1_ray[1, 1], M2_ray[1, 1], M2_ray[1, 1]]])

    fig = Figure()
    ax = fig.subplots()
    ax.plot(M1[:, 0], M1[:, 1], 'k', linewidth=3)
    ax.plot(M2[:, 0], M2[:, 1], 'k', linewidth=3)
    ax.plot(M1[1, 0], M1[1, 1], 'Xr')
    ax.plot(M2[1, 0], M2[1, 1], 'Xr')
    ax.plot(ray_low[0, :], ray_low[1, :], 'c', linewidth=1)
    ax.plot(ray_high[0, :], ray_high[1, :], 'c', linewidth=1)

    # ax.xlabel('Z [m]')
    # ax.ylabel('Y [m]')
    ax.grid(True, which="both")
    # ax.xticks(np.arange(-0.4, 1, 0.1))
    # ax.yticks(np.arange(-0.005, 0.026, 0.001))
    ax.set(xlim=(-0.5, 1), ylim=(-0.005, 0.026))

    with solara.Card("Load dataset", classes=["my-2"]):
        # solara.Title('BEATS DMM GUI')
        solara.Markdown("This is the home page")

        with solara.Columns():
            with solara.Card("Bragg-law for BEATS DMM"):
                solara.Image(theta_E_plot_path, width='600px')
            with solara.Card("BEATS DMM configuration", margin=0, classes=["my-2"]):
                solara.FigureMatplotlib(fig, dependencies=[energy])

        with solara.Card("Select settings"):
            with solara.Row():
                solara.Select("DMM stripe:", value=stripe, values=stripes)
                # solara.Button(label="PLot DMM position", on_click=plot_DMM_position(), style={"height": "40px", "width": "400px"})

            solara.SliderInt(label="Energy [eV]", value=energy, min=7000, max=60000, step=10, thumb_label='always', on_value=set_energy) # , on_value=plot_DMM_position()

        with solara.Columns():
            with solara.Card('Stripe 1'):
                solara.InputFloat(label='d-spacing [nm]', value=d_spacing1, disabled = True)
            with solara.Card('Stripe 2'):
                solara.InputFloat(label='d-spacing [nm]', value=d_spacing2, disabled = True)


@solara.component
def Layout(children):
    return solara.AppLayout(children=children)

Page()