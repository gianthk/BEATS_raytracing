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
import plotly.express as px
import solara.express as spx
# from matplotlib.figure import Figure

stripes = ["STRIPE 1 - [Ru/B4C]65 - d-spacing: 4.0 nm", "STRIPE 2 - [W/B4C]100 - d-spacing: 2.5 nm"]
stripe = solara.reactive("STRIPE 2 - [W/B4C]100 - d-spacing: 2.5 nm")
energy = solara.reactive(20000.)
energies = np.arange(7000., 60000., 10.)
d_spacing1 = solara.reactive(2.5)
d_spacing2 = solara.reactive(4.0)

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

@solara.component
def Page():
    with solara.Card("Load dataset", classes=["my-2"]):
        # solara.Title('BEATS DMM GUI')
        solara.Markdown("This is the home page")

        with solara.Columns():
            with solara.Card("Bragg-law for BEATS DMM"):
                solara.Image(theta_E_plot_path, width='600px')
            with solara.Card("BEATS DMM configuration", margin=0, classes=["my-2"]):
                solara.Image(theta_E_plot_path, width='600px')

        with solara.Card("Select settings"):
            solara.Select("DMM stripe:", value=stripe, values=stripes)
            solara.SliderInt(label="Energy [eV]", value=energy, min=7000, max=60000, step=10, thumb_label='always')

        with solara.Columns():
            with solara.Card('Stripe 1'):
                solara.InputFloat(label='d-spacing [nm]', value=d_spacing1, disabled = True)
            with solara.Card('Stripe 2'):
                solara.InputFloat(label='d-spacing [nm]', value=d_spacing2, disabled = True)


@solara.component
def Layout(children):
    return solara.AppLayout(children=children)

Page()