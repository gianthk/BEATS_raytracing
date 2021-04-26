#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Plot .H5 beam profiles.

For more information, call this script with the help option::
    h5plot.py -h

"""

"""
2DO:

"""

__author__ = 'Gianluca Iori'
__date_created__ = '2021-04-12'
__date__ = '2021-04-24'
__docformat__ = 'restructuredtext en'
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = 'Gianluca Iori'
__email__ = "gianthk.iori@gmail.com"


import logging
import numpy as np
from scipy.stats import norm
# import pandas as pd
import matplotlib.pyplot as plt
# from pandas_ods_reader import read_ods
import argparse
import textwrap
import silx.io
import os, glob
from astropy import modeling


def main():
    description = textwrap.dedent('''\
            Convert all .H5 beam profile plot files in given folder to .CSV data files.
            ''')

    epilog = textwrap.dedent('''\
            EXAMPLES:            
            * Convert given .h5 file to .csv:
                h5plot.py input_file.h5 -o output_file.csv
            
            * Convert all .h5 files in a given folder:
                h5plot.py /.../source_folder
            ''')

    parser = argparse.ArgumentParser(description = description, epilog=epilog, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('filein', type=str, default=None, help='Input filename (.h5).')
    parser.add_argument('-o', '--fileout', type=str, default=None, help='Output image filename.')
    parser.add_argument('-x', '--xsize', type=int, default=None, help='Crop size X (around center).')
    parser.add_argument('-y', '--ysize', type=int, default=None, help='Crop size Y (around center).')
    parser.add_argument('-w', '--avgwidth', type=int, default=None, help='Width for vertical Avg profile.')
    parser.add_argument('-v', '--verbose', type=bool, default=False, help='Verbose output')

    args = parser.parse_args()

    # verbose output
    if args.verbose is True:
        logging.basicConfig(level=logging.INFO)

    # check if input arg is file or folder
    # checks if path is a file
    isFile = os.path.isfile(args.filein)

    # checks if path is a directory
    isDirectory = os.path.isdir(args.filein)

    # get all .h5 files in dir
    if isDirectory:
        files = list(glob.glob(os.path.join(args.filein, '*.hdf5')))
    else:
        files = [args.filein]

    # for all files in list
    for file in files:

        if isDirectory:
            filename_out = os.path.splitext(file)[0] + ".csv"
            filename_out_png = os.path.splitext(file)[0] + ".png"
            filename_out_profile = os.path.splitext(file)[0] + "_Yprofile.png"
        elif args.fileout is None:
            filename_out = os.path.splitext(args.filein)[0] + ".csv"
            filename_out_png = os.path.splitext(args.filein)[0] + ".png"
            filename_out_profile = os.path.splitext(args.filein)[0] + "_Yprofile.png"
        else:
            filename_out = args.fileout
            filename_out_png = os.path.splitext(args.fileout)[0] + ".png"
            filename_out_profile = os.path.splitext(args.fileout)[0] + "_Yprofile.png"

        # read with silx.io
        sf = silx.io.open(file)
        grp = sf["xy_plots"]
        dataset = grp["last_plot"]
        image = np.transpose(dataset["23: Total Intensity = |Eσ|² + |Eπ|²"][:,:])
        
        coors = sf["coordinates"]

        if args.xsize is not None:
            X_size = args.xsize  # pixels for extracting CSV profile (around center)
            X_0 = image.shape[1] // 2 - X_size // 2
            X_1 = X_0 + X_size
            x = coors["X"][X_0:X_1].astype("float16")

            if args.ysize is not None:
                Y_size = args.ysize  # pixels for extracting CSV profile (around center)
                Y_0 = image.shape[0] // 2 - Y_size // 2
                Y_1 = Y_0 + Y_size

                image = image[Y_0:Y_1, X_0:X_1].astype("float16")
                # the first dim of image is Y, the second is X
                # image[477, 538] = 0.006057
                y = coors["Y"][Y_0:Y_1].astype("float16")

            else:
                Y_size = image.shape[0]
                image = image[:, X_0:X_1].astype("float16")
                y = coors["Y"][:].astype("float16")

        else:
            X_size = image.shape[1]
            x = coors["X"][:].astype("float16")
            if args.ysize is not None:
                Y_size = args.ysize  # pixels for extracting CSV profile (around center)
                Y_0 = image.shape[0] // 2 - Y_size // 2
                Y_1 = Y_0 + Y_size

                image = image[Y_0:Y_1, :].astype("float16")
                # the first dim of image is Y, the second is X
                # image[477, 538] = 0.006057
                y = coors["Y"][Y_0:Y_1].astype("float16")

            else:
                Y_size = image.shape[0]
                image = image[:, :].astype("float16")
                y = coors["Y"][:].astype("float16")

        # print the profile
        f, ax = plt.subplots()
        plt.imshow(image.astype("float"), extent=[x.min(),x.max(),y.min(),y.max()], aspect=4*(X_size/Y_size))
        plt.xlabel('X [mm]')
        plt.ylabel('Y [mm]')
        f.savefig(filename_out_png, bbox_inches='tight', dpi=600)

        # vertical average profile
        if args.avgwidth is not None:
            W_size = args.avgwidth
            W_0 = image.shape[1] // 2 - W_size // 2
            W_1 = W_0 + W_size
            profile = np.mean(image[:, W_0:W_1], axis=1)
        else:
            profile = image[:, image.shape[1] // 2]

        # normalize
        profile = profile/np.mean(profile)

        # Gaussian fit
        # mean, std = norm.fit(profile)
        fitter = modeling.fitting.LevMarLSQFitter()
        model = modeling.models.Gaussian1D(mean = 0.0)
        fitted_model = fitter(model, y, profile)

        RMS = np.mean((profile[20:-20]-fitted_model(y[20:-20]))**2) # Mean Squared Error
        print("RMS: ", RMS)
        RMS_norm = RMS/fitted_model.amplitude.value

        f, ax = plt.subplots()
        plt.plot(y, profile.astype("float"))
        plt.plot(y, fitted_model(y))
        plt.xlabel('Y [mm]')
        # leg = 'stddev: %.2f' %fitted_model.stddev.value
        # leg = 'normalized RMS: %.2f' % RMS_norm
        leg = 'RMS: %.3f' % RMS
        plt.legend([leg], loc='upper right')
        f.savefig(filename_out_profile, bbox_inches='tight', dpi=600)

    return

if __name__ == '__main__':
    main()
