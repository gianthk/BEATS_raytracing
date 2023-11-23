#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
To run this macro:

ipython --pylab

In[1]: run DMM_run_edge_scan
In[2]: run_edge_scan(...)

"""
import matplotlib.pyplot as plt
import numpy as np
import os
import sys

from epics import Motor 

from lmfit.models import (LinearModel, StepModel, ConstantModel)

sys.path.insert(0, '/sls/X02DA/applications/tomoalign/')
from tomoalign.profile_fit_functions import fit_profile
from tomoalign.scan import *



def run_edge_scan(data_dir, scan_name, pitch_range=0.5, nsteps=40,
                  filter_position_in=0, filter_position_out=-17000):

    plt.interactive(True)

    mot_y = Motor('X02DA-ES1-SMP1:TRY')

    # scan without the filter
    outfile_template = os.path.join(data_dir, scan_name,
        "{}_noFilter".format(scan_name),
        "{}_noFilter.tif".format(scan_name))
    mot_y.move(filter_position_out, wait=True)
    dat_nf = motor2_scan('X02DA-OP-DMM:ROTX_C1', -pitch_range, pitch_range,
                'X02DA-OP-DMM:ROTX_C2', -pitch_range, pitch_range,
                nsteps, outfile_template, relative=True)
    plt.figure()
    plt.plot(dat_nf['pos_readback1'], dat_nf['avg_cts'], label="no filter")
    plt.pause(0.05)
    plt.show()
    
    # scan with filter
    outfile_template = os.path.join(data_dir, scan_name,
        "{}_withFilter".format(scan_name),
        "{}_withFilter.tif".format(scan_name))
    mot_y.move(filter_position_in, wait=True)
    dat_wf = motor2_scan('X02DA-OP-DMM:ROTX_C1', -pitch_range, pitch_range,
                'X02DA-OP-DMM:ROTX_C2', -pitch_range, pitch_range,
                nsteps, outfile_template, relative=True)
    plt.plot(dat_wf['pos_readback1'], dat_wf['avg_cts'], label="with filter")
    plt.pause(0.05)
    plt.show()

    edge_profile = np.array(dat_wf['avg_cts']) / np.array(dat_nf['avg_cts'])
    plt.figure()
    plt.plot(dat_wf['pos_readback1'], edge_profile, label="edge_profile")
    plt.legend()
    plt.xlabel('Pitch 1 [urad]')
    plt.ylabel('Intensity [a.u.]')
    plt.title('Edge scan: {}'.format(scan_name))
    plt.show()

    outfile = os.path.join(data_dir, scan_name, "edgescan.dat")
    with open(outfile, 'w') as fh:
        for (pos, val) in zip(dat_wf['pos_readback1'], edge_profile):
            fh.write("{:12.6g} {:12.6g}\n".format(pos, val))
    
    fit = fit_energy_edge(dat_wf['pos_readback1'], edge_profile)

def fit_energy_edge(x, y):

    # Fit the energy edge
    #plt.figure()

    diffy = np.diff(y)
    center_0 = x[np.argmin(diffy)]

    prf = StepModel(form='erf')
    p0 = prf.make_params(
        amplitude=y.min()-y.max(),
        center=center_0,
        sigma=(x.max()-x.min())/10.0
    )
    model = prf

    bkgr = LinearModel()
    model += bkgr
    p0 += bkgr.make_params(slope=0, intercept=0)

    fit_result = model.fit(y, p0, x=x)
    fit_result.plot()
    print(fit_result.fit_report())
    return fit_result

