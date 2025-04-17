#!/bin/bash
cd /home/beatsbs/PycharmProjects/BEATS_raytracing/DMM

## activate reconstruction environment
eval "$(conda shell.bash hook)"
conda activate tomopy

## launch solara app
solara run BEATS_DMM_solara_GUI.py --host localhost --port 4567
