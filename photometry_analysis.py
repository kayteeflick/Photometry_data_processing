
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from photometry_functions import get_zdFF
from doricFunctions import read_doric, crop_exp, plot_data
      
## Load Data
folder = "C:\\Users\\kayte\Dropbox (MIT)\\Tonegawa Lab\\Projects\\BLA Dopamine\\Photometry\\Recordings\\photometry_raw_data\\rspoda_photometry\\rspoda_cfc\\"
filename = "20210621-cfc-rspoda25_0000.csv"
data = read_doric(folder, filename)

## Crop Data
crop_exp(data)

# compute z-score
data["zdff"] = get_zdFF(data["reference"],data["signal"],smooth_win=10,remove=0,lambd=5e4,porder=1,itermax=50)

## visualize data
shocks = [180,250,320]
figure = plot_data(data, shocks)

## Remove background/bleaching curve for the signal and reference

## Compute df/F for signal and reference

## Find the slope of signal_dF/F versus ref_dF/F

## Correct signal dF/F by subtracting the ref_dF/F by slope

## Lowpass filter

## Create peri-stimulus time histograms

## Z-score