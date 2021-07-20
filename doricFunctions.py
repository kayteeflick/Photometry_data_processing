# -*- coding: utf-8 -*-
"""
Created on Sun Jul 18 20:07:05 2021

@author: kayte
"""

import numpy as np
import pandas as pd
import os
import glob
from scipy import signal as ss
import matplotlib.pyplot as plt

# function to read doric csv file
def read_doric(folder, filename):
    print("### Reading csv file....")
    names1 = ['time_s', 'reference', 'signal', 'raw','dio1']
    df = pd.read_csv(folder+filename, header = 1, index_col = False, names = names1)
    data = {}
    for name in names1:
        data[name] = np.array(df[name])
    date, exp, name = filename.split("-")
    data["filename"] = filename
    data["name"] = name.split("_")[0]
    data["data"] = date
    data["exp_type"] = exp
    data["timestep"] = np.median(np.diff(data["time_s"]))
    data['freq'] = 1/data['timestep']
    
    print("Data from csv file fetched....")

    return data


#function to crop to beginning and end of experiment
def crop_exp(data):
    temp = np.where(data["dio1"] == 1)
    names1 = ['time_s', 'reference', 'signal', 'raw','dio1']
    for name in names1:
        data[name] = data[name][temp[0].min():temp[0].max()]
    
    # normalize timesteps   
    data['time_s'] = np.arange(1, len(data['time_s'])+1) * data['timestep']

def add_stim_to_plot(ax):
    ax.axvspan(start_stim, end_stim, alpha=shade_alpha,
               color='gray')
    ax.axvline(start_stim, alpha=lines_alpha, color='gray', ls='--')
    ax.axvline(end_stim, alpha=lines_alpha, color='gray', ls='--')
    
def plot_data(data, cues):
    '''

    Parameters
    ----------
    data : numpy array
        Numpy array containing formatted doric photometry data
    cues : iterable list
        Cue times

    Returns
    -------
    None.

    '''
    fig = plt.figure(figsize=(11, 8.5), dpi = 300)
    fig.suptitle(data['name'] + " " + data['exp_type'], fontsize = 16)
    ax1 = fig.add_subplot(311)
    ax1.plot(data["time_s"], data['signal']*100,'blue',linewidth=1.5)
    ax1.set_title("Raw Signal")
    ax1.set_xlabel("Time (s)")
    ax1.set_ylabel("mV")
    
    ax2 = fig.add_subplot(312)
    ax2.plot(data["time_s"], data["reference"]*100,'purple',linewidth=1.5)
    ax2.set_title("Raw Control")
    ax2.set_xlabel("Time (s)")
    ax2.set_ylabel("mV")
    
    ax3 = fig.add_subplot(313)
    ax3.plot(data["time_s"], data["zdff"], 'blue')
    ax3.set_xlabel("Time (s)")
    ax3.set_ylabel("Z-score")
    ax3.set_title('Z-Scored dFF')
    for xc in cues:
        ax3.axvline(xc, alpha=.8, color='gray', ls = '--')
    plt.subplots_adjust(hspace = .5)
    plt.show()
    
    return fig


# function to compute deltaF/F using fitted control channel and filtered signal channel
def deltaFF(signal, control):
    
	res = np.subtract(signal, control)
	normData = np.divide(res, control)
	#deltaFF = normData
	normData = normData*100

	return normData

# function to fit control channel to signal channel
def controlFit(control, signal):
    
	p = np.polyfit(control, signal, 1)
	arr = (p[0]*control)+p[1]
	return arr


# function to filter control and signal channel, also execute above two function : controlFit and deltaFF
# function will also take care if there is only signal channel and no control channel
# if there is only signal channel, z-score will be computed using just signal channel
def execute_controlFit_dff(control, signal):

	b = np.divide(np.ones((100,)), 100)
	a = 1

	if (control==0).all()==True:
		signal_smooth = ss.filtfilt(b, a, signal)
		return signal_smooth
	else:
		control_smooth = ss.filtfilt(b, a, control)
		signal_smooth = ss.filtfilt(b, a, signal)
		control_fit = controlFit(control_smooth, signal_smooth)
		norm_data = deltaFF(signal_smooth, control_fit)
		return norm_data
    
# helper function to compute z-score and deltaF/F

def helper_z_score(data):     #helper_z_score(control_smooth, signal_smooth):
	
	z_score_arr, norm_data_arr = np.array([]), np.array([])

	norm_data = execute_controlFit_dff(data["reference"], data["signal"])
	res = np.subtract(norm_data, np.nanmean(norm_data))
	z_score = np.divide(res, np.nanstd(norm_data))
	z_score_arr = np.concatenate((z_score_arr, z_score))
	norm_data_arr = np.concatenate((norm_data_arr, norm_data))

	return z_score_arr, norm_data_arr


# compute z-score and deltaF/F and save it to hdf5 file
def compute_z_score(filepath, inputParameters):

	print("Computing z-score for each of the data...")
	remove_artifacts = inputParameters['removeArtifacts']


	path_1 = c(os.path.join(filepath, 'control*'))
	path_2 = glob.glob(os.path.join(filepath, 'signal*'))
	
	path = sorted(path_1 + path_2)


	b = np.divide(np.ones((100,)), 100)
	a = 1

	if len(path)%2 != 0:
		raise Exception('There are not equal number of Control and Signal data')

	path = np.asarray(path).reshape(2,-1)

	for i in range(path.shape[1]):
		name_1 = ((os.path.basename(path[0,i])).split('.')[0]).split('_')
		name_2 = ((os.path.basename(path[1,i])).split('.')[0]).split('_')
		#dirname = os.path.dirname(path[i])
		
		if name_1[-1]==name_2[-1]:
			name = name_1[-1]
			control = read_hdf5('', path[0,i], 'data').reshape(-1)
			signal = read_hdf5('', path[1,i], 'data').reshape(-1)
			#control_smooth = ss.filtfilt(b, a, control)
			#signal_smooth = ss.filtfilt(b, a, signal)
			#_score, dff = helper_z_score(control_smooth, signal_smooth)
			z_score, dff = helper_z_score(control, signal, filepath, name, inputParameters)
			if remove_artifacts==True:
				write_hdf5(z_score, 'z_score_'+name, filepath, 'data')
				write_hdf5(dff, 'dff_'+name, filepath, 'data')
			else:
				write_hdf5(z_score, 'z_score_'+name, filepath, 'data')
				write_hdf5(dff, 'dff_'+name, filepath, 'data')
		else:
			raise Exception('Error in naming convention of files or Error in storesList file')


	print("z-score for the data computed.")
    
# function to compute z-score and deltaF/F using functions : compute_z_score and/or processTimestampsForArtifacts
def execute_zscore(folderNames, inputParameters, timeForLightsTurnOn, remove_artifacts, plot_zScore_dff, combine_data):

	storesListPath = []
	for i in range(len(folderNames)):
		if combine_data==True:
			storesListPath.append([folderNames[i][0]])
		else:
			filepath = folderNames[i]
			storesListPath.append(glob.glob(os.path.join(filepath, '*_output_*')))
	
	storesListPath = np.concatenate(storesListPath)
	
	for j in range(len(storesListPath)):
		filepath = storesListPath[j]
		storesList = np.genfromtxt(os.path.join(filepath, 'storesList.csv'), dtype='str', delimiter=',')

		if remove_artifacts==True:
			print("Removing Artifacts from the data and correcting timestamps...")
			compute_z_score(filepath, inputParameters)
			processTimestampsForArtifacts(filepath, timeForLightsTurnOn, storesList)
			visualizeControlAndSignal(filepath, remove_artifacts)
			print("Artifacts from the data are removed and timestamps are corrected.")
		else:
			compute_z_score(filepath, inputParameters)
			visualizeControlAndSignal(filepath, remove_artifacts)

		if plot_zScore_dff=='z_score':
			visualize_z_score(filepath)
		if plot_zScore_dff=='dff':
			visualize_dff(filepath)
		if plot_zScore_dff=='Both':
			visualize_z_score(filepath)
			visualize_dff(filepath)

	print("Signal data and event timestamps are extracted.")


def extractTsAndSignal(inputParametersPath):

	print("Extracting signal data and event timestamps...")

	with open(inputParametersPath) as f:	
		inputParameters = json.load(f)

	#storesList = np.genfromtxt(inputParameters['storesListPath'], dtype='str', delimiter=',')

	folderNames = inputParameters['folderNames']
	timeForLightsTurnOn = inputParameters['timeForLightsTurnOn']
	remove_artifacts = inputParameters['removeArtifacts']
	plot_zScore_dff = inputParameters['plot_zScore_dff']
	combine_data = inputParameters['combine_data']

	print("Remove Artifacts : ", remove_artifacts)
	print("Combine Data : ", combine_data)
	#print(type(remove_artifacts))

	if combine_data==False:
		execute_timestamp_correction(folderNames, timeForLightsTurnOn)
		execute_zscore(folderNames, inputParameters, timeForLightsTurnOn, remove_artifacts, plot_zScore_dff, combine_data)
	else:
		execute_timestamp_correction(folderNames, timeForLightsTurnOn)
		storesList = check_storeslistfile(folderNames)
		op_folder = combineData(folderNames, timeForLightsTurnOn, storesList)
		execute_zscore(op_folder, inputParameters, timeForLightsTurnOn, remove_artifacts, plot_zScore_dff, combine_data)
	