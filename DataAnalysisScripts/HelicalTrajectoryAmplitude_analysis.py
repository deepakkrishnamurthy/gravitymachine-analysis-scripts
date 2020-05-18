# -*- coding: utf-8 -*-
"""
Created on Sat Mar 28 21:20:42 2020

@author: Deepak
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 15:56:35 2019

@author: Deepak
"""
'''
Created on Mon Sep 23 16:28:46 2019

@author: Deepak
'''

import GravityMachineTrack
import imp
import msdanalyzer
import pandas as pd
import numpy as np
import os
imp.reload(GravityMachineTrack)
imp.reload(msdanalyzer)

import matplotlib.pyplot as plt

    

def get_HelixAmplitude(track, avg_window):
    Time = track.T
    X_pos = track.X
    X_pos_slow = track.smoothSignal(track.X, avg_window)
    X_pos_fast = X_pos - X_pos_slow
    
    plt.figure()

    plt.plot(Time, X_pos, color = 'r', linestyle = '-',label = 'raw signal')
    plt.plot(Time, X_pos_slow, color = 'm', linestyle = '--',label = 'slow component')
    plt.plot(Time, X_pos_fast, color = 'k', linestyle = '-', label = 'fast component')
    plt.xlabel('Time (s)')
    plt.ylabel('X position (mm)')
    
    signal_std = np.std(X_pos_fast)
    
    signal_amp = 2*(2**(1/2))*signal_std

    print('Std dev in fast signal: {} mm'.format(signal_std))

    print('Estimated amplitude: {} mm'.format(signal_amp))
    
    return signal_amp
    
#------------------------------------------------------------------------------------------------------------------------------    
# Define some constants for the data analysis:

minTrackDuration = 0
    
    
# Folder in which to save analysis results
saveFolder = 'C:/Users/deepak/Dropbox/GravityMachine/ExperimentResults/HelicalTrajectory_Analysis'

if(not os.path.exists(saveFolder)):
    os.makedirs(saveFolder)
    
    
#analysis_file = 'C:/Users/Deepak/Dropbox/GravityMachine/GravityMachineAnalysis_Scripts/TracksUsedForAnalysis/EntropyAnalysis_TracksUsed/Starfish_LongestTrack.csv'

analysis_file = 'C:/Users/Deepak/Dropbox/GravityMachine/GravityMachineAnalysis_Scripts/TracksUsedForAnalysis/HelicalAmplitude_analysis/BatchProcess.csv'

analysis_df = pd.read_csv(analysis_file)


#pairwise_conditions = [(x,y) for x in analysis_df['Organism'] for y in analysis_df['Condition']]
#
#pairwise_conditions = []
#
#for ii, x in enumerate(analysis_df['Organism']):
#    
#    for jj, y in enumerate(analysis_df['Condition']):
        
        
        
        
    
    

used = set()
unique_conditions = [x for x in analysis_df['Organism'] if x not in used and (used.add(x) or True)]
print(50*'-')
print('Unique condition : {}'.format(unique_conditions))
print(50*'-')

ensemble_tracks = []
TrackArray = []

nUniqueConditions = len(unique_conditions)

# Assemble the list of unique conditions i.e unique combinations of Organisms, Condition. Each Unique Condition will in general have several tracks

OrgDims = []
overwrite = False


for ii in range(nUniqueConditions):
    df = pd.DataFrame({'trackFolder':[], 'Organism':[], 'Condition':[], 'Size':[],'Helix amplitude':[], 'Chamber depth':[]})

    curr_Organism = unique_conditions[ii]
    
    print(curr_Organism)
  
    
    bool_Org = analysis_df['Organism'] == curr_Organism
    
    saveFile = os.path.join(saveFolder, curr_Organism +'_HelixParameters.csv')
    # Check if the analysis data is already present on disk for this Organism, Condition
    if(overwrite is True or not os.path.exists(saveFile)):
        
        tracks_sameCondition = analysis_df[bool_Org]
        
        nTracksSameCondition = len(tracks_sameCondition)
        
        
        
        for jj in range(nTracksSameCondition):
            
            Track_df = tracks_sameCondition.iloc[jj]
            
            avg_window = Track_df['Averaging window']
        
            chamber_depth = Track_df['Chamber depth']
            
            Condition = Track_df['Condition']
            
            Track_folder = Track_df['trackFolder']
            
            use_postprocessed_data = Track_df['Postprocessed']
            
            Pixelpermm = Track_df['Pixelpermm']
            
            print(Track_df)
            
            full_path = os.path.join(Track_df['rootFolder'],Track_df['trackFolder'], Track_df['trackFile'])
            
            
            track = GravityMachineTrack.gravMachineTrack(trackFile = full_path , organism = Track_df['Organism'], condition = Track_df['Condition'], Tmin = Track_df['Tmin'], Tmax = Track_df['Tmax'], findDims = True, pixelPermm = Pixelpermm, use_postprocessed = use_postprocessed_data)
            
            try:
                OrgDims.append(track.OrgDim)
                print('Organism dimensions: {} mm'.format(track.OrgDim))
            except:
                temp = float(input('Enter Organism dimensions in mm'))
                track.OrgDim = temp
                    
            # Filter tracks based on min Track Duration
            if(track.trackDuration >= minTrackDuration):
                TrackArray.append(track)
            
            print(avg_window)
            helix_amp = get_HelixAmplitude(track, avg_window)
            
            print(helix_amp)
            
            df = df.append(pd.DataFrame({'trackFolder':[Track_folder], 'Organism':[curr_Organism], 'Condition':[Condition], 'Size':[track.OrgDim],'Helix amplitude':[helix_amp], 'Chamber depth':[chamber_depth]}))
            
        
        print(df)
            
        df.to_csv(saveFile)
   
#
#    
#    
#    
#        
#    
#    
#      
