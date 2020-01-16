#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 27 15:02:45 2018

@author: deepak
"""

import csv as csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.markers as markers
from mpl_toolkits.mplot3d import Axes3D
import scipy.interpolate as interpolate
import scipy.signal as signal
import scipy.ndimage as ndimage
import matplotlib.patches as patches
from matplotlib.lines import Line2D
from numpy.polynomial import polynomial as Poly
import cmocean
import pickle
import math
import os
import pandas as pd
import math
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import cycler
plt.close("all")
import imp
import Track
import GravityMachineTrack
import FigureParameters
imp.reload(Track)
import seaborn as sns
import scipy.stats as stats
#plt.style.use(os.path.join(os.getcwd(),'gravmachinestyle.mplstyle'))
#==============================================================================
#                              Plot Parameters and Functions 
#==============================================================================
from matplotlib import rcParams
from matplotlib import rc
#rcParams['axes.titlepad'] = 20 
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
### for Palatino and other serif fonts use:
##rc('font',**{'family':'serif','serif':['Palatino']})
#rc('text', usetex=False)
#plt.rc('font', family='serif')

rc('font', family='sans-serif') 
rc('font', serif='Helvetica') 
rc('text', usetex='false') 
rcParams.update({'font.size': 18})
#==============================================================================
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

#dataFolder ='/Volumes/GRAVMACH1/VolvoxPhotoresponse/vvx25'

dataFolder = 'E:/VolvoxPhotoresponse/vvx25'


file = "track.csv"

*rest,orgName = os.path.split(dataFolder)
saveFolder = '/Users/deepak/Dropbox/GravityMachine/GravityMachineManuscript/EnsembleTrackStatistics'




print(50*'*')
print('Folder in which to save data: {}'.format(dataFolder))
print(50*'*')


#======================================================================
# Plots of Tracks
#======================================================================
Tmin = 0
Tmax = 4030
#OrgTrack1 = Track.GravMachineTrack(dataFolder,file,Tmin,Tmax)
OrgTrack1 = GravityMachineTrack.gravMachineTrack(trackFile = os.path.join(dataFolder, file), Tmin = Tmin, Tmax = Tmax)

#ImageFolder = '/Volumes/GRAVMACH1/VolvoxPhotoresponse/vvx25/images00012'
#saveFolder = '/Users/deepak/Dropbox/GravityMachine/GravityMachineManuscript/Figures/Figure6/Volvox_RisingEdge'

ImageFolder = 'E:/VolvoxPhotoresponse/vvx25/images00012'


#if (not os.path.exists(saveFolder)):
#    os.makedirs(saveFolder)
    
#OrgTrack1.makeMovie(ImageFolder, saveFolder)



custom_line = [Line2D([0], [0], color='gold', lw=4), Line2D([0], [0], color='k', lw=4) ]


RisingEdge = []
FallingEdge = []
counter = 0
for ii in range(OrgTrack1.trackLen - 1):
    if(OrgTrack1.LED_intensity[ii]==0 and OrgTrack1.LED_intensity[ii+1]>0):
       
        RisingEdge.append(ii)
        counter += 1
            
    elif(OrgTrack1.LED_intensity[ii]>0 and OrgTrack1.LED_intensity[ii+1]==0):
        FallingEdge.append(ii)
        
RisingEdge = np.array(RisingEdge)
FallingEdge = np.array(FallingEdge)

Index = []
Index_falling  = []
for ii in range(len(RisingEdge)-1):
    
    if(RisingEdge[ii+1] - RisingEdge[ii] > 100):
        Index.append(ii+1)
        
        
RisingEdge = RisingEdge[Index]

for ii in range(len(FallingEdge)):
    
    nearestRisingEdge = find_nearest(RisingEdge, FallingEdge[ii])
    
    if(abs(FallingEdge[ii]-nearestRisingEdge) > 200):
        Index_falling.append(ii)

        
FallingEdge = FallingEdge[Index_falling]

nBlinks = len(RisingEdge)
# Mid piint indcices between Rising and Falling edges where there is no trabnsition
noTransitionIndex = [int((FallingEdge[ii+1] + RisingEdge[ii])/2) for ii in range(nBlinks-1)]

print(noTransitionIndex)



Z_movingAvg = OrgTrack1.smoothSignal(OrgTrack1.ZobjWheel,20)


## Moving avg of the vertical velocity
#Tmin = 0
#Tmax = 4030
#mask = np.array(OrgTrack1.Time>=Tmin, dtype='bool') & np.array(OrgTrack1.Time<=Tmax, dtype = 'bool')

#------------------------------------------------------------------------------
# Volvox phototaxis plots
#------------------------------------------------------------------------------
# Plots of displacement vs time
#------------------------------------------------------------------------------
#fig, ax1 = plt.subplots(figsize=(6,4))
#
#ax1.fill_between(OrgTrack1.Time[mask],y1=np.max(OrgTrack1.ZobjWheel[mask])+0.1,y2 = np.min(OrgTrack1.ZobjWheel[mask])-0.1, where = OrgTrack1.LED_intensity[mask]>0,facecolor = 'gold', alpha = 0.4)
#ax1.fill_between(OrgTrack1.Time[mask],y1=np.max(OrgTrack1.ZobjWheel[mask])+0.1,y2 = np.min(OrgTrack1.ZobjWheel[mask])-0.1,where = OrgTrack1.LED_intensity[mask]==0,facecolor = 'k', alpha = 0.4)
#ax1.plot(OrgTrack1.Time[mask], OrgTrack1.ZobjWheel[mask] , 'b-', linewidth = 2)
##ax1.plot(OrgTrack1.Time[0:-1], OrgTrack1.Vz , 'b', linewidth = 1)
#
#
#
#    
##ax1.set_xlabel('Time (s)')
##ax1.set_ylabel('Z (mm)', color='b')
#
##ax1.legend(custom_line, ['Light ON', 'Light OFF'])
#ax1.set_aspect(2.5)
##ax1.set_xticks(range(Tmin, Tmax,5))
#ax1.set_xlim([np.min(OrgTrack1.Time[mask]), np.max(OrgTrack1.Time[mask])])
#ax1.set_ylim([np.min(OrgTrack1.ZobjWheel[mask])-0.1, np.max(OrgTrack1.ZobjWheel[mask])+0.1])
#plt.show()
#    
#------------------------------------------------------------------------------

Z_fast = OrgTrack1.ZobjWheel - Z_movingAvg



#------------------------------------------------------------------------------
# Plot of the Signal - Moving average
#------------------------------------------------------------------------------
fig, ax1 = plt.subplots()

  
ax1.fill_between(OrgTrack1.T,y1=np.max(OrgTrack1.ZobjWheel),y2 =np.min(OrgTrack1.ZobjWheel), where = OrgTrack1.LED_intensity>0,facecolor = 'gold', alpha = 0.4)
ax1.fill_between(OrgTrack1.T,y1=np.max(OrgTrack1.ZobjWheel),y2 =np.min(OrgTrack1.ZobjWheel),where = OrgTrack1.LED_intensity==0,facecolor = 'k', alpha = 0.4)
ax1.plot(OrgTrack1.T, OrgTrack1.ZobjWheel , 'b-', linewidth = 2)
ax1.plot(OrgTrack1.T[FallingEdge], Z_fast[FallingEdge], 'ro')
ax1.plot(OrgTrack1.T[RisingEdge], Z_fast[RisingEdge], 'go')
ax1.plot(OrgTrack1.T[noTransitionIndex], Z_fast[noTransitionIndex], 'bo')

ax1.set_xlabel('Time (s)')
ax1.set_ylabel('Z (mm)', color='b')

ax1.legend(custom_line, ['Light ON', 'Light OFF'])
ax1.set_xlim([0, np.max(OrgTrack1.T)])
ax1.set_ylim([np.min(OrgTrack1.ZobjWheel), np.max(OrgTrack1.ZobjWheel)])
plt.show()
#
##------------------------------------------------------------------------------
## Plot of the Vertical velocity vs Time
##------------------------------------------------------------------------------

#Vz_smooth = OrgTrack1.smoothSignal(OrgTrack1.Vz,OrgTrack1.window_time)
#
#Vz_smooth = np.insert(Vz_smooth,0,np.nan)
#
#fig, ax1 = plt.subplots(figsize=(6,4))
#
#
#  
#ax1.fill_between(OrgTrack1.Time[mask],y1=np.nanmax(Vz_smooth[mask])+0.1,y2 = np.nanmin(Vz_smooth[mask])-0.1, where = OrgTrack1.LED_intensity[mask]>0,facecolor = 'gold', alpha = 0.4)
#ax1.fill_between(OrgTrack1.Time[mask],y1=np.nanmax(Vz_smooth[mask])+0.1,y2 = np.nanmin(Vz_smooth[mask])-0.1,where = OrgTrack1.LED_intensity[mask]==0,facecolor = 'k', alpha = 0.4)
#ax1.plot(OrgTrack1.Time[mask], Vz_smooth[mask] , 'k', linewidth = 2)
#ax1.plot(OrgTrack1.Time[0:-1], OrgTrack1.Vz , 'b', linewidth = 1)
#
#
#
#    
#ax1.set_xlabel('Time (s)')
#ax1.set_ylabel('Vz (mm)', color='b')
#
##ax1.legend(custom_line, ['Light ON', 'Light OFF'])
#ax1.set_xlim([np.nanmin(OrgTrack1.Time[mask]), np.nanmax(OrgTrack1.Time[mask])])
#ax1.set_ylim([np.nanmin(Vz_smooth[mask])-0.1, np.nanmax(Vz_smooth[mask])+0.1])
##ax1.set_xticks(range(Tmin, Tmax,5))
#plt.show()

#------------------------------------------------------------------------------
# Create a dataframe with all the behavioral events aligned at the Rising Edge and Falling Edge
#------------------------------------------------------------------------------
T_window = 10

nWindow = math.ceil(int(T_window*OrgTrack1.samplingFreq)/2.)*2
T_blink = np.zeros(nWindow)



print('No:of blinks: {}'.format(nBlinks))
color = plt.cm.inferno(np.linspace(0, 1,nBlinks))
    
#    color = cmocean.cm.algae(np.linspace(0, 1,nBlinks))

plt.rcParams['axes.prop_cycle'] = cycler.cycler('color', color)


##------------------------------------------------------------------------------
## Create a dataframe with all the behavioral events aligned at the Rising Edge
##------------------------------------------------------------------------------

Vz_blink = np.zeros((nBlinks,nWindow),dtype='float')
# Save the mean and std of velocity in the neighbhourhood of the blink
Vz_DL_mean = np.zeros(nBlinks)
Vz_DL_std = np.zeros(nBlinks)
plt.figure(3)
for i,peakIndex in enumerate(RisingEdge):
    
    lower_index = peakIndex
    upper_index = peakIndex + int(nWindow)
    if(lower_index>=0 and upper_index<OrgTrack1.trackLen):
        T_blink = OrgTrack1.Time[lower_index:upper_index] - OrgTrack1.Time[peakIndex]
#        Z_blink[i,:] = Z_fast[lower_index:upper_index] - Z_fast[peakIndex]
        Vz_blink[i,:] = OrgTrack1.Vz[lower_index:upper_index] - np.nanmean(OrgTrack1.Vz[lower_index:upper_index])
        
        Vz_DL_std[i] = np.nanstd(Vz_blink[i,:])


df = pd.DataFrame({'Condition':[], 'Time': [], 'Event':[], 'Vz':[]})
        
for ii in range(len(RisingEdge)):
    df = df.append(pd.DataFrame({'Condition':np.repeat('D-L', len(T_blink), axis = 0), 'Time':T_blink, 'Event':np.repeat(ii, len(T_blink),axis=0),'Vz':np.transpose(Vz_blink[ii,:])}))
    

# Save the data frame for future use
    
df.to_csv(os.path.join(saveFolder, 'VolvoxVerticalVelocity_D_L.csv'))
    
# calculate the mean and std of the velocity variability
Vz_variability_DL_mean, Vz_variability_DL_std = np.nanmean(Vz_DL_std) , np.nanstd(Vz_DL_std) 

print('Mean +- Std in velocity variability:{}+- {}'.format(Vz_variability_DL_mean, Vz_variability_DL_std))

palette = sns.color_palette("Blues", len(RisingEdge))


print(df)

plt.figure()


ax = sns.lineplot(x = 'Time',y = 'Vz', ci=95,data = df, lw = 3, palette = palette, legend = False)
ax.set_title('D-L Transition')
#ax.set_aspect(1)

plt.show()
#
##------------------------------------------------------------------------------
## Create a dataframe with all the behavioral events aligned at the Falling Edge
##------------------------------------------------------------------------------
nWindow = math.ceil(int(T_window*OrgTrack1.samplingFreq)/2.)*2
T_blink = np.zeros(nWindow)

nBlinks = len(FallingEdge)

print('No:of blinks: {}'.format(nBlinks))
color = plt.cm.inferno(np.linspace(0, 1,nBlinks))
    
#    color = cmocean.cm.algae(np.linspace(0, 1,nBlinks))

plt.rcParams['axes.prop_cycle'] = cycler.cycler('color', color)

Vz_blink = np.zeros((nBlinks,nWindow),dtype='float')


Vz_LD_std = np.zeros(nBlinks)

plt.figure(3)
for i,peakIndex in enumerate(FallingEdge):
    
    lower_index = peakIndex
    upper_index = peakIndex + int(nWindow)
    if(lower_index>=0 and upper_index<OrgTrack1.trackLen):
        T_blink = OrgTrack1.Time[lower_index:upper_index] - OrgTrack1.Time[peakIndex]
#        Z_blink[i,:] = Z_fast[lower_index:upper_index] - Z_fast[peakIndex]
#        Z_blink[i,:] = OrgTrack1.ZobjWheel[lower_index:upper_index] - OrgTrack1.ZobjWheel[peakIndex]
        
        # Calculate the fluctuating componet of the velocity
        Vz_blink[i,:] = OrgTrack1.Vz[lower_index:upper_index] - np.nanmean(OrgTrack1.Vz[lower_index:upper_index])
        Vz_LD_std[i] = np.nanstd(Vz_blink[i,:])


df1 = pd.DataFrame({'Time': [], 'Blink Number': [], 'Vz':[]})
        
for ii in range(1,len(FallingEdge)-1):
    df1 = df1.append(pd.DataFrame({'Condition':np.repeat('L-D', len(T_blink), axis = 0), 'Time':T_blink, 'Event':np.repeat(ii, len(T_blink),axis=0),'Vz':np.transpose(Vz_blink[ii,:])}))
    
df1.to_csv(os.path.join(saveFolder, 'VolvoxVerticalVelocity_L_D.csv'))

print(Vz_LD_mean)
print(Vz_LD_std)
# calculate the mean and std of the velocity variability
Vz_variability_LD_mean , Vz_variability_LD_std = np.nanmean(Vz_LD_std) , np.nanstd(Vz_LD_std) 

print('Mean +- Std in velocity variability:{}+- {}'.format(Vz_variability_LD_mean, Vz_variability_LD_std))
    
palette = sns.color_palette("Blues", len(FallingEdge)-2)

plt.figure()


ax = sns.lineplot(x = 'Time',y = 'Vz',ci=95,data = df1, lw = 3, palette = palette, legend = False)
ax.set_title('L-D Transition')
#ax.set_aspect(1)

plt.show()

#------------------------------------------------------------------------------------------------------
# As a control create a data frame in neighbhourhoods where there is no transition in light intensity
#------------------------------------------------------------------------------------------------------
Vz_noTrans = np.zeros((len(noTransitionIndex),nWindow),dtype='float')

Vz_noTrans_std = np.zeros(len(noTransitionIndex))

for i,peakIndex in enumerate(noTransitionIndex):
    
    lower_index = peakIndex
    upper_index = peakIndex + int(nWindow)
    if(lower_index>=0 and upper_index<OrgTrack1.trackLen):
        T_blink = OrgTrack1.Time[lower_index:upper_index] - OrgTrack1.Time[peakIndex]
#        Z_blink[i,:] = Z_fast[lower_index:upper_index] - Z_fast[peakIndex]
#        Z_blink[i,:] = OrgTrack1.ZobjWheel[lower_index:upper_index] - OrgTrack1.ZobjWheel[peakIndex]
        Vz_noTrans[i,:] = OrgTrack1.Vz[lower_index:upper_index] - np.nanmean(OrgTrack1.Vz[lower_index:upper_index])
        Vz_noTrans_std[i] = np.nanstd(Vz_noTrans[i,:])


# calculate the mean and std of the velocity variability
Vz_variability_noTrans_mean , Vz_variability_noTrans_std = np.nanmean(Vz_noTrans_std) , np.nanstd(Vz_noTrans_std) 

print('Mean +- Std in velocity variability:{}+- {}'.format(Vz_variability_noTrans_mean, Vz_variability_noTrans_std))

df_noTrans = pd.DataFrame({'Time': [], 'sub track': [], 'Vz':[]})
        
for ii in range(len(noTransitionIndex)):
    df_noTrans = df_noTrans.append(pd.DataFrame({'Condition':np.repeat('No Transition', len(T_blink), axis = 0), 'Time':T_blink, 'Event':np.repeat(ii, len(T_blink),axis=0),'Vz':np.transpose(Vz_noTrans[ii,:])}))
    
df_noTrans.to_csv(os.path.join(saveFolder, 'VolvoxVerticalVelocity_noTransition.csv'))

plt.figure()


ax = sns.lineplot(x = 'Time',y = 'Vz',ci=95,data = df_noTrans, lw = 3, palette = palette, legend = False)

ax.set_title('No Transition')
#ax.set_aspect(1)

plt.show()


df_combined = pd.concat([df, df1, df_noTrans])

df_combined.to_csv(os.path.join(saveFolder, 'VolvoxVerticalVelocity_Combined.csv'))





##------------------------------------------------------------------------------
# Compare the velocity distributions during D-L, L-D and no Transition cases
##------------------------------------------------------------------------------


my_pal = {'D-L': 'gold' ,'V': 'r'}
color_list = sns.color_palette("RdBu_r", 7)
#my_pal = {'VelocityZ_noWall': color_list[0] ,'VelocityX_noWall': color_list[6]}

xlim1 = -2
xlim2 = 2
decimals = 1

plt.figure(figsize=(4.5,4))
ax0 = sns.distplot(df_combined.loc[df_combined["Condition"] == 'D-L',"Vz"],  kde = True , color = 'gold', norm_hist = True, hist_kws={"histtype": "bar","edgecolor":'w', "linewidth": 0.2, "alpha": 0.5, "label":'D-L'})
#ax0 = sns.distplot(dataFrame_full.loc[dataFrame_full["Organism"] == Organism2,"VelocityZ"],  kde = True , color = 'k', norm_hist = True, hist_kws={"histtype": "bar","edgecolor":'w', "linewidth": 0.2, "alpha": 0.5, "label":'Vz'})

ax1 = sns.distplot(df_combined.loc[df_combined["Condition"] == 'L-D',"Vz"],  kde = True , color = 'k', norm_hist = True, hist_kws={"histtype": "bar","edgecolor":'w', "linewidth": 0.2, "alpha": 0.5, "label":'L-D'})

ax2 = sns.distplot(df_combined.loc[df_combined["Condition"] == 'No Transition',"Vz"],  kde = True , color = 'b', norm_hist = True, hist_kws={"histtype": "bar","edgecolor":'w', "linewidth": 0.2, "alpha": 0.5, "label":'No Transition'})

plt.xlim([xlim1, xlim2])
plt.xticks(np.round(np.linspace(xlim1,xlim2,5), decimals=decimals))
#plt.savefig(os.path.join(saveSubFolder,Organism+'VelocityDistribution_FINAL.svg'))
plt.xlabel('Velocity')
plt.ylabel('PDF')
plt.legend()
ax = plt.gca()
ax.set_aspect(0.6)
plt.show()


# Do a KS test between D-L and L-D and compare it between D-L and no Transition and L-D no Transition

D_DLvsLD, p_value_DLvsLD = stats.ks_2samp(df_combined.loc[df_combined["Condition"] == 'D-L',"Vz"], df_combined.loc[df_combined["Condition"] == 'L-D',"Vz"])

D_LDvsNoTrans, p_value_LDvsNoTrans = stats.ks_2samp(df_combined.loc[df_combined["Condition"] == 'L-D',"Vz"], df_combined.loc[df_combined["Condition"] == 'No Transition',"Vz"])

D_DLvsNoTrans, p_value_DLvsNoTrans = stats.ks_2samp(df_combined.loc[df_combined["Condition"] == 'D-L',"Vz"], df_combined.loc[df_combined["Condition"] == 'No Transition',"Vz"])

print('D-L vs L-D')
print('KS test statistic: {}, and p-value : {}'.format(D_DLvsLD, p_value_DLvsLD))

print('L-D vs No Transition')
print('KS test statistic: {}, and p-value : {}'.format(D_LDvsNoTrans, p_value_LDvsNoTrans))

print('D-L vs No Transition')

print('KS test statistic: {}, and p-value : {}'.format(D_DLvsNoTrans, p_value_DLvsNoTrans))
##------------------------------------------------------------------------------
## Align the responses for all Dark to Light transitions
##------------------------------------------------------------------------------
#
#for ii in range(len(RisingEdge)):
#    
#  
#    
#    fig, ax1 = plt.subplots()
#
#  
##    ax1.fill_between(OrgTrack1.Time,y1=np.max(Z_fast),y2 =np.min(Z_fast), where = OrgTrack1.LED_intensity>0,facecolor = 'gold', alpha = 0.5)
##    ax1.fill_between(OrgTrack1.Time,y1=np.max(Z_fast),y2 =np.min(Z_fast),where = OrgTrack1.LED_intensity==0,facecolor = 'k', alpha = 0.2)
##    ax1.plot(OrgTrack1.Time, OrgTrack1.ZobjWheel , 'b-', linewidth = 2)
#    
#    ax1.plot(OrgTrack1.Time, OrgTrack1.ZobjWheel , 'b-', linewidth = 2)
#
#    ax2 = ax1.twinx()
#    
#    ax2.plot(OrgTrack1.Time, OrgTrack1.LED_intensity, color = 'gold')
#    
#        
#    ax1.set_xlabel('Time (s)')
#    ax1.set_ylabel('V_Z (mm/s)', color='b')
#    ax2.set_ylabel('LED intensity')
#    
##    ax1.legend(custom_line, ['Light ON', 'Light OFF'])
##    ax1.set_xlim([0, np.max(OrgTrack1.Time)])
##    ax1.set_ylim([np.min(Z_fast), np.max(Z_fast)])
#    plt.show()
#    
#    fig, ax1 = plt.subplots()
#
#  
##    ax1.fill_between(OrgTrack1.Time,y1=np.max(Z_fast),y2 =np.min(Z_fast), where = OrgTrack1.LED_intensity>0,facecolor = 'gold', alpha = 0.5)
##    ax1.fill_between(OrgTrack1.Time,y1=np.max(Z_fast),y2 =np.min(Z_fast),where = OrgTrack1.LED_intensity==0,facecolor = 'k', alpha = 0.2)
##    ax1.plot(OrgTrack1.Time, OrgTrack1.Xobj , 'r-', linewidth = 2)
#    
#    ax1.plot(OrgTrack1.Time, OrgTrack1.Vx , 'b-', linewidth = 2)
#
#    ax2 = ax1.twinx()
#    
#    ax2.plot(OrgTrack1.Time, OrgTrack1.LED_intensity, color = 'gold')
#    
#        
#    ax1.set_xlabel('Time (s)')
#    ax1.set_ylabel('V_x (mm/s)', color='b')
#    ax2.set_ylabel('LED intensity')
#    
##    ax1.legend(custom_line, ['Light ON', 'Light OFF'])
##    ax1.set_xlim([0, np.max(OrgTrack1.Time)])
##    ax1.set_ylim([np.min(Z_fast), np.max(Z_fast)])
#    plt.show()