# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 18:12:49 2019

@author: Deepak
"""

import GravityMachineTrack
import imp
import msdanalyzer
import pandas as pd
import numpy as np
import os
imp.reload(GravityMachineTrack)
imp.reload(msdanalyzer)

from matplotlib import rcParams
from matplotlib import rc
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats


rc('font', family='sans-serif') 
rc('font', serif='Helvetica') 
rc('text', usetex='false') 
rcParams.update({'font.size': 16})


#File1 = '/Users/deepak/Dropbox/GravityMachine/ExperimentResults/MSD_Analysis/VelocityTimeSeries/Polychaete_Day__VelocityTimeSeries.csv'
#File2 = '/Users/deepak/Dropbox/GravityMachine/ExperimentResults/MSD_Analysis/VelocityTimeSeries/Polychate_Night__VelocityTimeSeries.csv'
#
#data1 = pd.read_csv(os.path.join(File1))
#data2 = pd.read_csv(os.path.join(File2))
#
#Vz_Day = data1['VelocityX']
#Vz_Night = data2['VelocityX']
#
#D, p_value = stats.ks_2samp(Vz_Day, Vz_Night)
##
#print('KS test statistic: {}, and p-value : {}'.format(D, p_value))

#------------------------------------------------------------------------------
# Testing the KS test on random sub-samples of data drawn from the original data
#------------------------------------------------------------------------------
#data_1 = data.sample(frac = 0.1)
#data_2 = data.sample(frac = 0.1)
#
#Vz_1 = data_1['VelocityZ']
#Vz_2 = data_2['VelocityZ']
#
#Vx_1 = data_1['VelocityX']
#Vx_2 = data_2['VelocityX']
#
#D, p_value = stats.ks_2samp(Vx_1, Vx_2)
#
#print('KS test statistic: {}, and p-value : {}'.format(D, p_value))
#------------------------------------------------------------------------------

File1 = '/Users/deepak/Dropbox/GravityMachine/ExperimentResults/MSD_Analysis/VelocityTimeSeries/Akashiwo_Ambient light_Down__VelocityTimeSeries.csv'
File2 = '/Users/deepak/Dropbox/GravityMachine/ExperimentResults/MSD_Analysis/VelocityTimeSeries/Akashiwo_Ambient light_Up__VelocityTimeSeries.csv'

data1 = pd.read_csv(os.path.join(File1))
data2 = pd.read_csv(os.path.join(File2))


xlim1 = -2
xlim2 = 2
decimals = 1

plt.figure(figsize=(4.5,4))
ax0 = sns.distplot(data1["VelocityZ"],  kde = True , color = 'b', norm_hist = True, hist_kws={"histtype": "bar","edgecolor":'w', "linewidth": 0.2, "alpha": 0.5, "label":'Vz'})
ax0 = sns.distplot(data2["VelocityZ"],  kde = True , color = 'c', norm_hist = True, hist_kws={"histtype": "bar","edgecolor":'w', "linewidth": 0.2, "alpha": 0.5, "label":'Vz'})

#ax0 = sns.distplot(dataFrame_full.loc[dataFrame_full["Organism"] == Organism2,"VelocityZ"],  kde = True , color = 'k', norm_hist = True, hist_kws={"histtype": "bar","edgecolor":'w', "linewidth": 0.2, "alpha": 0.5, "label":'Vz'})

ax1 = sns.distplot(data1["VelocityX"],  kde = True , color = 'r', norm_hist = True, hist_kws={"histtype": "bar","edgecolor":'w', "linewidth": 0.2, "alpha": 0.5, "label":'Vx'})  
ax1 = sns.distplot(data2["VelocityX"],  kde = True , color = 'm', norm_hist = True, hist_kws={"histtype": "bar","edgecolor":'w', "linewidth": 0.2, "alpha": 0.5, "label":'Vx'})  


plt.xlim([xlim1, xlim2])
plt.xticks(np.round(np.linspace(xlim1,xlim2,5), decimals=decimals))
#plt.savefig(os.path.join(saveSubFolder,Organism+'VelocityDistribution_FINAL.svg'))
plt.xlabel('Velocity')
plt.ylabel('PDF')
plt.legend()
ax = plt.gca()
ax.set_aspect(0.6)
plt.show()

my_pal = {'VelocityZ': 'b' ,'VelocityX': 'r'}

#plt.savefig(self.Organism+'_'+self.Condition+'.svg', dpi = 150)

