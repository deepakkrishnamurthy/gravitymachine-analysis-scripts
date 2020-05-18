# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 12:43:23 2019

@author: Deepak
"""

import GravityMachineTrack
import imp
import msdanalyzer
import pandas as pd
import numpy as np
import os
import seaborn as sns
import cmocean
import PlotFunctions.PlotUtils as PlotUtils

from matplotlib import rcParams
from matplotlib import rc
import matplotlib.pyplot as plt


rc('font', family='sans-serif') 
rc('font', serif='Helvetica') 
rc('text', usetex='false') 
rcParams.update({'font.size': 24})

def TaylorFunction(t, v, tau):
    """ Taylor function:
    """
    return (2*v**2*tau*(t - tau*(1 - np.exp(-t/tau))))**(1/2)

# Color schemes for plotting
Organisms = ['Sea cucumber','Dendraster','Brittlestar','Acorn worm','Sea urchin','Starfish','Snail','Polychaete','Volvox','Noctiluca','Pyrocystis','Akashiwo', 'Ceratium_sp_noFork', 'Ceratium_sp_Fork', 'Euglena', 'Stentor']

OrganismScientific = {'Acorn worm':'S. californicum','Sea urchin':'S. purpuratus','Sea cucumber':'P. parvimensis','Brittlestar':'O. spiculata','Dendraster':'D. excentricus','Polychaete':'Owenia spp.','Starfish':'P. miniata','Snail':'C. fornicata','Noctiluca':'N. scintillans','Pyrocystis':'P. noctiluca','Volvox':'V. aureus', 'Akashiwo':'A. sanguinea', 'Ceratium_sp_noFork':'Ceratium sp.', 'Ceratium_sp_Fork': 'C. furca', 'Euglena' : 'E. gracilis', 'Stentor': 'S. coeruleus'}
MarkerStyle = {'Acorn worm':'o','Sea urchin':'s','Sea cucumber':'p','Brittlestar':'1','Dendraster':'D','Polychaete':'P','Starfish':'*','Snail':'8','Noctiluca':'x','Pyrocystis':'+','Volvox':'v','Akashiwo':'h', 'Ceratium_sp_noFork':'s', 'Ceratium_sp_Fork': '>', 'Euglena':'<', 'Stentor':'H'}

cmap = plt.get_cmap("tab20")
# cmap = cmocean.cm.matter
# cmap = plt.get_cmap("Accent", 18)

print(cmap.N)

# cmap = sns.hls_palette(len(Organisms), l=.3, s=.8)


cmap_new = []
ColorStyle={}
for ii in np.linspace(int(0),int(cmap.N),len(Organisms),dtype='int'):
    print(ii)
    cmap_new.append(cmap(ii))


for ii, org in enumerate(Organisms):
    ColorStyle[org] = cmap_new[ii]

save_folder = 'C:/Users/Deepak/Dropbox/GravityMachine/ExperimentResults/CombinedDatasets'

os.chdir(save_folder)

#data_file = 'C:/Users/Deepak/Dropbox/GravityMachine/ExperimentResults/CombinedDatasets/MSE_calculation_combined.csv'
data_file = 'C:/Users/Deepak/Dropbox/GravityMachine/ExperimentResults/Entropy_Analysis/MSE_calculation/Starfish_No light_MSE.csv'
    

df = pd.read_csv(data_file)

#df_larvae = df[df['Organism'].isin(['Acorn worm', 'Starfish', 'Dendraster', 'Sea urchin', 'Sea cucumber', 'Snail'])]
#
#df_singlecell = df[df['Organism'].isin(['Akashiwo', 'Ceratium_sp_Fork', 'Stentor'])]


## Plot the data using seaborn
#saveFile = 'Starfish_longTime_Entropy'
#plt.figure(figsize=(16,12))
#sns.lineplot(x="Scale factor times", y="MSE_Vz",
#             hue="Organism", style = 'Organism', markers = True, palette = sns.color_palette("Set2", 2), linewidth = 3, 
#             data=df , markersize = 20)
#
#sns.lineplot(x="Scale factor times", y="MSE_noise",
#             hue="Organism", style = 'Organism', markers = True, palette = sns.color_palette("Set2", 2), linewidth = 3, 
#             data = df)
#ax = plt.gca()
#ax.legend(markerscale=4)
#plt.title('MSE analysis')
#plt.savefig(saveFile + '_' + 'MSEanalysis.svg', dpi = 150)
#plt.savefig(saveFile + '_' + 'MSEanalysis.png', dpi = 150)


# Plot the data using seaborn
saveFile = 'Starfish_longTime_Entropy'
plt.figure(figsize=(16,12))
sns.lineplot(x="Scale factor times", y="MSE_Vz",
             hue="Organism", style = 'Organism', markers = '^', color = 'y', linewidth = 3, 
             data=df , markersize = 20)

sns.lineplot(x="Scale factor times", y="MSE_noise",
             hue="Organism", style = 'Organism', color = 'k', linewidth = 3, 
             data = df)
ax = plt.gca()
ax.legend(markerscale=4)
plt.title('MSE analysis')
plt.savefig(saveFile + '_' + 'MSEanalysis.svg', dpi = 150)
plt.savefig(saveFile + '_' + 'MSEanalysis.png', dpi = 150)




#plt.figure(figsize=(16,12))
#sns.lineplot(x="Scale factor times", y="MSE_Vz",
#             hue="Organism", style = 'Organism', markers = True, palette = sns.color_palette("Set2", 3), linewidth = 3, 
#             data =  df_larvae[df_larvae['Organism'].isin(['Sea urchin', 'Sea cucumber', 'Snail'])], markersize = 20)
#ax = plt.gca()
#ax.legend(markerscale=4)
#plt.ylim([0, 2.0])
#plt.title('Low complexity trajectories')
#plt.savefig('LowComplexityTrajectories.svg', dpi = 150)
#plt.savefig('LowComplexityTrajectories.png', dpi = 150)
#
#
#plt.figure(figsize=(16,12))
#sns.lineplot(x="Scale factor times", y="MSE_Vz",
#             hue="Organism", style = 'Organism', markers = True, palette = sns.color_palette("Set2", 3), linewidth = 3, 
#             data =  df_larvae[df_larvae['Organism'].isin(['Acorn worm', 'Starfish', 'Dendraster'])], markersize = 20)
#ax = plt.gca()
#ax.legend(markerscale=4)
#plt.title('High complexity trajectories')
#plt.ylim([0, 2.0])
#plt.savefig('HighComplexityTrajectories.svg', dpi = 150)
#plt.savefig('HighComplexityTrajectories.png', dpi = 150)



#plt.figure(figsize=(16,12))
#sns.lineplot(x="Scale factor times", y="MSE_Vz",
#             hue="Organism", style = 'Organism', markers = True, palette = sns.color_palette("Set2", 3), linewidth = 3, 
#             data =  df_singlecell, markersize = 20)
#
#sns.lineplot(x="Scale factor times", y="MSE_noise",
#             hue="Organism", style = 'Organism', markers = False, palette = sns.color_palette("Set2", 1), linewidth = 2, 
#             data = df_singlecell[df_singlecell['Organism'].isin(['Stentor'])] )
#ax = plt.gca()
#ax.legend(markerscale=4)
#plt.title('Single-cell trajectories')
#
#
#plt.ylim([0, 2.0])
#plt.savefig('SingleCellTrajectories.svg', dpi = 150)
#plt.savefig('SinglecellTrajectories.png', dpi = 150)


# Calculate the entropy averaged across time scales

#data_folder_larva = 'C:/Users/Deepak/Dropbox/GravityMachine/ExperimentResults/Entropy_Analysis/MSE_calculation_larva'
#
#files = os.listdir(data_folder_larva)
#
#print(files)
#
#mse_mean = np.zeros(len(files))
#mse_std = np.zeros(len(files))
#Organism =[]
#
#for ii, file in enumerate(files):
#    
#    data = pd.read_csv(os.path.join(data_folder_larva, file))
#    
#    Organism.append(data['Organism'][0])
#    
#    print(Organism)
#    
#    TrackNames = np.unique(data['Track name'])
#    
#    mse_integral = []
#    
#    for Tracks in TrackNames:
#        
#        mse_track = data['MSE_Vz'][data['Track name'] == Tracks]
#        
#        scaleFactor_track = data['Scale factor times'][data['Track name'] == Tracks]
#        
#        
#        bool_mask = scaleFactor_track>=1
#        
#        mse_track = mse_track[bool_mask]
#        scaleFactor_track = scaleFactor_track[bool_mask]
#        
#        # Integrate over the scales to get an average entropy value
#        
#        mse_integral.append(np.trapz(y = mse_track, x = scaleFactor_track))
#        
#    
#    mse_mean[ii] = np.nanmean(mse_integral)
#    mse_std[ii] = np.nanstd(mse_integral)
#    
#
#print(Organism)
#print(mse_mean)
#print(mse_std)
#    
    
        
        
        
        
        
        
        

