# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 16:22:13 2019

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
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import rc



rc('font', family='sans-serif') 
rc('font', serif='Helvetica') 
rc('text', usetex='false') 
rcParams.update({'font.size': 24})

def TaylorFunction(t, v, tau):
    """ Taylor function:
    """
    return (2*v**2*tau*(t - tau*(1 - np.exp(-t/tau))))**(1/2)

# Color schemes for plotting
Organisms = ['Sea cucumber','Dendraster','Brittlestar','Acorn worm','Sea urchin','Starfish','Snail','Polychaete_Day','Polychaete','Volvox','Noctiluca','Pyrocystis','Akashiwo', 'Akashiwo_Up', 'Ceratium_sp_noFork', 'Ceratium_sp_Fork', 'Euglena', 'Stentor']

OrganismScientific = {'Acorn worm':'S. californicum','Sea urchin':'S. purpuratus','Sea cucumber':'P. parvimensis','Brittlestar':'O. spiculata','Dendraster':'D. excentricus','Polychaete':'Owenia spp.','Starfish':'P. miniata','Snail':'C. fornicata','Noctiluca':'N. scintillans','Pyrocystis':'P. noctiluca','Volvox':'V. aureus', 'Akashiwo':'A. sanguinea', 'Ceratium_sp_noFork':'Ceratium sp.', 'Ceratium_sp_Fork': 'C. furca', 'Euglena' : 'E. gracilis', 'Stentor': 'S. coeruleus'}
MarkerStyle = {'Acorn worm':'o','Sea urchin':'s','Sea cucumber':'p','Brittlestar':'1','Dendraster':'D','Polychaete':'P', 'Polychaete_Day':'P','Starfish':'*','Snail':'8','Noctiluca':'x','Pyrocystis':'+','Volvox':'v','Akashiwo':'h','Akashiwo_Up':'h', 'Ceratium_sp_noFork':'s', 'Ceratium_sp_Fork': '>', 'Euglena':'<', 'Stentor':'H'}

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


# Load all the MSD trajectories

msd_folder = 'C:/Users/Deepak/Dropbox/GravityMachine/ExperimentResults/MSD_Analysis/MSD_trajectories'

# Load the Taylor Function fitted parameters

df_taylor = pd.read_csv('C:/Users/Deepak/Dropbox/GravityMachine/ExperimentResults/CombinedDatasets/TaylorFittingAnalysis_calculation_combined.csv')

print(df_taylor)
# os.chdir(msd_folder)

# msd_file = 'msd_combined_csv.csv'

# df = pd.read_csv(msd_file)

files = os.listdir(msd_folder)

fig, axes = plt.subplots(1, 2, figsize=(16, 8))


    

for file in files:
    
    if(file.endswith('csv')):
    
        print(file)
        data = pd.read_csv(os.path.join(msd_folder, file))

        Organism = data['Organism'][0]
        Condition = data['Condition'][0]
        
    

        print(Organism)



        delays = data['delays']

        RMSD_X = (data['MSD_X'])**(1/2)
        RMSD_Z = (data['MSD_Z'])**(1/2)

        stdev_X = (data['stdev_X'])**(1/2)
        stdev_Z = (data['stdev_Z'])**(1/2)
        
        OrgSize = data['OrgSize mean']
        
        RMSD_X_bySize = RMSD_X/(OrgSize)
        RMSD_Z_bySize = RMSD_Z/(OrgSize)

        stdev_X_bySize = data['stdev_X']/(OrgSize)
        stdev_Z_bySize = data['stdev_Z']/(OrgSize)

        print(MarkerStyle[Organism])
        print(ColorStyle[Organism])

    #     PlotUtils.errorfill_sns(delays, MSD_Z, stdev_Z, ax = axes[0], color = ColorStyle[Organism], marker = MarkerStyle[Organism], label = 'Vertical (Z)')
    #     PlotUtils.errorfill_sns(delays, MSD_X, stdev_X, ax = axes[1], color = ColorStyle[Organism], marker = MarkerStyle[Organism], label = 'Horizontal (X)')

    #     sns.relplot(x = delays, y = MSD_Z, ax = axes[0], color = ColorStyle[Organism], marker = MarkerStyle[Organism], label = 'Vertical (Z)')
    #     sns.relplot(x = delays, y = MSD_X, ax = axes[1], color = ColorStyle[Organism], marker = MarkerStyle[Organism], label = 'Horizontal (X)')

        # Plots of Markers
#        axes[0].plot(delays[::20],MSD_Z[::20], color = ColorStyle[Organism], marker = MarkerStyle[Organism], MarkerSize = 7, label = Organism, alpha = 0.7)
#        axes[1].plot(delays[::20], MSD_X[::20], color = ColorStyle[Organism], marker = MarkerStyle[Organism], MarkerSize = 7, label = Organism, alpha = 0.7)
#        
        # Plots of Lines
#         axes[0].plot(delays,MSD_Z, color = ColorStyle[Organism], linewidth = 2, label = Organism, alpha = 0.7)
#         axes[1].plot(delays, MSD_X, color = ColorStyle[Organism], linewidth = 2, label = Organism, alpha = 0.7)
        
        axes[0].plot(delays[::60],RMSD_Z_bySize[::60], color = ColorStyle[Organism], marker = MarkerStyle[Organism], MarkerSize = 10, label = Organism + Condition, alpha = 0.7)
        axes[1].plot(delays[::60], RMSD_X_bySize[::60], color = ColorStyle[Organism], marker = MarkerStyle[Organism], MarkerSize = 10, label = Organism + Condition, alpha = 0.7)
        
        
        
        axes[0].set_title('Size-scaled Root-Mean-Squared-Displacement - Vertical (Z)')
        axes[1].set_title('Size-scaled Root-Mean-Squared-Displacement - Horizontal (X)')

        legend = plt.legend(prop={'size': 10}, loc = 'best')
        

axes[0].set_xlabel('Time (s)')
axes[1].set_xlabel('Time (s)')

axes[0].set_ylabel('RMSD (Body lengths)')
axes[1].set_ylabel('RMSD (Body lengths)')

plt.savefig(os.path.join(msd_folder, 'Combined_RMSD_SizeScaled.svg'), dpi = 150)
plt.savefig(os.path.join(msd_folder,'Combined_RMSD_SizeScaled.png'), dpi = 150)
    

plt.show()


# Plots of Root-Mean-squared-Displacement

fig, axes = plt.subplots(1, 2, figsize=(16, 8))


    

for file in files:
    
    if(file.endswith('csv')):
    
        print(file)
        data = pd.read_csv(os.path.join(msd_folder, file))

        Organism = data['Organism'][0]
        Condition = data['Condition'][0]
        
        print(Organism)
        print(Condition)
        
        
#        print(df_taylor.loc[(df_taylor['Organism'] == Organism) & (df_taylor['Condition'] == Condition)])
        
#        tau_X = float(df_taylor.loc[(df_taylor['Organism'] == Organism) & (df_taylor['Condition'] == Condition)]['tau_X'])
        
#        v_X = float(df_taylor.loc[(df_taylor['Organism'] == Organism) & (df_taylor['Condition'] == Condition)]['v_X'])
        
        
#        print(tau_X)
#        print(v_X)


        



        delays = data['delays']

        RMSD_X = (data['MSD_X'])**(1/2)
        RMSD_Z = (data['MSD_Z'])**(1/2)

        stdev_X = (data['stdev_X'])**(1/2)
        stdev_Z = (data['stdev_Z'])**(1/2)
        
        
        OrgSize = data['OrgSize mean']


        print(MarkerStyle[Organism])
        print(ColorStyle[Organism])

    #     PlotUtils.errorfill_sns(delays, MSD_Z, stdev_Z, ax = axes[0], color = ColorStyle[Organism], marker = MarkerStyle[Organism], label = 'Vertical (Z)')
    #     PlotUtils.errorfill_sns(delays, MSD_X, stdev_X, ax = axes[1], color = ColorStyle[Organism], marker = MarkerStyle[Organism], label = 'Horizontal (X)')

    #     sns.relplot(x = delays, y = MSD_Z, ax = axes[0], color = ColorStyle[Organism], marker = MarkerStyle[Organism], label = 'Vertical (Z)')
    #     sns.relplot(x = delays, y = MSD_X, ax = axes[1], color = ColorStyle[Organism], marker = MarkerStyle[Organism], label = 'Horizontal (X)')

        # Plots of Markers
        axes[0].plot(delays[::30],RMSD_Z[::30], color = ColorStyle[Organism], marker = MarkerStyle[Organism], MarkerSize = 10, label = Organism + Condition, alpha = 0.7)
        axes[1].plot(delays[::30], RMSD_X[::30], color = ColorStyle[Organism], marker = MarkerStyle[Organism], MarkerSize = 10, label = Organism + Condition, alpha = 0.7)
#        
        # Plots of Lines
#         axes[0].plot(delays,MSD_Z, color = ColorStyle[Organism], linewidth = 2, label = Organism, alpha = 0.7)
#         axes[1].plot(delays, MSD_X, color = ColorStyle[Organism], linewidth = 2, label = Organism, alpha = 0.7)
        
        
#        print(TaylorFunction(t = delays, v =  v_X, tau =  tau_X))
        
#        axes[0].scatter(delays[::30],RMSD_Z[::30], 10, color = ColorStyle[Organism], marker = MarkerStyle[Organism], MarkerSize = 10, label = Organism, alpha = 0.7)
#        axes[1].scatter(delays[::30], RMSD_X[::30], 10, color = ColorStyle[Organism], marker = MarkerStyle[Organism], MarkerSize = 10, label = Organism, alpha = 0.7)
#        
#        axes[1].plot(delays, TaylorFunction(t = delays, v =  v_X, tau =  tau_X), linestyle = '-', linewidth = 2,  color = ColorStyle[Organism])
        
        axes[0].set_title('Root-Mean-Squared-Displacement - Vertical (Z)')
        axes[1].set_title('Root-Mean-Squared-Displacement - Horizontal (X)')

        legend = plt.legend(prop={'size': 10}, loc = 'best')
        
axes[1].hlines(7.5, min(delays), max(delays), linestyle = '--', color = 'k')

axes[0].set_xlabel('Time (s)')
axes[1].set_xlabel('Time (s)')

axes[0].set_ylabel('RMSD (mm)')
axes[1].set_ylabel('RMSD (mm)')

plt.savefig(os.path.join(msd_folder, 'Combined_RMSD.svg'), dpi = 150)
plt.savefig(os.path.join(msd_folder,'Combined_RMSD.png'), dpi = 150)
    

plt.show()

    
    

    
    
    
