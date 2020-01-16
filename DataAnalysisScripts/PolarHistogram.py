# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 16:28:46 2019

Plots the velocity orientations of a 3D track as a polar histogram.

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

# Color schemes for plotting
Organisms = ['Sea cucumber','Dendraster','Brittlestar','Acorn worm','Sea urchin','Starfish','Snail','Polychaete','Polychate', 'Volvox','Noctiluca','Pyrocystis','Akashiwo', 'Akashiwo_Up', 'Ceratium_sp_noFork', 'Ceratium_sp_Fork', 'Euglena', 'Stentor']

OrganismScientific = {'Acorn worm':'S. californicum','Sea urchin':'S. purpuratus','Sea cucumber':'P. parvimensis','Brittlestar':'O. spiculata','Dendraster':'D. excentricus','Polychaete':'Owenia spp.','Starfish':'P. miniata','Snail':'C. fornicata','Noctiluca':'N. scintillans','Pyrocystis':'P. noctiluca','Volvox':'V. aureus', 'Akashiwo':'A. sanguinea', 'Ceratium_sp_noFork':'Ceratium sp.', 'Ceratium_sp_Fork': 'C. furca', 'Euglena' : 'E. gracilis', 'Stentor': 'S. coeruleus'}
MarkerStyle = {'Acorn worm':'o','Sea urchin':'s','Sea cucumber':'p','Brittlestar':'1','Dendraster':'D','Polychaete':'P','Polychate':'.','Starfish':'*','Snail':'8','Noctiluca':'x','Pyrocystis':'+','Volvox':'v','Akashiwo':'h', 'Ceratium_sp_noFork':'s', 'Ceratium_sp_Fork': '>', 'Euglena':'<', 'Stentor':'H'}

cmap = plt.get_cmap("tab20")
# cmap = cmocean.cm.matter
# cmap = plt.get_cmap("Accent", 18)

print(cmap.N)

# cmap = sns.hls_palette(len(Organisms), l=.3, s=.8)





# DataFolder 
dataFolder = '/Users/deepak/Dropbox/GravityMachine/ExperimentResults/MSD_Analysis/VelocityTimeSeries'

files = os.listdir(dataFolder)

Organisms_to_plot = ['Acorn worm', 'Dendraster', 'Starfish', 'Polychaete', 'Polychate', 'Snail','Brittlestar','Sea cucumber','Sea urchin']

SingeCells = ['Akashiwo','Akashiwo_Up', 'Ceratium_sp_noFork','Euglena', 'Stentor']

Organisms_to_plot += SingeCells

fig = plt.figure(figsize=(16,16))

count = 1

print(files)

Theta_df = pd.DataFrame({'Organism':[], 'Condition':[], 'Theta mean':[], 'Theta std':[]})

cmap_new = []
ColorStyle={}
for ii in np.linspace(int(0),int(cmap.N),len(Organisms),dtype='int'):
    print(ii)
    cmap_new.append(cmap(ii))


for ii, org in enumerate(Organisms):
    ColorStyle[org] = cmap_new[ii]
    
    
for file in files:
    
    if(file.endswith(".csv")):
        data = pd.read_csv(os.path.join(dataFolder, file))
    
    
        Organism = data['Organism'][0]
        Condition = data['Condition'][0]
        
        if(Organism in Organisms_to_plot):
            
            Vx = data['VelocityX']
            Vy = data['VelocityY']
            Vz = data['VelocityZ']
            
            
            if(Organism == 'Stentor'):
                Vz = -Vz
            
            
        
            
            
            
            vector_magnitude = (Vx**2 + Vy**2 + Vz**2)**(1/2)
            
            Orientation_vectors = np.zeros((3, len(Vx)))
            
            Orientation_vectors[0,:] = Vx/vector_magnitude
            Orientation_vectors[1,:] = Vy/vector_magnitude
            Orientation_vectors[2,:] = Vz/vector_magnitude
            
            # Extract the orientation angle from the vertical
            
            Z_gravity = [0, 0, 1]
            
            cos_theta = Orientation_vectors[0,:]*Z_gravity[0] + Orientation_vectors[1,:]*Z_gravity[1] + Orientation_vectors[2,:]*Z_gravity[2]
    
            # Theta value in degrees
            theta_array = np.arccos(cos_theta)*(180/(np.pi))
            
            theta_mean, theta_std = (np.nanmean(theta_array), np.nanstd(theta_array))
            
            Theta_df = Theta_df.append(pd.DataFrame({'Organism':[Organism], 'Condition':[Condition], 'Theta mean':[theta_mean], 'Theta std':[theta_std]}))
            
            
            # Plot a polar histogram of the angles
            
            bin_size = 5
            
            a , b = np.histogram(theta_array, bins=np.arange(0, 180 + bin_size, bin_size), density = True)
            
            centers = np.deg2rad(np.ediff1d(b)//2 + b[:-1])
            
            
            
            
           
            ax = fig.add_subplot(4,4,count, projection='polar')
            if(Organism == 'Akashiwo_Up'):
                ax.bar(centers, a, width=np.deg2rad(bin_size), bottom=0.0, color = ColorStyle[Organism], edgecolor='r', linestyle = '--', alpha = 0.85)
            else:
                ax.bar(centers, a, width=np.deg2rad(bin_size), bottom=0.0, color = ColorStyle[Organism], edgecolor='k', linestyle = '-', alpha = 0.85)

            ax.set_theta_zero_location("N")
            ax.set_theta_direction(-1)
            if(count>1):
                ax.set_xticklabels([])
            ax.set_yticklabels([])
            ax.set_title(Organism, pad = 0, fontdict = {'fontSize':12})
            
            plt.show()
            
            count += 1
            
            
#
plt.savefig('VelocityOrientationPolarPlots_SingleCell_binsize_5.svg',dpi = 150)
plt.savefig('VelocityOrientationPolarPlots_SingleCell_binsize_5.png',dpi = 150)

Theta_df.to_csv('Theta statistic.csv')    

plt.figure()
sns.scatterplot(x = 'Organism', y = 'Theta mean', data = Theta_df) 
    

    
KS_test_df = pd.DataFrame({'Organism':[], 'Condition':[], 'KS test statistic':[], 'p-value':[]})

velocity_df = pd.DataFrame({'Organism':[], 'Condition':[], 'Vx mean':[], 'Vx std':[], 'Vz mean':[], 'Vz std':[]})
# Velcity distributions
for file in files:
    
    if(file.endswith(".csv")):
        data = pd.read_csv(os.path.join(dataFolder, file))
    
        
        
        Organism = data['Organism'][0]
        Condition = data['Condition'][0]
        
        if(Organism in Organisms_to_plot):
            
            Vx = data['VelocityX']
            Vy = data['VelocityY']
            Vz = data['VelocityZ']
            
            if(Organism == 'Stentor'):
                Vz = -Vz
                
            # calculate mean and SD of the velocities
            Vx_mean, Vx_std = (np.nanmean(Vx), np.nanstd(Vx))
            Vz_mean, Vz_std = (np.nanmean(Vz), np.nanstd(Vz))
            
            velocity_df = velocity_df.append(pd.DataFrame({'Organism':[Organism], 'Condition':[Condition], 'Vx mean':[Vx_mean], 'Vx std':[Vx_std], 'Vz mean':[Vz_mean], 'Vz std':[Vz_std]}))
            
            
            # Calculate the Kolmogorv-Smirnov distance between the distributions of Vx and Vy for each organism
            D, p_value = stats.ks_2samp(Vx, Vz)
            
            print('KS test statistic: {}, and p-value : {}'.format(D, p_value))
            
            D_subsample, p_value_subsample = stats.ks_2samp(Vx[:], Vz)
            
            KS_test_df = KS_test_df.append(pd.DataFrame({'Organism':[Organism], 'Condition':[Condition], 'KS test statistic':[D], 'p-value':[p_value]}))
            
            my_pal = {'VelocityZ_noWall': 'b' ,'VelocityX_noWall': 'r'}
            color_list = sns.color_palette("RdBu_r", 7)
            #my_pal = {'VelocityZ_noWall': color_list[0] ,'VelocityX_noWall': color_list[6]}

            xlim1 = -2
            xlim2 = 2
            decimals = 1

            plt.figure(figsize=(4,4))
            ax0 = sns.distplot(Vz,kde = True , color = 'b', norm_hist = True, hist_kws={"histtype": "bar","edgecolor":'w', "linewidth": 0.2, "alpha": 0.5, "label":'Vz'})
            #ax0 = sns.distplot(dataFrame_full.loc[dataFrame_full["Organism"] == Organism2,"VelocityZ"],  kde = True , color = 'k', norm_hist = True, hist_kws={"histtype": "bar","edgecolor":'w', "linewidth": 0.2, "alpha": 0.5, "label":'Vz'})

            ax1 = sns.distplot(Vx,kde = True , color = 'r', norm_hist = True, hist_kws={"histtype": "bar","edgecolor":'w', "linewidth": 0.2, "alpha": 0.5, "label":'Vx'})  

            plt.xlim([xlim1, xlim2])
            plt.xticks(np.round(np.linspace(xlim1,xlim2,5), decimals=decimals))
            #plt.savefig(os.path.join(saveSubFolder,Organism+'VelocityDistribution_FINAL.svg'))
            plt.xlabel('Velocity')
            plt.ylabel('PDF')
            plt.legend()
            ax = plt.gca()
#            ax.set_aspect(0.6)
            plt.title(Organism)
            plt.show()

            my_pal = {'VelocityZ': 'b' ,'VelocityX': 'r'}

            plt.savefig(Organism+'_'+Condition+'VelocityDist.svg', dpi = 150)
            plt.savefig(Organism+'_'+Condition+'VelocityDist.png', dpi = 150)

velocity_df.to_csv('VelocityStatistic_combined.csv')
KS_test_df.to_csv('KS_test statistic.csv') 