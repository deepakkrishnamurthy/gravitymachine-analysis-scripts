# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 14:34:26 2019

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
import seaborn as sns
import cmocean

import matplotlib.pyplot as plt



# Color schemes for plotting
Organisms = ['SeaCucumber','Dendraster','BrittleStar','AcornWorm','SeaUrchin','Starfish','Snail','Polychaete','Volvox','Noctiluca','Pyrocystis','Akashiwo', 'Ceratium_sp_noFork', 'Ceratium_sp_Fork']

OrganismScientific = {'AcornWorm':'S. californicum','SeaUrchin':'S. purpuratus','SeaCucumber':'P. parvimensis','BrittleStar':'O. spiculata','Dendraster':'D. excentricus','Polychaete':'Owenia spp.','Starfish':'P. miniata','Snail':'C. fornicata','Noctiluca':'N. scintillans','Pyrocystis':'P. noctiluca','Volvox':'V. aureus', 'Akashiwo':'A. sanguinea', 'Ceratium_sp_noFork':'Ceratium sp.', 'Ceratium_sp_Fork': 'C. furca'}
MarkerStyle = {'AcornWorm':'o','SeaUrchin':'s','SeaCucumber':'p','BrittleStar':'1','Dendraster':'d','Polychaete':'P','Starfish':'*','Snail':'8','Noctiluca':'X','Pyrocystis':'+','Volvox':'v','Akashiwo':'hexagon1', 'Ceratium_sp_noFork':'triangle_up', 'Ceratium_sp_Fork': 'tri_up'}

cmap = plt.get_cmap("Dark2")
#    cmap = cmocean.cm.deep
    
cmap_new = []
ColorStyle={}
for ii in np.linspace(int(0),int(cmap.N),len(Organisms),dtype='int'):
    cmap_new.append(cmap(ii))


for ii, org in enumerate(Organisms):
    ColorStyle[org] = cmap_new[ii]

# Load all the MSD trajectories

msd_directory = 'C:/Users\Deepak/Dropbox/GravityMachine/ExperimentResults/MSD_Analysis/MSD_trajectories'

files = os.listdir(msd_directory)


fig, (ax1, ax2) = plt.subplots(nrows = 1, ncols = 2 , figsize = (16,8))

for file in files:
    
    data = pd.read_csv(file)
    delays = data['delays']
    
    MSD_X = data['MSD_X']
    MSD_Z = data['MSD_Z']
    stdev_X = data['stdev_X']
    stdev_Z = data['stdev_Z']
    
    PlotUtils.errorfill_sns(self.delays, self.MSD_X, self.stdev_X, ax = ax1, color = 'r', label = 'Horizontal (X)')

    PlotUtils.errorfill_sns(self.delays, self.MSD_Z, self.stdev_Z, ax = ax2, color = 'b', label = 'Vertical (Z)')
    
    ax1.set_xlabel('Time (s)')
    ax1.set_ylabel('MSD')
    
    ax2.set_xlabel('Time (s)')
    ax2.set_ylabel('MSD')