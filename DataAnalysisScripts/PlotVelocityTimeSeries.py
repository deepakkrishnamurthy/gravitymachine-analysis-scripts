# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 14:35:57 2020
Plot velocity time series
@author: Deepak
"""
import os
import imp
import numpy as np
#import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pickle
from IPython import get_ipython
#For plots in a separate window
get_ipython().run_line_magic('matplotlib', 'qt')
# For inline plot
#get_ipython().run_line_magic('matplotlib', 'inline')
# Analysis file path
sns.set(context='paper', style='white', palette='deep', font='sans-serif', font_scale=2, color_codes=True, rc=None)

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
rcParams.update({'font.size': 26})

path1 = 'E:/HOT317_GravityMachine_1/Day1/withoutBlock_655PM'

track1 = 'track000_0_0_analysis.csv'

path2 = 'E:/HOT317_GravityMachine_1/Day2/PIV_afterEddingTwoLayeredBlock_track1_chamberCenter'

track2 = 'track000_0_10_analysis.csv'

vel_file_1 = 'fluidVelocityTimeSeries_0_0.pkl'

vel_file_2 = 'fluidVelocityTimeSeries_0_10.pkl'



df_1 = pd.read_csv(os.path.join(path1, track1))
df_2 = pd.read_csv(os.path.join(path2, track2))


with open(os.path.join(path1, vel_file_1), 'rb') as f:  # Python 3: open(..., 'wb')
    imageIndex_array_1, u_avg_array_1, v_avg_array_1, u_std_array_1, v_std_array_1 = pickle.load(f)

print(u_avg_array_1)
    
with open(os.path.join(path2, vel_file_2), 'rb') as f:  # Python 3: open(..., 'wb')
    imageIndex_array_2, u_avg_array_2, v_avg_array_2, u_std_array_2, v_std_array_2 = pickle.load(f)
    
    

# Plot the two velocity time series
    

plt.figure()

ax1 = plt.plot(df_1['Time'], (u_avg_array_1**2 + v_avg_array_1**2)**(1/2), 'r-', linewidth = 2)

ax2 = plt.plot(df_2['Time'], (u_avg_array_2**2 + v_avg_array_2**2)**(1/2), 'b--', linewidth = 2)

plt.xlabel('Time (s)')
plt.ylabel('Spatially averaged velocity magnitude (mm/s)')
plt.show()

