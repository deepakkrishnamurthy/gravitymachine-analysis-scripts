# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 17:56:27 2020

@author: Deepak
"""

import GravityMachineTrack
import imp
import pandas as pd
import numpy as np
import os
imp.reload(GravityMachineTrack)
import matplotlib.pyplot as plt
import FigureParameters

# Track file
#trackFile = 'F:/HopkinsEmbroyologyCourse_GoodData/2018_06_12/Polychaete_4D/Polychaete4/track_mod.csv'
trackFile = 'F:/Hopkins_2018_08_31_MarineSnow/MarSno2/track_cropped.csv'

track = GravityMachineTrack.gravMachineTrack(trackFile = trackFile, Tmin=0, Tmax = 5400, findDims = False)


# Organism depth vs Time
Time = track.T
Y_pos = track.df['Yobjet']



plt.figure()

plt.plot(Time, Y_pos, color = 'g', linestyle = '--')
plt.xlabel('Time (s)')
plt.ylabel('Y position (mm)')



# Vertical velocity vs Time
plt.figure()

plt.plot(track.T, track.Vz_smooth, color = 'b', linestyle = '--')

plt.xlabel('Time (s)')
plt.ylabel('Z velocity (mm/s)')

# Shift the track positions based on the wall position
Y_pos = -(Y_pos - max(Y_pos))

# Plot vertical velocity and distance from wall on the same plot
fig, ax1 = plt.subplots()

ax2 = ax1.twinx()
ax1.plot(track.T, track.Vz_smooth, 'b-')
ax2.plot(track.T, Y_pos, 'g-')

ax1.set_xlabel('Time (s)')
ax1.set_ylabel('Vertical velocity (mm/s)', color='b')
ax2.set_ylabel('Y-position (mm)', color='g')

plt.show()




# Distance to wall vs Vertical velocity
plt.figure()

plt.scatter(Y_pos, track.Vz_smooth, 20, 'b')

plt.xlabel('Distance to wall (mm)')
plt.ylabel('Z-velocity (mm/s)')

# Distance to wall vs X-velocity
plt.figure()

plt.scatter(Y_pos, track.Vx_smooth, 20, 'r')

plt.xlabel('Distance to wall (mm)')
plt.ylabel('X-velocity (mm/s)')


# Distance to wall vs Y-velocity
plt.figure()

plt.scatter(Y_pos, track.Vy_smooth, 20, 'g')

plt.xlabel('Distance to wall (mm)')
plt.ylabel('Y-velocity (mm/s)')