# -*- coding: utf-8 -*-
"""
Created on Wed Nov  6 16:15:46 2019

@author: Deepak
"""
import numpy as np

def velMag(u,v):
    return np.sqrt((u**2+v**2))

def calcVorticity(x,y,u,v):
    
    
    u_copy = u.copy()
    v_copy = v.copy()
    
    u[np.isnan(u)] = 0
    
    v[np.isnan(v)] = 0
    
    
    dx = x[0,1]-x[0,0]
    
    dy = y[1,0] - y[0,0]
    
    u_dy = np.gradient(u, dy, axis = 1)
    
    v_dx = np.gradient(v, dx, axis = 0)
    
    vorticity = v_dx - u_dy
    
    vorticity_smooth = smoothData(vorticity,1.5)
    
    vorticity_smooth[np.isnan(u_copy)] = np.nan
    
    return vorticity_smooth