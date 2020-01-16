# -*- coding: utf-8 -*-
"""
Created on Wed Nov  6 16:13:14 2019
@author: Deepak
General PIV analysis for sequence of image pairs

"""
import openpiv.tools
import openpiv.scaling
import openpiv.process
import numpy as np
import openpiv.validation
import openpiv.filters
import csv as csv
import openpiv.tools
import openpiv.process
import openpiv.scaling
import openpiv.validation
import openpiv.filters
import cv2
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import cmocean
import os
import time
import scipy
import scipy.ndimage as ndimage
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import pickle
plt.close("all")
from PIL import Image
import imp
import GravityMachineTrack
imp.reload(GravityMachineTrack)
import FigureParameters
from mpl_toolkits.axes_grid1 import make_axes_locatable
from smoothn import *
from rangeslider_functions import *
from matplotlib.patches import Circle, Wedge, Polygon
import PIVanalysis.FlowFunctions
import PIVanalysis.PIV_Functions as PIV_Functions
import pandas as pd


trackFile = 'E:/HopkinsEmbroyologyCourse_GoodData/2018_06_12/Starfish/StarFish6/track_mod.csv'

Track = GravityMachineTrack.gravMachineTrack(trackFile = trackFile, findDims = True)
# PIV frames based on time
Tmin = 0
Tmax = 10

# PIV frames based on Image indices
ImgMin = 11372
nImages = 10
Imgmax = 0


def doPIVBatch(Track = None, ImagePairs = None, dT_array = None, win_size = 64, overlap = 32, searchArea = 64, save = False, masking = True):
    
    PIV_Folder = os.path.join(Track.path, 'PIV_data')
    if (save is True):
        
        if(not os.path.exists(PIV_Folder)):
            os.makedirs(PIV_Folder)
            
        # Save the configrations used for the PIV analysis
        # things to save: Window size, overlap, Search Area, Image Size, Pixel Size
        config_df = pd.DataFrame({'Window size':[win_size], 'Overlap':[overlap], 'Search area':[searchArea],'Pixels per mm':[Track.pixelPermm],'ImW':[Track.imW], 'imH':[Track.imH]})
    
        config_df.to_csv(os.path.join(PIV_Folder, 'PIV_configs.csv'))
            
    
   
    
    
    for ii, images in enumerate(ImagePairs):
        
        image_a, image_b = images
        
        print(image_a)
        print(image_b)
        #--------------------------------------------------------------------------
        # Load the frame-pair into memory
        #--------------------------------------------------------------------------
        frame_a_color = cv2.imread(os.path.join(Track.path, Track.image_dict[image_a], image_a))
        frame_b_color = cv2.imread(os.path.join(Track.path, Track.image_dict[image_b], image_b))
        
#        cv2.imshow('frame a',frame_a_color)
#        cv2.waitKey(0)
#        cv2.imshow('frame_b', frame_b_color)
#        cv2.waitKey(0)
        
        saveFile = os.path.join(PIV_Folder,'PIV_' + image_a[:-4]+'.pkl')
        
        deltaT = dT_array[ii]
        # Do the PIV calculation
        x,y,u,v, sig2noise = PIV_Functions.doPIV(frame_a_color,frame_b_color, dT = deltaT, win_size = win_size, overlap = overlap, searchArea = searchArea, apply_clahe = False)
        
        
#        u, v = PIV_Functions.pivPostProcess(u,v,sig2noise, sig2noise_min = 1.5, smoothing_param = 0)

        
        u,v = (PIV_Functions.data2RealUnits(data = u,scale = 1/(Track.pixelPermm)), PIV_Functions.data2RealUnits(data = v,scale = 1/(Track.pixelPermm)))
            
            #--------------------------------------------------------------------------
            # Threshold the image to extract the object regions
            #--------------------------------------------------------------------------
        if(masking is True):
            Contours = PIV_Functions.findContours(frame_a_color,Track.threshLow,Track.threshHigh,'largest')
        else:
            Contours = np.nan
                
        # Plot the PIV data
        PIV_Functions.plotPIVdata(frame_a_color,x,y,u,v, orgContour=Contours, Centroids = None, pixelPermm = Track.pixelPermm)
        
        if(save is True):
            with open(saveFile, 'wb') as f:  # Python 3: open(..., 'wb')
                pickle.dump((x, y , u, v, Contours), f)
        
    

def createImagePairArray(Track, startImage = ImgMin, nImages = nImages, stopImage = None):
    
    startImage_str = 'IMG_'+'{:05d}'.format(startImage)+'.tif'
    stopImage = startImage + nImages
    stopImage_str = 'IMG_'+'{:05d}'.format(stopImage)+'.tif'
    
    
    
    
#        print(Track1.ImageName[imageIndex] == startImage_str)
    
    startImageIndex = np.int(np.array(np.nonzero(Track.df['Image name'][Track.imageIndex]== startImage_str)))
    stopImageIndex = np.int(np.array(np.nonzero(Track.df['Image name'][Track.imageIndex] == stopImage_str)))
    
    

    print(startImageIndex)
    print(stopImageIndex)
#        
#        print(Track1.ImageName[imageIndex[startImageIndex]])
#        print(Track1.ImageName[imageIndex[stopImageIndex]])
    
    # Find the global index corresponding to the image indices
    
    # Contains the global indices in the data corresponding to the images
    
    dT_array = [Track.df['Time'][Track.imageIndex[ii+1]] - Track.df['Time'][Track.imageIndex[ii]] for ii in range(startImageIndex,stopImageIndex)]
        
    ImagePairs = [[Track.df['Image name'][Track.imageIndex[ii]], Track.df['Image name'][Track.imageIndex[ii+1]]]  for ii in range(startImageIndex,stopImageIndex)]
    ImagePairs = np.array(ImagePairs)
    
    print(ImagePairs)
    
    return ImagePairs, dT_array
    


ImagePairs, dT_array = createImagePairArray(Track, startImage=ImgMin, nImages = 10)

# PIV analysis of batch of images

doPIVBatch(Track = Track, ImagePairs = ImagePairs, dT_array = dT_array, save = False)






