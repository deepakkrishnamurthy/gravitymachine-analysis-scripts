#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  1 01:04:22 2019
@author: deepak
"""
import os
import numpy as np
import pandas as pd
import scipy
import scipy.interpolate as interpolate
from scipy.ndimage.filters import uniform_filter1d
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
plt.close("all")
import cv2
import GravityMachine.imageprocessing.imageprocessing_utils as ImageProcessing
from GravityMachine._def import VARIABLE_MAPPING, units


try:
    import rangeslider_functions
except:
    print('Range slider functions not found')

try:
    import PIVanalysis.PIV_Functions as PIV_Functions
except:
    print('PIV modules not found')


class GravityMachineAnalysis:

    def __init__(self, track_folder = None, track_file = None, organism = 'Organism', condition = 'Condition', Tmin=0, Tmax=0, 
        image_streams = None, pixelPermm = None, flip_z = False, rescale_time = True):
        """ Gravity Machine Analysis object for loading, interacting, analyzing and plotting Gravity Machine data.

        """
        self.data = {}
        self.Organism = organism
        self.Condition = condition

        self.Tmin = Tmin
        self.Tmax = Tmax
        
        self.track_file = track_file  # Full, absolute path to the track csv file being analyzed
        self.track_folder = track_folder

        self.image_streams = image_streams  # Image streams whose images are used in the analysis. This can be a list for mutiple image streams. 

        self.pixelPermm = None
      
        self.load_metadata()
     
        # Read the CSV file as a pandas dataframe
        self.df = pd.read_csv(os.path.join(self.track_folder, self.track_file))

        if(rescale_time == True):
            self.df[VARIABLE_MAPPING['Time']] = self.df[VARIABLE_MAPPING['Time']] - self.df[VARIABLE_MAPPING['Time']][0]
        
        if(Tmax == 0 or Tmax is None):
            Tmax = np.max(self.df[VARIABLE_MAPPING['Time']])
        
        if(Tmin is not None and Tmax is not None):
            # Crop the trajectory
            self.df = self.df.loc[(self.df[VARIABLE_MAPPING['Time']]>=Tmin) & (self.df[VARIABLE_MAPPING['Time']] <= Tmax)]
        
        # Internal representation of the track data
        for key in VARIABLE_MAPPING:
            self.data[key] = np.array(self.df[VARIABLE_MAPPING[key]])
            
        self.trackLen = len(self.data)

        self.regularize_time_intervals()

    def set_data(self, data):
        """ Alternative means of loading data from a user supplied data-frame.
            data must be in the form of a data-frame

        """
        # data must be a dataframe. 

        assert(isinstance(data, pd.DataFrame) == True)

        # The data frame should contain the expected data. 
        for key in VARIABLE_MAPPING:
            assert key in data.columns

        self.data = data

    def regularize_time_intervals(self):
        """ Resample the data on a regular time grid

        """
        # Find the average sampling frequency
        self.T = np.linspace(self.data['Time'][0],self.data['Time'][self.trackLen-1], self.trackLen)  # Create a equi-spaced (in time) vector for the data.

        #Sampling Interval
        self.dT = self.T[1]-self.T[0]
        self.samplingFreq = 1/float(self.dT)

        if(self.T is not None):
            func_X = interpolate.interp1d(self.data['Time'], self.data['X'], kind = 'linear')
            func_Y = interpolate.interp1d(self.data['Time'],self.data['Y'], kind = 'linear')
            func_Z = interpolate.interp1d(self.data['Time'],self.data['Z'], kind = 'linear')

            self.X = func_X(self.T)
            self.Y = func_Y(self.T)
            self.Z = func_Z(self.T)
            
    def load_metadata(self):
        """ Load track metadata if available.
        
        """
        metadata = {'Local time': None, 'PixelPermm': None, 'Objective': None}

        if(self.track_folder is not None):
            
            metadata_file = os.path.join(self.track_folder, 'metadata.csv')
            
            if(os.path.exists(metadata_file)):
                print(50*'*')
                print('Loading metadata file ...')
                metadata_df = pd.read_csv(metadata_file)
                
                for key in metadata:
                    if(key in metadata_df.columns):
                        metadata[key] = metadata_df[key][0]

            
            print('Loaded metadata...')
            print(metadata)
            print(50*'*')
    
    def build_image_dict(self):
        """
            Builds a dictionary with image names for each image stream so the image can be loaded from the appropriate sub-directory just based in the image name.
        """
        if(self.image_streams is not None):
            self.image_dicts = {key:{} for key in self.image_streams} # Create a nested dictionary to store the paths for image stream

         
            if os.path.exists(self.track_folder):

                for image_stream_name in self.image_streams:
                # Walk through the folders and identify ones that contain images
                    for dirs, subdirs, files in os.walk(self.track_folder, topdown=False):
                       
                        root, subFolderName = os.path.split(dirs)
                        print('Dir: ',dirs)
                        print('Subdir: ',subdirs)
                        for file in files: 
                            # Find folders that contain the image stream we want to load. For backward compatibility this also checks all folders starting with 'images' 
                            if(file.lower().endswith('tif') and (((image_stream_name in dirs) or ('images' in dirs)))):
                                key = file
                                value = dirs
                                self.image_dict[image_stream_name][key] = value   
        else:
            prinr('No image streams found')
            
          


    def saveAnalysisData(self, overwrite = True):

        if(overwrite or os.path.exists(self.analysis_save_path)==False):
            self.df_analysis = pd.DataFrame({'Organism':[],'Condition':[],'Size':[],'Local time':[],'Track description':[],'Time':[], 'Image name':[], 'Xpos_raw':[],'Zpos_raw':[],'Xobj':[],'Yobj':[],'ZobjWheel':[],'Xvel':[],'Yvel':[],'Zvel':[]})
            
            analysis_len = len(self.imageIndex_array)
            
            print('Local time {}'.format(self.localTime))
            print('Track description {}'.format(self.track_desc))
            
            try:
                self.df_analysis = self.df_analysis.append(pd.DataFrame({'Organism':np.repeat(self.Organism,analysis_len,axis = 0),
                                 'Condition':np.repeat(self.Condition,analysis_len,axis = 0),
                                 'Size': np.repeat(self.OrgDim,analysis_len,axis = 0),
                                 'Local time':np.repeat(self.localTime,analysis_len, axis = 0),
                                 'Track description':np.repeat(self.track_desc, analysis_len, axis=0),
                                 'Time':self.data['Time'][self.imageIndex_array], 
                                 'Image name':self.data['Image name'][self.imageIndex_array], 
                                 'Xpos_raw':self.data['Xobj'][self.imageIndex_array],
                                 'Zpos_raw':self.data['ZobjWheel'][self.imageIndex_array],
                                 'Xobj':self.X_objFluid,'Yobj':self.data['Yobj'][self.imageIndex_array], 
                                 'ZobjWheel':self.Z_objFluid,'Xvel':self.Vx_objFluid,
                                 'Yvel':self.Vy[self.imageIndex_array],'Zvel':self.Vz_objFluid}))
            except:
                self.df_analysis = self.df_analysis.append(pd.DataFrame({'Organism':np.repeat(self.Organism,analysis_len,axis = 0),
                                 'Condition':np.repeat(self.Condition,analysis_len,axis = 0),
                                 'Size': np.repeat(self.OrgDim,analysis_len,axis = 0),
                                 'Local time':np.repeat(self.localTime,analysis_len, axis = 0),
                                 'Track description':np.repeat(self.track_desc, analysis_len, axis=0),
                                 'Time':self.data['Time'][self.imageIndex_array], 
                                 'Image name':self.data['Image name'][self.imageIndex_array], 
                                 'Xpos_raw':self.data['Xobj'][self.imageIndex_array],
                                 'Zpos_raw':self.data['ZobjWheel'][self.imageIndex_array],
                                 'Xobj':self.X_objFluid,'Yobj':self.data['Yobj'][self.imageIndex_array], 
                                 'ZobjWheel':self.Z_objFluid,'Xvel':self.Vx[self.imageIndex_array],
                                 'Yvel':self.Vy[self.imageIndex_array],'Zvel':self.Vz_objFluid}))
                
            self.df_analysis.to_csv(self.analysis_save_path)

        
    def loadAnalysisData(self):

        if(os.path.exists(self.analysis_save_path)):

            self.data =  pd.read_csv(self.analysis_save_path)

        else:

            print('Analysis data does not exist!')
            
    
    def createImageIndex(self):
        # Create an index of all time points for which an image is available
        self.imageIndex = []
        for ii in range(self.trackLen):
            
            if(self.data['Image name'][ii] is not np.nan):
#                print(self.data['Image name'][ii])
                self.imageIndex.append(ii)
                
#        print(self.imageIndex)
                
        # Open the first image and save the image size
        imageName = self.data['Image name'][self.imageIndex[0]]
        
        image_a = cv2.imread(os.path.join(self.path,self.image_dict[imageName],imageName))
        
      
        self.imH, self.imW, *rest = np.shape(image_a)
        
        
        
    def setColorThresholds(self, overwrite = False):
        '''
        Displays an image and allows the user to choose the threshold values so the object of interest is selected
        '''
        
        saveFile = 'colorThresholds.pkl'

        image_a = None
        
        if(not os.path.exists(os.path.join(self.root, saveFile)) or overwrite):
            # If a color threshold does not exist on file then display an image and allow the user to choose the thresholds
            
            
        
            imageName = self.data['Image name'][self.imageIndex[0]]
            
            image_a = cv2.imread(os.path.join(self.path,self.image_dict[imageName],imageName))
            
          
            self.imH, self.imW, *rest = np.shape(image_a)
    
            
        
            print('Image Width: {} px \n Image Height: {} px'.format(self.imW, self.imH))
            
            print(os.path.join(self.path, self.image_dict[imageName],imageName))
            v1_min,v2_min,v3_min,v1_max,v2_max,v3_max = rangeslider_functions.getColorThreshold(os.path.join(self.path,self.image_dict[imageName],imageName))
            threshLow = (v1_min,v2_min,v3_min)
            threshHigh = (v1_max,v2_max,v3_max)
            
            # Save this threshold to file
            with open(os.path.join(self.root,saveFile), 'wb') as f:  # Python 3: open(..., 'wb')
                pickle.dump((threshLow, threshHigh), f)
        else:
            # If available just load the threshold
            print('Color thresholds available! \n Loading file {} ...'.format(os.path.join(self.root,saveFile)))
            with open(os.path.join(self.root,saveFile), 'rb') as f:
                threshLow, threshHigh = pickle.load(f)
                
        self.threshLow = threshLow
        self.threshHigh = threshHigh
        
        print('Color thresholds for segmentation: \n LOW: {}, HIGH : {}'.format(self.threshLow, self.threshHigh))
                
            
    def findOrgDims(self, circle=0, overwrite = False):
        # Finds the maximum dimensions of the organism 
        saveFile = 'orgDims.csv'

        OrgMajDim = []
        OrgMinDim = []
        OrgDim = []

        size_df = pd.DataFrame({'Organism':[],'Condition':[],'Track':[],'OrgDim':[],'OrgMajDim':[],'OrgMinDim':[]})
    
        print(self.path)
        if(not os.path.exists(os.path.join(self.path,saveFile)) or overwrite):
            
            nImages = 100
            nTotal = len(self.imageIndex)

            print(nTotal)
         

            fileList = self.data['Image name'][self.imageIndex[:nImages]]
            print(fileList)
            print(type(fileList))
            # Calculate based on 100 images

            # Enter an event loop where the use can choose which image is used for calculating organism size
            img_num = 0
            roiFlag = False

            while True:
                file = fileList.iloc[img_num]
                print(file)
                image = cv2.imread(os.path.join(self.path,self.image_dict[file],file))

                if(roiFlag is True):
                    r = cv2.selectROI('Select ROI', image)
                    image = image[int(r[1]):int(r[1]+r[3]), int(r[0]):int(r[0]+r[2])]
                    roiFlag = False
                    cv2.destroyWindow('Select ROI')


                orgContour = ImageProcessing.colorThreshold(image = image, threshLow = self.threshLow, threshHigh = self.threshHigh)

                if(orgContour is not None):
                    
                    if(circle):
                        (x_center,y_center), Radius = cv2.minEnclosingCircle(orgContour)
                        center = (int(x_center), int(y_center))
                        cv2.circle(image,center, int(Radius),(0,255,0),2)
                        cv2.imshow('Press ESC to exit, SPACEBAR for next img, R for ROI',image)
                        key = cv2.waitKey(0)
                    else:
                        ellipse = cv2.fitEllipse(orgContour)
                        cv2.ellipse(image,box=ellipse,color=[0,1,0])
                        cv2.imshow('Press ESC to exit, SPACEBAR for next img, R for ROI',image)
                        key = cv2.waitKey(0)
            #                 
                    
                    if(key == 27):
                        # ESC. If the organism is detected correctly, then store the value and break from the loop
                        if(circle):
                            OrgMajDim.append(self.mmPerPixel*2*Radius)
                            OrgMinDim.append(self.mmPerPixel*2*Radius)
                            OrgDim.append(self.mmPerPixel*2*Radius)

                        else:
                            OrgMajDim.append(self.mmPerPixel*ellipse[1][0])
                            OrgMinDim.append(self.mmPerPixel*ellipse[1][1])
                            OrgDim.append(self.mmPerPixel*(ellipse[1][1] + ellipse[1][0])/float(2))

                        cv2.destroyAllWindows()
                        break
                    elif(key == 32):
                        # Spacebar: If Organism is not found in given frame then show the next frame
                        img_num += 1
                        if(img_num >= nImages):
                            img_num = 0
                        continue
                    elif(key == ord('r')):
                        # Press 'r'. If the object is present but is not the only bright object in the frame
                        roiFlag = True
                    elif(key == ord('c')):
                        self.setColorThresholds(overwrite = True)
                        
                        
                else:
                    # Select new color thresholds
                    img_num += 1
                    key = ord('c')
                        



            OrgDim_mean = np.nanmean(np.array(OrgDim))
            OrgMajDim_mean = np.nanmean(np.array(OrgMajDim))
            OrgMinDim_mean = np.nanmean(np.array(OrgMinDim))

            self.OrgDim = OrgDim_mean
            self.OrgMajDim = OrgMajDim_mean
            self.OrgMinDim = OrgMinDim_mean
            
            # with open(os.path.join(self.path,saveFile), 'wb') as f:  # Python 3: open(..., 'wb')
            #     pickle.dump((OrgDim_mean, OrgMajDim_mean, OrgMinDim_mean), f)
            # Save the Organism dimensions to file

            size_df = size_df.append(pd.DataFrame({'Organism':[self.Organism],'Condition':[self.Condition],'Track':[self.trackFolder],'OrgDim':[self.OrgDim],'OrgMajDim':[self.OrgMajDim],'OrgMinDim':[self.OrgMinDim]}))
                        
            size_df.to_csv(os.path.join(self.path, saveFile))

        else:
            # Load the Organism Size data
            print('Loading organism size from memory ...')
            
            # with open(os.path.join(self.path,saveFile), 'rb') as f:
            #     OrgDim_mean, OrgMajDim_mean, OrgMinDim_mean = pickle.load(f)
            size_df = pd.read_csv(os.path.join(self.path,saveFile))
        
            self.OrgDim = size_df['OrgDim'][0]
            self.OrgMajDim = size_df['OrgMajDim'][0]
            self.OrgMinDim = size_df['OrgMinDim'][0]
            
        print('*'*50)
        print('Organism dimension {} mm'.format(self.OrgDim))
        print('Organism Major dimension {} mm'.format(self.OrgMajDim))
        print('Organism Minor dimension {} mm'.format(self.OrgMinDim))
        print('*'*50)
                        
    def velocity_central_diff(self):

        for i in range(0,self.trackLen):
                
                
            if(i==0):
                # Forward difference at the start points
                self.Vx[i] = (self.data[self.Xobj_name][i+1]-self.data[self.Xobj_name][i])/(self.data['Time'][i+1]-self.data['Time'][i])
                self.Vy[i] = (self.data[self.Yobj_name][i+1]-self.data[self.Yobj_name][i])/(self.data['Time'][i+1]-self.data['Time'][i])
                self.Vz[i] = (self.data['ZobjWheel'][i+1]-self.data['ZobjWheel'][i])/(self.data['Time'][i+1]-self.data['Time'][i])
                
                self.Vz_objLab[i] = (self.data[self.Zobj_name][i+1]-self.data[self.Zobj_name][i])/(self.data['Time'][i+1]-self.data['Time'][i])
                
                self.Theta_dot[i] = (self.data['ThetaWheel'][i+1]-self.data['ThetaWheel'][i])/(self.data['Time'][i+1]-self.data['Time'][i])

                if self.XposImageAvailable:
                    # Note: This will be Vx_objLab for a setup where the optical FOV is fixed in the lab reference
                    # > GM v2.0 and higher
#                    self.Vx_objStage[i] = (self.data[self.XobjImage_name][i+1]-self.data[self.XobjImage_name][i])/(self.data['Time'][i+1]-self.data['Time'][i])
                    self.Vx_objLab[i] = (self.data[self.XobjImage_name][i+1]-self.data[self.XobjImage_name][i])/(self.data['Time'][i+1]-self.data['Time'][i])

            elif(i==self.trackLen-1):
                # Backward difference at the end points
                self.Vx[i] = (self.data[self.Xobj_name][i]-self.data[self.Xobj_name][i-1])/(self.data['Time'][i]-self.data['Time'][i-1])
                self.Vy[i] = (self.data[self.Yobj_name][i]-self.data[self.Yobj_name][i-1])/(self.data['Time'][i]-self.data['Time'][i-1])
                self.Vz[i] = (self.data['ZobjWheel'][i]-self.data['ZobjWheel'][i-1])/(self.data['Time'][i]-self.data['Time'][i-1])
                
                self.Vz_objLab[i] = (self.data[self.Zobj_name][i]-self.data[self.Zobj_name][i-1])/(self.data['Time'][i]-self.data['Time'][i-1])
                
                self.Theta_dot[i] = (self.data['ThetaWheel'][i]-self.data['ThetaWheel'][i-1])/(self.data['Time'][i]-self.data['Time'][i-1])

                if self.XposImageAvailable:
#                    self.Vx_objStage[i] = (self.data[self.XobjImage_name][i]-self.data[self.XobjImage_name][i-1])/(self.data['Time'][i]-self.data['Time'][i-1])
                    self.Vx_objLab[i] = (self.data[self.XobjImage_name][i]-self.data[self.XobjImage_name][i-1])/(self.data['Time'][i]-self.data['Time'][i-1])

            else:
                # Central difference for all other points
                self.Vx[i] = (self.data[self.Xobj_name][i+1]-self.data[self.Xobj_name][i-1])/(self.data['Time'][i+1]-self.data['Time'][i-1])
                self.Vy[i] = (self.data[self.Yobj_name][i+1]-self.data[self.Yobj_name][i-1])/(self.data['Time'][i+1]-self.data['Time'][i-1])
                self.Vz[i] = (self.data['ZobjWheel'][i+1]-self.data['ZobjWheel'][i-1])/(self.data['Time'][i+1]-self.data['Time'][i-1])
                
                self.Vz_objLab[i] = (self.data[self.Zobj_name][i+1]-self.data[self.Zobj_name][i-1])/(self.data['Time'][i+1]-self.data['Time'][i-1])
                
                
                self.Theta_dot[i] = (self.data['ThetaWheel'][i+1]-self.data['ThetaWheel'][i-1])/(self.data['Time'][i+1]-self.data['Time'][i-1])

                if self.XposImageAvailable:
#                    self.Vx_objStage[i] = (self.data[self.XobjImage_name][i + 1]-self.data[self.XobjImage_name][i-1])/(self.data['Time'][i+1]-self.data['Time'][i-1])
                    self.Vx_objLab[i] = (self.data[self.XobjImage_name][i + 1]-self.data[self.XobjImage_name][i-1])/(self.data['Time'][i+1]-self.data['Time'][i-1])




    def computeVelocity(self):
        
        self.Vx = np.zeros(self.trackLen)
        self.Vy = np.zeros(self.trackLen)
        self.Vz = np.zeros(self.trackLen)
        self.Vz_objLab = np.zeros(self.trackLen)
        # Velocity of the object in the reference frame of the X-stage (GM v<2.0)
#        self.Vx_objStage = np.zeros(self.trackLen)

        # Velocity of the object in reference frame of the fixed optical FOV (GM v>2.0)
        self.Vx_objLab = np.zeros(self.trackLen)
        self.Theta_dot = np.zeros(self.trackLen)
        
        # If using post-procssed data, then load velocities from file
        if(self.use_postprocessed):
            self.Vx = self.data['Xvel']
            self.Vy = self.data['Yvel']
            self.Vz = self.data['Zvel']
            
            
        else:
            # Try to calculate the velocities
            self.velocity_central_diff()
            
        # Smooth the velocity data to only keep frequencies 10 times lower than the sampling frequency (low-pass filter)
        self.Vx_smooth = self.smoothSignal(self.Vx, window_time = self.window_time)
        self.Vy_smooth = self.smoothSignal(self.Vy, window_time = self.window_time)
        self.Vz_smooth = self.smoothSignal(self.Vz, window_time = self.window_time)
        self.Vz_objLab_smooth = self.smoothSignal(self.Vz_objLab, window_time = self.window_time)
        self.Theta_dot_smooth = self.smoothSignal(self.Theta_dot, window_time = self.window_time)
        self.Speed = (self.Vx_smooth**2 + self.Vy_smooth**2 + self.Vz_smooth**2)**(1/2)
        self.Speed_z = (self.Vz_smooth**2)**(1/2)
        
    
    def computeAccln(self):
        
        self.Theta_ddot = np.zeros(self.trackLen)
        
        self.a_z = np.zeros(self.trackLen)
        
        for i in range(0,self.trackLen):
            
            if(i==0):
                # Forward difference at the start points
                
                self.a_z[i] = (self.Vz[i+1]-self.Vz[i])/(self.data['Time'][i+1]-self.data['Time'][i])
                self.Theta_ddot[i] = (self.Theta_dot[i+1]-self.Theta_dot[i])/(self.data['Time'][i+1]-self.data['Time'][i])

            
            elif(i==self.trackLen-1):
                # Backward difference at the end points
                self.a_z[i] = (self.Vz[i]-self.Vz[i-1])/(self.data['Time'][i]-self.data['Time'][i-1])
                self.Theta_ddot[i] = (self.Theta_dot[i]-self.Theta_dot[i-1])/(self.data['Time'][i]-self.data['Time'][i-1])
            else:
                # Central difference for all other points
                self.a_z[i] = (self.Vz[i+1]-self.Vz[i-1])/(self.data['Time'][i+1]-self.data['Time'][i-1])
                self.Theta_ddot[i] = (self.Theta_dot[i+1]-self.Theta_dot[i-1])/(self.data['Time'][i+1]-self.data['Time'][i-1])
    
        self.a_z = self.smoothSignal(self.a_z, self.window_time)
        self.Theta_ddot = self.smoothSignal(self.Theta_ddot, self.window_time)

    def computeVelocityOrientation(self):
        
        
        vector_magnitude = self.Speed
        
        Orientation_vectors = np.zeros((3, len(self.Vx_smooth)))
        
        Orientation_vectors[0,:] = self.Vx_smooth/vector_magnitude
        Orientation_vectors[1,:] = self.Vy_smooth/vector_magnitude
        Orientation_vectors[2,:] = self.Vz_smooth/vector_magnitude
        
        # Extract the orientation angle from the vertical
        
        Z_gravity = [0, 0, 1]
        
        cos_theta = Orientation_vectors[0,:]*Z_gravity[0] + Orientation_vectors[1,:]*Z_gravity[1] + Orientation_vectors[2,:]*Z_gravity[2]

        # Theta value in degrees
        self.orientation_theta = np.arccos(cos_theta)*(180/(np.pi))
        
        
        
    def computeDisplacement(self, x_data = None, y_data = None):
        '''
        Compute the displacement by integrating a velocity time series.
        '''
        
        disp = scipy.integrate.cumtrapz(y = y_data, x = x_data, initial = 0)
        
        return disp
#        self.disp_z_computed = scipy.integrate.cumtrapz(y = self.Vz, x = self.data['Time'])
        
        

    def computeFluidVelocity(self, image_a, image_b, deltaT = 1, overwrite_piv = False, overwrite_velocity = False, masking = False, obj_position = None, obj_size = 0.1):
        '''
        Computes the mean fluid velocity given a pair of images, 
        far away from objects (if any objects are present).
        
        '''        
    
        #--------------------------------------------------------------------------
        # Load the frame-pair into memory
        #--------------------------------------------------------------------------
        frame_a_color = cv2.imread(os.path.join(self.path, self.image_dict[image_a], image_a))
        frame_b_color = cv2.imread(os.path.join(self.path, self.image_dict[image_b], image_b))
    
        
       
        
        # Plot the object's position on the image to verify it is correct
#        print('Circle diameter: {}'.format(int(2*self.scaleFactor*self.OrgDim*self.pixelPermm)))
#        frame_a_color_copy = np.copy(frame_a_color)
#        
#        cv2.circle(frame_a_color_copy, (int(obj_position[0]), int(obj_position[1])), int(self.OrgDim*self.pixelPermm*self.scaleFactor), [255,255,255])
#        
#        cv2.imshow('Frame', frame_a_color_copy)
#        cv2.waitKey(1)

        #--------------------------------------------------------------------------
        # Perform the PIV computation
        #--------------------------------------------------------------------------
        saveFile = os.path.join(self.PIVfolder,'PIV_' + image_a[:-4]+'.pkl')
        
        
        if(not os.path.exists(saveFile) or overwrite_piv):
            print('-'*50)
            print('Analyzing Frame pairs: {} and {} \n'.format(image_a,image_b))
            print('-'*50)
            x,y,u,v, sig2noise = PIV_Functions.doPIV(frame_a_color,frame_b_color, dT = deltaT, win_size = self.window_size, overlap = self.overlap, searchArea = self.searchArea, apply_clahe = False)
            
            
            
            u, v = PIV_Functions.pivPostProcess(u,v,sig2noise, sig2noise_min = 1.0, smoothing_param = 0)
            
            
            
            u,v = (PIV_Functions.data2RealUnits(data = u,scale = 1/(self.pixelPermm)), PIV_Functions.data2RealUnits(data = v,scale = 1/(self.pixelPermm)))
            
            #--------------------------------------------------------------------------
            # Threshold the image to extract the object regions
            #--------------------------------------------------------------------------
            if(masking and obj_position is None):
                Contours = PIV_Functions.findContours(frame_a_color,self.threshLow,self.threshHigh,'largest')
            else:
                Contours = np.nan
            
            
            with open(saveFile, 'wb') as f:  # Python 3: open(..., 'wb')
                pickle.dump((x, y , u, v, Contours), f)
            
        else:
            #--------------------------------------------------------------------------
            # Read the PIV data 
            #--------------------------------------------------------------------------
            print('-'*50)
            print('Loading PIV data for: {} and {} \n'.format(image_a,image_b))
            print('-'*50)
            pklFile = saveFile
            x,y,u,v,Contours = PIV_Functions.readPIVdata(pklFile)
            print(np.nanmean(u+v))
           
        
        
        
        
        
#        Centroids = PIV_Functions.findCentroids(Contours)
        
       
        
#        plt.figure()
#        plt.imshow(maskInsideCircle)
#        plt.show()
        
#        cv2.circle(frame_a_gs,(int(x_cent), int(y_cent)),int(scale_factor*radius), color = (0,255,0))
#        cv2.imshow('frame',frame_a_gs)
#        cv2.waitKey(1)
        
#        cv2.drawContours(frame_a_color, [Contours],0,(0,255,0),3)
#        cv2.imshow('frame',frame_a_color)
#        cv2.waitKey(1)

        
        
#        PIV_Functions.plotPIVdata(frame_a_color,x,y,u,v, orgContour=Contours)
        
        
        if(masking is True):    
            if(obj_position is None):
                x_cent, y_cent, radius = PIV_Functions.findCircularContour(Contours)
            else:
                x_cent, y_cent = obj_position
                radius = round((self.OrgDim/2.0)*self.pixelPermm)
                
            maskInsideCircle = PIV_Functions.pointInCircle(x,y,x_cent,y_cent,self.scaleFactor*radius)

            u[maskInsideCircle] = np.nan
            v[maskInsideCircle] = np.nan 
                
            
#            assert(np.nanmean(u+v) is not np.nan)
         
            
            
            
#            print(x_cent, y_cent)
#            print(radius)
#            
#            print(maskInsideCircle)
            
#            plt.figure(1)
##            plt.scatter(x_cent, y_cent, 'ro')
#            plt.imshow(maskInsideCircle)
#            plt.pause(0.001)
#            plt.show()
            
            
#            plt.figure(2)
#            cv2.circle(frame_a_gs,(int(x_cent), int(y_cent)),int(scaleFactor*radius), color = (0,255,0))
#            cv2.imshow('frame',frame_a_gs)
#            cv2.waitKey(1)
       
        else:
            pass
          
           # Plot the vectors
        PIV_Functions.plotPIVdata(frame_a_color, x, y, u, v, Centroids = obj_position)
        # Find the mean velocity
        u_avg, v_avg = (np.nanmean(u), np.nanmean(v))
        u_std, v_std = (np.nanstd(u), np.nanstd(v))
        
        print(u_avg)
        print(v_avg)
        return u_avg, v_avg, u_std, v_std
 
      
    def FluidVelTimeSeries(self, overwrite_velocity = False, masking = False):
        # 
        '''
        Computes the Fluid velocity at each time point during which an image is available 
        and stores the result.
        
        For each pair of time-points with images:
            > Generate a mask for each image as necessary
            > Do PIV on the pair of images
            > Calculate the average velocity of the fluid in regions far from the object
            > Store the average velocity as a finction of time
        '''
        if(not os.path.exists(self.FluidVelocitySavePath) or overwrite_velocity):
            
            print("calculating fluid velocity time series ...")
        
            nImages = len(self.imageIndex)
#            nImages = 100
            
            n = min(nImages, len(self.imageIndex)-1)
            
            self.u_avg_array = np.zeros(n)
            self.v_avg_array = np.zeros(n)
            self.u_std_array = np.zeros(n)
            self.v_std_array = np.zeros(n)
            self.imageIndex_array = np.zeros(n, dtype='int')
            
            
            for ii in range(n):
                
               
                                    
                imageindex_a = self.imageIndex[ii] 
                imageindex_b = self.imageIndex[ii + 1]
                
                image_a = self.data['Image name'][imageindex_a]
                image_b = self.data['Image name'][imageindex_b]
                
                try:
                    obj_position = (self.imW/2 - round(self.data['Xobj_image'][imageindex_a]*self.pixelPermm), self.imH/2 - round(self.data['Zobj'][imageindex_a]*self.pixelPermm))
                except:
                    obj_position = (self.imW/2, self.imH/2 - round(self.data['Zobj'][imageindex_a]*self.pixelPermm))
                
                
                # First check if both these images exist in memory
                try:
                    image_a_exists = os.path.exists(os.path.join(self.path, self.image_dict[image_a], image_a))
                except:
                    image_a_exists = False
                try:
                    image_b_exists = os.path.exists(os.path.join(self.path, self.image_dict[image_b], image_b))
                except:
                    image_b_exists = False
                                
                image_a_num = int(image_a[4:-4])
                image_b_num = int(image_b[4:-4])
                
                frame_gap = image_b_num - image_a_num
                
                print(frame_gap)
                
                if(image_a_exists and image_b_exists and frame_gap == 1):
                    print('Consequtive images found ...')
                
                    print(image_a)
                    print(image_b)
                    
                    dT = self.data['Time'][imageindex_b] - self.data['Time'][imageindex_a]
                    
                    self.u_avg_array[ii], self.v_avg_array[ii], self.u_std_array[ii], self.v_std_array[ii] = self.computeFluidVelocity(image_a,image_b,deltaT = dT, masking = masking, obj_position = obj_position, obj_size = self.OrgDim, overwrite_piv = self.overwrite_piv)
                    self.imageIndex_array[ii] = imageindex_a
                    
                # If either of those images do not exist, assume that the velocity remains constant over the missing frames
                elif(not image_a_exists or not image_b_exists):
                    print('One or more of image pair not found...')
                    print('Checking for next image index...')
                    self.u_avg_array[ii], self.v_avg_array[ii], self.u_std_array[ii], self.v_std_array[ii] = self.u_avg_array[ii-1], self.v_avg_array[ii-1], self.u_std_array[ii-1], self.v_std_array[ii-1] 
                    self.imageIndex_array[ii] = imageindex_a
                    continue
               
            self.u_avg_array, self.v_avg_array = (self.smoothSignal(self.u_avg_array, self.window_time),self.smoothSignal(self.v_avg_array, self.window_time))
            
            with open(self.FluidVelocitySavePath, 'wb') as f:  # Python 3: open(..., 'wb')
                    pickle.dump((self.imageIndex_array, self.u_avg_array, self.v_avg_array, self.u_std_array, self.v_std_array), f)
                
            
        else:
            print("Fluid time series found! Loading ...")
            with open(self.FluidVelocitySavePath, 'rb') as f:  # Python 3: open(..., 'wb')
                    self.imageIndex_array, self.u_avg_array, self.v_avg_array, self.u_std_array, self.v_std_array = pickle.load(f)
                
                
  
    def correctedDispVelocity(self, overwrite_flag = False):
        
        
        self.FluidVelTimeSeries(overwrite_velocity = overwrite_flag)
        '''
             Vector operations to calculate V_objFluid which is what we want.
             
             Note: Vz is actually the object velocity relative to the stage V_objStage
             V_objLab: is measured from the displacement of the object centroid in the image. 
             V_objStage = V_objLab - V_stageLab
             Therefore, 
                 V_stageLab = V_objLab - VobjStage   ---- (1)
             
             The measured fluid velocity using PIV is:
                 V_measured = V_stageLab + V_fluidStage
                 Therefore, 
                 V_fluidStage = V_measured - V_stageLab ---- (2)
                 We can substitute for V_stageLab from (1) to get V_fluidStage
                 
             Now, 
                 V_objFluid = V_objStage - V_fluidStage
                 
                
        '''
        
        Vz_stageLab = -self.Vz[self.imageIndex_array] + self.Vz_objLab[self.imageIndex_array]
        
        Vz_fluidStage = self.v_avg_array - Vz_stageLab
        
        self.Vz_objFluid = self.Vz[self.imageIndex_array] - Vz_fluidStage
        
        self.Z_objFluid =  self.computeDisplacement(x_data = self.data['Time'][self.imageIndex_array], 
                                                    y_data = self.Vz_objFluid)
      
        # Correcting X-velocity
        
        #----------------------------------------------------------------------------------------
        # For GM v > 2.0 data (Fixed Optical System, Moving Stage in X,Y, Theta)
        
        # Note that if the X-centroid of the object is not available then the velocity contribution of V_objLab 
        # is assumed to be zero.
        Vx_stageLab = -self.Vx[self.imageIndex_array] + self.Vx_objLab[self.imageIndex_array]
    
        Vx_fluidStage = self.u_avg_array - Vx_stageLab
    
        self.Vx_objFluid = self.Vx[self.imageIndex_array] - Vx_fluidStage
    
        self.X_objFluid =  self.computeDisplacement(x_data = self.data['Time'][self.imageIndex_array], 
                                                y_data = self.Vx_objFluid)
        #----------------------------------------------------------------------------------------
        # Uncomment block below for GM < v2.0 data
#            self.Vx_objFluid = self.Vx_objStage[self.imageIndex_array] - self.u_avg_array
#        
#            self.X_objFluid =  self.computeDisplacement(x_data = self.data['Time'][self.imageIndex_array], 
#                                                        y_data = self.Vx_objFluid)  
        #----------------------------------------------------------------------------------------
       
            


    #--------------------------------------------------------------------------       
    # Signal Processing Functions
    #--------------------------------------------------------------------------  
    # @@@ Need to modift and use the pandas based method @@@     
    def smoothSignal(self, data, window_time):      # Window is given in seconds
            
            avgWindow = int(window_time*self.samplingFreq)
            return uniform_filter1d(data, size = avgWindow, mode="reflect")
#            data = pd.Series(data)
#            rolling_mean = np.array(data.rolling(window = avgWindow, center = True).mean())
##            try:
#            return rolling_mean
#            except:
#                return pd.rolling_mean(data, avgWindow, min_periods = 1, center = True)
          
    #--------------------------------------------------------------------------
    # Plotting functions
    #--------------------------------------------------------------------------

    def plot_displacement_timeseries(self, save = False):

        title = 'Displacement time series'
        f, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True, figsize = (8,12))
        ax1.set_title(title)
        sns.lineplot(x = self.data['Time'], y = self.data['X'], color = 'r', linewidth = 1, label = 'X', ax = ax1, ci = None)
        ax1.set_ylabel('X '+units['X'])
        sns.lineplot(x = self.data['Time'], y = self.data['Y'], color = 'g', linewidth = 1, label = 'Y', ax = ax2, ci = None)
        ax2.set_ylabel('Y '+units['Y'])
        sns.lineplot(x = self.data['Time'], y = self.data['Z'], color = 'b', linewidth = 1, label = 'Z', ax = ax3, ci = None)
        ax3.set_ylabel('Z '+units['Z'])
        ax3.set_xlabel('Time' + units['Time'])
        plt.show()



