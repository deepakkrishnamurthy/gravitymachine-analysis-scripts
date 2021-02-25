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
import GravityMachine.utils as utils
import imp
imp.reload(utils)
from tqdm import tqdm

try:
	import GravityMachine.rangeslider_functions as rangeslider_functions 
except:
	print('Range slider functions not found')

try:
	import GravityMachine.PIVanalysis.PIV_Functions as PIV_Functions
except:
	print('PIV modules not found')

units = {'Time':'(s)', 'X':'(mm)','Y':'(mm)','Z':'(mm)', 'V_x': '(mm/s)','V_y': '(mm/s)','V_z': '(mm/s)', 'Theta':'(rad)'}
imgFormats = ['.png', '.svg']
# Map between the header names in the CSV file and the internal variable names
# X_objStage  Y_objStage  Z_objStage  Theta_stage X_image Z_image
# VARIABLE_MAPPING = {'Time':'Time', 'X':'X_objStage','Y':'Y_objStage','Z':'Z_objStage','Image name':'DF1', 'X_image':'X_image', 'Z_image':'Z_image'}

VARIABLE_MAPPING = {'Time':'Time', 'X':'Xobj','Y':'Yobj','Z':'ZobjWheel','Image name':'Image name', 'X_image':'Xobj_image', 'Z_image':'Zobj'}


class GravityMachineTrack:

	def __init__(self, track_folder = None, track_file = None, organism = 'Organism', condition = 'Condition', Tmin=0, Tmax=0, 
		image_streams = 'Image name', pixelPermm = 314, flip_z = False, rescale_time = True):
		""" Gravity Machine Analysis object for loading, interacting, analyzing and plotting Gravity Machine data.

		"""
		self.data = pd.DataFrame({})
		self.Organism = organism
		self.Condition = condition

		self.Tmin = Tmin
		self.Tmax = Tmax
		
		self.track_file = track_file  # Full, absolute path to the track csv file being analyzed
		self.track_folder = track_folder

		self.image_stream = VARIABLE_MAPPING['Image name']  # Image streams whose images are used in the analysis. This can be a list for mutiple image streams. 

		self.pixelPermm = pixelPermm
	  
		self.load_metadata()
	 
		# Read the CSV file as a pandas dataframe
		self.df = pd.read_csv(os.path.join(self.track_folder, self.track_file))

		self.ColumnNames = list(self.df.columns.values)

		print(self.ColumnNames)
		print(VARIABLE_MAPPING)

		# The data MUSt have timestamps
		assert('Time' in self.ColumnNames)

		if(rescale_time == True):
			self.df[VARIABLE_MAPPING['Time']] = self.df[VARIABLE_MAPPING['Time']] - self.df[VARIABLE_MAPPING['Time']][0]
		
		if(Tmax == 0 or Tmax is None):
			Tmax = np.max(self.df[VARIABLE_MAPPING['Time']])
		
		if(Tmin is not None and Tmax is not None):
			# Crop the trajectory
			self.df = self.df.loc[(self.df[VARIABLE_MAPPING['Time']]>=Tmin) & (self.df[VARIABLE_MAPPING['Time']] <= Tmax)]

		# Internal representation of the track data
		for key in VARIABLE_MAPPING:
			if(VARIABLE_MAPPING[key] in self.ColumnNames):
				self.data[key] = np.array(self.df[VARIABLE_MAPPING[key]])
			else:
				print('Warning {} not found in input data'.format(key))
				self.data[key] = None
			
		self.trackLen = len(self.data['Time'])

		self.regularize_time_intervals()
		self.initialize_data()
		self.build_image_dict()
		self.create_image_index()

		self.plots_folder = os.path.join(self.track_folder, 'Plots')


	def initialize_data(self):
		"""
			Initialize core and derived data.
		"""
		self.data['V_x'] = None
		self.data['V_y'] = None
		self.data['V_z'] = None

		df_fluid_velocity = pd.DataFrame({})

		self.derived_data = {}

		self.threshLow = None
		self.threshHigh = None

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

	def initialize_piv_analysis(self):
		self.PIVfolder = os.path.join(self.track_folder, 'PIV_data')
		if(not os.path.exists(self.PIVfolder)):
			os.makedirs(self.PIVfolder)

		self.piv_settings = {'window size': 64, 'overlap': 32, 'search area':64}

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
			Builds a dictionary with image names for each image stream so the image can be loaded from the appropriate sub-directory just based on the image name.
		"""
		if(self.image_stream is not None):
			self.image_dict = {} 

		 
			if os.path.exists(self.track_folder):

				# Walk through the folders and identify ones that contain images
				for dirs, subdirs, files in os.walk(self.track_folder, topdown=False):
				   
					root, subFolderName = os.path.split(dirs)
					print('Dir: ',dirs)
					print('Subdir: ',subdirs)
					for file in files: 
						# Find folders that contain the image stream we want to load. For backward compatibility this also checks all folders starting with 'images' 
						if(file.lower().endswith('tif') and (((self.image_stream in dirs) or ('images' in dirs)))):
							key = file
							value = dirs
							self.image_dict[key] = value   
		else:
			print('No image streams found')

	def compute_velocity(self, smooth_window = None):
		""" Computes velocity of tracked object using a central-difference scheme. 
			Other methods can also be implemented.

			smooth_window: Sliding time window (in s) over which to average velocity.

		"""
		self.data['V_x'] = utils.velocity_central_diff(self.data['Time'][:], self.data['X'][:])
		self.data['V_y'] = utils.velocity_central_diff(self.data['Time'][:], self.data['Y'][:])
		self.data['V_z'] = utils.velocity_central_diff(self.data['Time'][:], self.data['Z'][:])

		# Object velocity relative to image reference frame
		if(self.data['X_image'][0] is not None):
			self.data['V_x_image'] = utils.velocity_central_diff(self.data['Time'][:], self.data['X_image'][:])
		else:
			# If the object position realtive to image is unvailable then assume it is zero.
			self.data['V_x_image'] = np.zeros(self.trackLen)


		if(self.data['Z_image'][0] is not None):
			self.data['V_z_image'] = utils.velocity_central_diff(self.data['Time'][:], self.data['Z_image'][:])
		else:
			self.data['V_z_image'] = np.zeros(self.trackLen)


		if(smooth_window is not None):

			average_dt = np.nanmean(self.data['Time'][1:] - self.data['Time'][:-1])

			window_size = int(smooth_window/average_dt)

			self.data['V_x'] = np.array(self.data['V_x'].rolling(window = window_size, center = True).mean())


	
	def create_image_index(self):
		""" Create an index of all time points for which an image is available.
			All time points may not have an index available so this provides a convenient indexinf of all time points with an associated image.
		"""
		self.imageIndex = []
		for ii in range(self.trackLen):
			if(self.data['Image name'][ii] is not np.nan):
#                print(self.data['Image name'][ii])
				self.imageIndex.append(ii)
				
		# Open the first image and save the image size
		imageName = self.data['Image name'][self.imageIndex[0]]
		image_a = cv2.imread(os.path.join(self.track_folder,self.image_dict[imageName],imageName))
		self.imH, self.imW, *rest = np.shape(image_a)

		
	def set_color_thresholds(self, overwrite = False):
		"""
		Displays an image and allows the user to choose the threshold values so the object of interest is selected
		"""
		
		saveFile = 'colorThresholds.pkl'

		image_a = None
		
		if(not os.path.exists(os.path.join(self.track_folder, saveFile)) or overwrite):
			# If a color threshold does not exist on file then display an image and allow the user to choose the thresholds
			imageName = self.data['Image name'][self.imageIndex[0]]
			image_a = cv2.imread(os.path.join(self.track_folder,self.image_dict[imageName],imageName))
			self.imH, self.imW, *rest = np.shape(image_a)
	
			print('Image Width: {} px \n Image Height: {} px'.format(self.imW, self.imH))
			
			print(os.path.join(self.track_folder, self.image_dict[imageName],imageName))
			v1_min,v2_min,v3_min,v1_max,v2_max,v3_max = rangeslider_functions.getColorThreshold(os.path.join(self.track_folder,self.image_dict[imageName],imageName))
			threshLow = (v1_min,v2_min,v3_min)
			threshHigh = (v1_max,v2_max,v3_max)
			
			# Save this threshold to file
			with open(os.path.join(self.track_folder, saveFile), 'wb') as f:  # Python 3: open(..., 'wb')
				pickle.dump((threshLow, threshHigh), f)
		else:
			# If available just load the threshold
			print('Color thresholds available! \n Loading file {} ...'.format(os.path.join(self.track_folder, saveFile)))
			with open(os.path.join(self.track_folder, saveFile), 'rb') as f:
				threshLow, threshHigh = pickle.load(f)
				
		self.threshLow = threshLow
		self.threshHigh = threshHigh
		
		print('Color thresholds for segmentation: \n LOW: {}, HIGH : {}'.format(self.threshLow, self.threshHigh))
				
			
	def find_object_size(self, circle = False, overwrite = False):
		""" Finds the dimensions of the tracked object based on a color threshold.

		"""
		saveFile = 'orgDims.csv'

		OrgMajDim = []
		OrgMinDim = []
		OrgDim = []

		size_df = pd.DataFrame({'Organism':[],'Condition':[],'Track':[],'OrgDim':[],'OrgMajDim':[],'OrgMinDim':[]})
	
		if(self.threshLow is None or self.threshHigh is None):
			self.set_color_thresholds()

		if(not os.path.exists(os.path.join(self.track_folder,saveFile)) or overwrite):
			
			nImages = 100
			# generate a list of random indices so the user can find the best possible image
			random_ints = np.array(np.random.randint(0, len(self.imageIndex), nImages), dtype = 'int')
		 
			print(random_ints)
			fileList = np.array(self.data['Image name'][self.imageIndex])
			fileList = fileList[random_ints]
			
			# Calculate based on 100 images

			# Enter an event loop where the use can choose which image is used for calculating organism size
			img_num = 0
			roiFlag = False


			while True:
				try:
					file = fileList.iloc[img_num]
				except:
					file = fileList[img_num]
				print(file)
				image = cv2.imread(os.path.join(self.track_folder,self.image_dict[file],file))

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
							OrgMajDim.append((1/self.pixelPermm)*2*Radius)
							OrgMinDim.append((1/self.pixelPermm)*2*Radius)
							OrgDim.append((1/self.pixelPermm)*2*Radius)

						else:
							OrgMajDim.append((1/self.pixelPermm)*ellipse[1][0])
							OrgMinDim.append((1/self.pixelPermm)*ellipse[1][1])
							OrgDim.append((1/self.pixelPermm)*(ellipse[1][1] + ellipse[1][0])/float(2))

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
						self.set_color_thresholds(overwrite = True)
						
						
				else:
					# Select new color thresholds
					img_num += 1
					key = ord('c')
						



			OrgDim_mean = np.nanmean(np.array(OrgDim))
			OrgMajDim_mean = np.nanmean(np.array(OrgMajDim))
			OrgMinDim_mean = np.nanmean(np.array(OrgMinDim))

			self.obj_diameter = OrgDim_mean
			self.obj_diameter_maj = OrgMajDim_mean
			self.obj_diameter_min = OrgMinDim_mean
			
			# with open(os.path.join(self.track_folder,saveFile), 'wb') as f:  # Python 3: open(..., 'wb')
			#     pickle.dump((OrgDim_mean, OrgMajDim_mean, OrgMinDim_mean), f)
			# Save the Organism dimensions to file

			size_df = size_df.append(pd.DataFrame({'Organism':[self.Organism],'Condition':[self.Condition],'Track':[self.track_folder],'object diameter':[self.obj_diameter],'object diameter max':[self.obj_diameter_maj],'object diameter min':[self.obj_diameter_min]}))
						
			size_df.to_csv(os.path.join(self.track_folder, saveFile))

		else:
			# Load the Organism Size data
			print('Loading organism size from memory ...')
			
			# with open(os.path.join(self.track_folder,saveFile), 'rb') as f:
			#     OrgDim_mean, OrgMajDim_mean, OrgMinDim_mean = pickle.load(f)
			size_df = pd.read_csv(os.path.join(self.track_folder,saveFile))
		
			self.obj_diameter = size_df['object diameter'][0]
			self.obj_diameter_maj = size_df['object diameter max'][0]
			self.obj_diameter_min = size_df['object diameter min'][0]
			
		print('*'*50)
		print('object diameter {} mm'.format(self.obj_diameter))
		print('object diameter max {} mm'.format(self.obj_diameter_maj))
		print('object diameter min {} mm'.format(self.obj_diameter_min))
		print('*'*50)
					

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
		
		
		
	
		
	def compute_background_fluid_velocity(self, image_a, image_b, deltaT = 1, overwrite_piv = False, overwrite_velocity = False, masking = False, obj_position = None, obj_size = 0.1):
		"""
		Computes the mean fluid velocity given a pair of images, 
		far away from objects (if any objects are present).
		
		"""        
		#--------------------------------------------------------------------------
		# Load the frame-pair into memory
		#--------------------------------------------------------------------------
		frame_a_color = cv2.imread(os.path.join(self.track_folder, self.image_dict[image_a], image_a))
		frame_b_color = cv2.imread(os.path.join(self.track_folder, self.image_dict[image_b], image_b))

		# Plot the object's position on the image to verify it is correct
#        print('Circle diameter: {}'.format(int(2*self.scaleFactor*self.obj_diameter*self.pixelPermm)))
#        frame_a_color_copy = np.copy(frame_a_color)
#        
#        cv2.circle(frame_a_color_copy, (int(obj_position[0]), int(obj_position[1])), int(self.obj_diameter*self.pixelPermm*self.scaleFactor), [255,255,255])
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
			x,y,u,v, sig2noise = PIV_Functions.doPIV(frame_a_color,frame_b_color, dT = deltaT, win_size = self.piv_settings['window size'], overlap = self.piv_settings['overlap'], searchArea = self.piv_settings['search area'], apply_clahe = False)
			
			u, v, mask = PIV_Functions.pivPostProcess(u,v,sig2noise, sig2noise_min = 1.5, smoothing_param = 0)
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
			PIV_Functions.overlayPIVdata(frame_a_color, x, y, u, v, Centroids = obj_position)

		   
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
				radius = round((self.obj_diameter/2.0)*self.pixelPermm)
				
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
		PIV_Functions.overlayPIVdata(frame_a_color, x, y, u, v, Centroids = obj_position)
		# Find the mean velocity
		u_avg, v_avg = (np.nanmean(u), np.nanmean(v))
		u_std, v_std = (np.nanstd(u), np.nanstd(v))
		
		print(u_avg)
		print(v_avg)
		return u_avg, v_avg, u_std, v_std
 
	  
	def compute_fluid_velocity_timeseries(self, overwrite_velocity = False, masking = False):
		# 
		""" Computes the averge background fluid velocity at each time point during which an image is available and stores the result.

			For a tracked object, this velocity is equal and opposite to the object's velocity relative to the fluid. 
		
		For each pair of time-points with images:
			> Generate a mask for each image as necessary
			> Do PIV on the pair of images
			> Calculate the average velocity of the fluid in regions far from the object
			> Store the average velocity as a finction of time
		"""
		self.fluid_velocity_folder = os.path.join(self.track_folder, 'FluidVelocityTimeseries')
		self.fluid_velocity_file = 'fluid_velocity_timeseries_'+ str(self.Tmin)+'_' + str(self.Tmax) + '.csv'
		self.fluid_velocity_path = os.path.join(self.fluid_velocity_folder, self.fluid_velocity_file)

		self.initialize_piv_analysis()

		if(not os.path.exists(self.fluid_velocity_folder)):
			os.makedirs(self.fluid_velocity_folder)

		

		if(not os.path.exists(self.fluid_velocity_path) or overwrite_velocity):
			
			print("calculating fluid velocity time series ...")
		
			n_image_pairs = len(self.imageIndex)-1

			
			self.u_avg_array = np.zeros(n_image_pairs)
			self.v_avg_array = np.zeros(n_image_pairs)
			self.u_std_array = np.zeros(n_image_pairs)
			self.v_std_array = np.zeros(n_image_pairs)
			image_pairs = []
			self.imageIndex_array = np.zeros(n_image_pairs, dtype='int')
			Time_array = np.zeros(n_image_pairs)
			
			for ii in tqdm(range(n_image_pairs)):
						 
				imageindex_a = self.imageIndex[ii] 
				imageindex_b = self.imageIndex[ii + 1]
				
				Time_array[ii] = self.data['Time'][imageindex_a]
				self.imageIndex_array[ii] = imageindex_a

				image_a = self.data['Image name'][imageindex_a]
				image_b = self.data['Image name'][imageindex_b]

				image_pairs.append((image_a, image_b))
				
				try:
					obj_position = (self.imW/2 - round(self.data['X_image'][imageindex_a]*self.pixelPermm), self.imH/2 - round(self.data['Z_image'][imageindex_a]*self.pixelPermm))
				except:
					# If X-image data is unavailable assume it is 0
					obj_position = (self.imW/2, self.imH/2 - round(self.data['Z_image'][imageindex_a]*self.pixelPermm))
				
				
				# First check if both these images exist in memory
				try:
					image_a_exists = os.path.exists(os.path.join(self.track_folder, self.image_dict[image_a], image_a))
				except:
					image_a_exists = False
				try:
					image_b_exists = os.path.exists(os.path.join(self.track_folder, self.image_dict[image_b], image_b))
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
					self.u_avg_array[ii], self.v_avg_array[ii], self.u_std_array[ii], self.v_std_array[ii] = self.compute_background_fluid_velocity(image_a,image_b,deltaT = dT, masking = masking, obj_position = obj_position, obj_size = self.obj_diameter, overwrite_piv = False)
				 
				# If either of those images do not exist, assume that the velocity remains constant over the missing frames
				elif(not image_a_exists or not image_b_exists):
					print('One or more of image pair not found...')
					print('Checking for next image index...')
					self.u_avg_array[ii], self.v_avg_array[ii], self.u_std_array[ii], self.v_std_array[ii] = u_avg_array[ii-1], self.v_avg_array[ii-1], self.u_std_array[ii-1], self.v_std_array[ii-1] 
					continue
			   
			
			# with open(self.FluidVelocitySavePath, 'wb') as f:  # Python 3: open(..., 'wb')
			#         pickle.dump((self.imageIndex_array, self.u_avg_array, self.v_avg_array, self.self.u_std_array, self.self.v_std_array), f)
			
			# Save the fluid velocity time series into a CSV file
			df_fluid_velocity = pd.DataFrame({'Time': Time_array, 'Image pairs': image_pairs,  'fluid velocity X mean':self.u_avg_array, 'fluid velocity Z mean':self.v_avg_array, 'fluid velocity X std':self.u_std_array, 'fluid velocity Z std':self.v_std_array})

			df_fluid_velocity.to_csv('FluidVelocityTimeseries_{}_{}'.format(self.Tmin, self.Tmax))
			
		else:
			print("Fluid time series found! Loading ...")
			# Load data from the CSV file
			df_fluid_velocity = pd.read_csv(self.fluid_velocity_path)


			# with open(self.FluidVelocitySavePath, 'rb') as f:  # Python 3: open(..., 'wb')
			#         self.imageIndex_array, self.u_avg_array, self.v_avg_array, self.self.u_std_array, self.self.v_std_array = pickle.load(f)
				
				
  
# 	def compute_fluid_relative_disp(self, overwrite_flag = False, save = False):
# 		 """ Computes the object's displacement relative to the fluid.
# 			Prequisites: 
# 			- Stage velocity relative to lab.
# 			- Object's velocity with respect to the lab.
# 			- Fluid velocity with respect to the lab.
# 			Note that in gravity machine v2.0 and higher, the optical system is fixed relative to lab.
# 			Therefore displacements relative to the optical system are displacements relative to the lab.
# 			For earlier Gravity machine versions the optical system is fixed in Z but moving in X and Y. 

# 			Stage velocity relative to lab (recorded in track data: Z)
# 			Fluid velocity relative to lab (from PIV)
# 			Object's velocity relative to lab (from object motion in the image)

# 			We want to compute the object's velocity relative to the fluid.

# 			Vector operations to calculate V_objFluid which is what we want.
# 			 Note: Vz is the object velocity relative to the stage V_objStage
# 			 V_objLab: is measured from the displacement of the object centroid in the image. 
# 			 V_objStage = V_objLab - V_stageLab
# 			 Therefore, 
# 				 V_stageLab = V_objLab - VobjStage   ---- (1)
			 
# 			 The measured fluid velocity using PIV is:
# 				 V_measured = V_stageLab + V_fluidStage
# 				 Therefore, 
# 				 V_fluidStage = V_measured - V_stageLab ---- (2)
# 				 We can substitute for V_stageLab from (1) to get V_fluidStage
				 
# 			 Now, 
# 				 V_objFluid = V_objStage - V_fluidStage
				 
# 		"""
# 		# self.compute_fluid_velocity_timeseries(overwrite_velocity = overwrite_flag)

# 		# V_z is the object's velocity relative to the stage, V_z_image is object's velocity relative to the image
# 		Vz_stageLab = self.data['V_z_image'][self.imageIndex_array] - self.data['V_z'][self.imageIndex_array]
# 		Vz_fluidStage = self.v_avg_array - Vz_stageLab
# 		Vz_objFluid = self.data['V_z'][self.imageIndex_array] - Vz_fluidStage
# 		Z_objFluid =  utils.compute_displacement_from_velocity(x_data = self.data['Time'][self.imageIndex_array], y_data = Vz_objFluid)
# 		#----------------------------------------------------------------------------------------
# 		# Correcting X-velocity
# 		#----------------------------------------------------------------------------------------
# 		# For GM v > 2.0 data (Fixed Optical System, Moving Stage in X,Y, Theta)
		
# 		# Note that if the X-centroid of the object is not available then the velocity contribution of V_objLab 
# 		# is assumed to be zero.
# 		Vx_stageLab = -self.Vx[self.imageIndex_array] + self.Vx_objLab[self.imageIndex_array]
# 		Vx_fluidStage = self.u_avg_array - Vx_stageLab
# 		Vx_objFluid = self.Vx[self.imageIndex_array] - Vx_fluidStage
# 		X_objFluid =  utils.compute_displacement_from_velocity(x_data = self.data['Time'][self.imageIndex_array], y_data = Vx_objFluid)
# 		#----------------------------------------------------------------------------------------
# 		# Uncomment block below for GM < v2.0 data
# 		# Optical system fixed in Z, moving in X and Y

# #            Vx_objFluid = self.Vx_objStage[self.imageIndex_array] - self.u_avg_array
# #        
# #            self.X_objFluid =  self.utils.compute_displacement_from_velocity(x_data = self.data['Time'][self.imageIndex_array], 
# #                                                        y_data = Vx_objFluid)  
# 		#----------------------------------------------------------------------------------------
	   
# 		# Store the displacements relative to the fluid in the derived data (X and Z positions become the calculated positions relative to the fluid)
# 		self.derived_data['X'] = X_objFluid
# 		self.derived_data['Z'] = Z_objFluid  
# 		self.derived_data['V_x_objFluid'] = Vx_objFluid 
# 		self.derived_data['V_z_objFluid'] = Vz_objFluid

# 		if(save == True):
# 			# Save the corrected data in a separate CSV file
# 			self.save_corrected_track()


	def save_corrected_track(self, overwrite = True):

		# The remaining data is the same the orginal data (but only sampled at time points where an image exists)
		self.derived_data['Y'] = self.data['Y'][self.imageIndex_array]
		self.derived_data['Time'] = self.data['Time'][self.imageIndex_array]
		self.derived_data['Image name'] = self.data['Image name'][self.imageIndex_array]
		self.derived_data['X_image'] = self.data['X_image'][self.imageIndex_array]
		self.derived_data['Z_image'] = self.data['Z_image'][self.imageIndex_array]

		self.analysis_save_path = os.path.join(self.track_folder, 'track_corrected_{}_{}.csv'.format(self.Tmin, self.Tmax))

		df_analysis_dict = {}

		save_variables = ['Time', 'Image name', 'X', 'Y', 'Z']

		if(overwrite or os.path.exists(self.analysis_save_path)==False):

			for key in save_variables:
				df_analysis_dict[VARIABLE_MAPPING[key]] = self.derived_data[key]

			self.df_analysis = pd.DataFrame(df_analysis_dict)                
			self.df_analysis.to_csv(self.analysis_save_path)

		else:

			self.load_corrected_track()


		
	def load_corrected_track(self):

		if(os.path.exists(self.analysis_save_path)):

			df_analysis = pd.read_csv(self.analysis_save_path)

			for key in VARIABLE_MAPPING:
				if(VARIABLE_MAPPING[key] in list(df_analysis.columns)):
					self.derived_data[key] = np.array(df_analysis[VARIABLE_MAPPING[key]])
				else:
					print('Warning {} not found in input data'.format(key))
					self.derived_data[key] = None

		else:

			print('Analysis data does not exist!')

	#--------------------------------------------------------------------------       
	# Signal Processing Functions
	#--------------------------------------------------------------------------  
	# @@@ Need to modift and use the pandas based method @@@     
  
	#-------------------------------------------------------------------------
	# Plotting functions
	#--------------------------------------------------------------------------
	def plot_displacement_timeseries(self, save = False, save_folder = None):

		title = 'Displacement time series'
		fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True, figsize = (8,12))
		ax1.set_title(title)
		sns.lineplot(x = self.data['Time'], y = self.data['X'], color = 'r', linewidth = 1, label = 'X', ax = ax1, ci = None)
		ax1.set_ylabel('X '+units['X'])
		sns.lineplot(x = self.data['Time'], y = self.data['Y'], color = 'g', linewidth = 1, label = 'Y', ax = ax2, ci = None)
		ax2.set_ylabel('Y '+units['Y'])
		sns.lineplot(x = self.data['Time'], y = self.data['Z'], color = 'b', linewidth = 1, label = 'Z', ax = ax3, ci = None)
		ax3.set_ylabel('Z '+units['Z'])
		ax3.set_xlabel('Time' + units['Time'])

		if(save):
			self.save_plot(fig, title, save_folder)

		plt.show()

	def plot_velocity_timeseries(self, save = False, save_folder = None):

		title = 'Velocity time series'
		fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True, figsize = (8,12))
		ax1.set_title(title)
		sns.lineplot(x = self.data['Time'], y = self.data['V_x'], color = 'magenta', linewidth = 1, ax = ax1, ci = None)
		ax1.set_ylabel('X velocity '+units['V_x'])
		sns.lineplot(x = self.data['Time'], y = self.data['V_y'], color = 'darkolivegreen', linewidth = 1, ax = ax2, ci = None)
		ax2.set_ylabel('Y velocity'+units['V_y'])
		sns.lineplot(x = self.data['Time'], y = self.data['V_z'], color = 'darkblue', linewidth = 1, ax = ax3, ci = None)
		ax3.set_ylabel('Z velocity'+units['V_z'])
		ax3.set_xlabel('Time' + units['Time'])

		if(save):
			self.save_plot(fig, title, save_folder)
			

		plt.show()

	def save_plot(self, fig, title = None, save_folder = None):

		if(save_folder is not None):
				file_path = save_folder
		else:
			file_path = self.plots_folder

		if(not os.path.exists(file_path)):
			os.makedirs(file_path)

		file_name = title 

		for image_type in imgFormats:
			fig.savefig(os.path.join(file_path, file_name + image_type), dpi = 300, bbox_inches = 'tight')

			print('Saved ' + file_name + image_type + ' to disk')




