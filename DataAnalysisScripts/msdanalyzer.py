# -*- coding: utf-8 -*-
"""
Created on Fri Sep 21 14:45:49 2018
Mean-Squared-Displacement Analyzer for multiple tracks
@author: deepak90
"""
import numpy as np
# import Track
import os
import imp
import matplotlib.pyplot as plt
import scipy.interpolate as interpolate
from cycler import cycler
import pandas as pd
import GravityMachineTrack
import PlotFunctions.PlotUtils as PlotUtils
imp.reload(GravityMachineTrack)
from lmfit import minimize, Parameters

def squaredDisp(data):
    
    return (data - data[0])**2

class msdanalyzer:
    
    def __init__(self, Tracks = None, nDims = 3, testFlag = False):
        # Tracks is a list of gravity machine tracks
        self.precision = 12
        self.nDims = 3
        self.testFlag = testFlag
        if(self.testFlag == 1 or Tracks is None):
            self.maxDelay = 50
            self.timeWindow = 1
            # Length of the delays window over which we calculate MSD
            self.tLen = 100
            # self.genTestTraj()
            self.genRunTumbleTraj()
            
        else:
            self.Tracks = Tracks
            self.nTracks = len(Tracks)
             # Maximum delay over which we want to calculate correlations
            self.maxDelay = 50
            self.timeWindow = 5
            self.getMinDiff()
            self.getCommonDelays()
       
        
    def genTestTraj(self):
        
        Tracks = []
        
        self.nTracks = 100
        
        boxSize = 2
        
        timeSteps = 1000
        self.DiffCoeff = 100
        dT = 0.05
        
        k = (self.DiffCoeff*dT)**(1/2)
        
        TimeArray = np.array(range(0,timeSteps))*dT
        
        for ii in range(self.nTracks):
            
            
            X0 = boxSize*np.random.rand(self.nDims)
            dX = k*np.random.randn(self.nDims,timeSteps)
            
            dX[:,0] = X0;
            
            X = np.cumsum(dX,1)
            
            Tracks.append(GravityMachineTrack.gravMachineTrack())
            Tracks[ii].T = TimeArray
            Tracks[ii].X = X[0,:]
            Tracks[ii].Y = X[1,:]
            Tracks[ii].ZobjWheel = X[2,:]
            
        
        self.Tracks = Tracks
        
        self.delaysCommon = np.linspace(0,self.maxDelay,self.tLen)
#        print(Tracks[0].T)
            
        
    def genRunTumbleTraj(self, runTime = 1):


        Tracks = []
        
        self.nTracks = 100
        
        timeSteps = 10000

        boxSize = 2
        
        dT = 0.1


        X_array = np.zeros((self.nDims, self.nTracks, timeSteps))
        
        TimeArray = np.array(range(0,timeSteps))*dT

        self.velocity = 1

        # Mean rate of tumbles
        self.tumble_rate = 0.1

        # Initial particle positions
        X0 = boxSize*np.random.rand(self.nDims, self.nTracks)

        # Array to store particle orientation
        P0 = np.zeros((self.nDims, self.nTracks))

        P = np.zeros_like(P0)

        rand_array =np.random.rand(self.nDims - 1, self.nTracks)

        phi = 2*np.pi*rand_array[0,:]

        if(self.nDims>2):
            q = 2.*rand_array[1,:]-1;
        else:
            q = np.zeros((1,self.nTracks))

        cos_theta = q
        sin_theta = (1-q**2)**(1/2)
       
        P0[0,:] = sin_theta*np.cos(phi)
        P0[1,:] = sin_theta*np.sin(phi)
        P0[2,:] = cos_theta

        pdot = np.sum(P0*P0, axis = 0)**(1/2)

        P0[0,:] = P0[0,:]/pdot
        P0[1,:] = P0[1,:]/pdot
        P0[2,:] = P0[2,:]/pdot

        P = P0
        X_array[:,:,0] = X0

        print(X_array[:,:,0])

        for ii in range(timeSteps-1):


            tumble_rand = np.array(np.random.rand(self.nTracks) <= self.tumble_rate*dT, dtype = 'bool')

          

            N_tumble = np.sum(tumble_rand)
            print(N_tumble)

            rand_array =np.random.rand(self.nDims - 1, N_tumble)

            phi = 2*np.pi*rand_array[0,:]

            if(self.nDims>2):
                q = 2.*rand_array[1,:]-1
            else:
                q = np.zeros((1,N_tumble))

            cos_theta = q;
            sin_theta = (1-q**2)**(1/2)
            
            P_new = np.zeros((self.nDims, N_tumble))

            P_new[0,:] = sin_theta*np.cos(phi)
            P_new[1,:] = sin_theta*np.sin(phi)
            P_new[2,:] = cos_theta

            

            # Update the orientations of particles undergoing a tumble
            
            P[:, tumble_rand] = P_new

            pdot = np.sum(P*P, axis = 0)**(1/2)

            # print(pdot)

            P[0,:] = P[0,:]/pdot
            P[1,:] = P[1,:]/pdot
            P[2,:] = P[2,:]/pdot



            X_array[:,:,ii+1] = X_array[:,:,ii] + self.velocity*P*dT


        for ii in range(self.nTracks):

            Tracks.append(GravityMachineTrack.gravMachineTrack())
            Tracks[ii].T = TimeArray
            Tracks[ii].X = X_array[0,ii,:]
            Tracks[ii].Y = X_array[1,ii,:]
            Tracks[ii].ZobjWheel = X_array[2,ii,:]

            if(ii==1):
                print(Tracks[ii].ZobjWheel)
                plt.figure()
                plt.plot(Tracks[ii].T, Tracks[ii].X,'k')
                plt.show()



        self.Tracks = Tracks
        
        self.delaysCommon = np.linspace(0,self.maxDelay,self.tLen)







    def getMinDiff(self):
        
        
        minDiff = 10
        maxDiff = 0
        for ii in range(self.nTracks):
            TimeArray = self.Tracks[ii].T
#            time = np.hstack((time, TimeArray))
            
            diff = np.diff(TimeArray)
#            print('Track {}, Time Array : {}'.format(ii,TimeArray))
            
            minDiff = min(np.min(diff), minDiff)
            maxDiff = max(np.max(diff), maxDiff)
           
        
#        time = np.unique(time)
        
#        print(np.shape(time))
        
#        time = np.round(time, decimals = self.precision)
        
#        self.time = time
        self.minDiff = minDiff
        self.maxDiff = maxDiff
        
#        print(self.minDiff)
#        print(self.maxDiff)
        
    def getCommonDelays(self):
        # This generates a common time array so many different tracks with different sampling rates can be avearaged
        
        self.tLen = int(self.maxDelay/float(self.minDiff))
        self.delaysCommon = np.linspace(0,self.maxDelay,self.tLen)
        self.delayIndex = np.array(range(0,len(self.delaysCommon)))
        
    def DispAtCommonTimes(self):
        
        Disp_X = []
        Disp_Y = []
        Disp_Z = []
        
       
        maxTrackLen = 0
        
        for ii in range(self.nTracks):


            currTrack = self.Tracks[ii]
            
            TimeArray = np.arange(0,currTrack.T[-1],self.minDiff)
            
            print('...Analyzing Track {} of Duration {} s'.format(ii, np.max(currTrack.T)))

            TimeArray = currTrack.T
            
            func_X = interpolate.interp1d(currTrack.T,currTrack.X, kind = 'linear')
            func_Y = interpolate.interp1d(currTrack.T,currTrack.Y, kind = 'linear')
            func_Z = interpolate.interp1d(currTrack.T,currTrack.ZobjWheel, kind = 'linear')
            
            Disp_X.append(func_X(TimeArray))
            Disp_Y.append(func_Y(TimeArray))
            Disp_Z.append(func_Z(TimeArray))
            
            maxTrackLen = max(maxTrackLen, len(TimeArray))
            
        return maxTrackLen, Disp_X, Disp_Y, Disp_Z
    
    
    def computeMeanDisp(self):
        maxTrackLen, Disp_X, Disp_Y, Disp_Z = self.DispAtCommonTimes()
        
        nSubTracks = len(Disp_X)
        
        print('Total number of subtracks : {}'.format(nSubTracks))
        
        Disp_matrix_X = np.zeros((nSubTracks,maxTrackLen))
        Disp_matrix_Y = np.zeros((nSubTracks,maxTrackLen))
        Disp_matrix_Z = np.zeros((nSubTracks,maxTrackLen))
        
#        print(np.shape(MSD_matrix_X))
        
       
        for ii in range(nSubTracks):
            
#            print(np.shape(SquaredDisp_X[ii]))
            Disp_matrix_X[ii,:len(Disp_X[ii])] =  Disp_X[ii]
            Disp_matrix_Y[ii,:len(Disp_Y[ii])] =  Disp_Y[ii]
            Disp_matrix_Z[ii,:len(Disp_Z[ii])] =  Disp_Z[ii]
            
        
        Weights = np.sum(Disp_matrix_X!=0, axis = 0)
            
        Disp_matrix_X[Disp_matrix_X==0] = np.nan
        Disp_matrix_Y[Disp_matrix_Y==0] = np.nan
        Disp_matrix_Z[Disp_matrix_Z==0] = np.nan
        
        
        
        
        
        
        self.meanDisp_X = np.nanmean(Disp_matrix_X,axis = 0)
        self.meanDisp_Y = np.nanmean(Disp_matrix_Y,axis = 0)
        self.meanDisp_Z = np.nanmean(Disp_matrix_Z,axis = 0)
        
        self.meanDisp_XY = np.nanmean(Disp_matrix_X + Disp_matrix_Y, axis = 0)
        
        self.stdDisp_X = np.nanstd(Disp_matrix_X,axis = 0)
        self.stdDisp_Y = np.nanstd(Disp_matrix_Y,axis = 0)
        self.stdDisp_Z = np.nanstd(Disp_matrix_Z,axis = 0)
        
        self.stdDisp_XY = np.nanstd(Disp_matrix_X + Disp_matrix_Y, axis = 0)
        
    
    def computeSqDisp(self):
        
        
        SquaredDisp_X = []
        SquaredDisp_Y = []
        SquaredDisp_Z = []
        
        
        for ii in range(self.nTracks):
            
            
            currTrack = self.Tracks[ii]
        
            print('...Analyzing Track {} of Duration {} s'.format(ii, np.max(currTrack.T)))
            TimeArray = currTrack.T
            
            func_X = interpolate.interp1d(currTrack.T,currTrack.X, kind = 'linear')
            func_Y = interpolate.interp1d(currTrack.T,currTrack.Y,kind = 'linear')
            func_Z = interpolate.interp1d(currTrack.T,currTrack.ZobjWheel,kind = 'linear')
            
            
            
            
            T_start = TimeArray[0]
            startIndex = 0
            counter = 0
            while (T_start < np.max(TimeArray)-self.maxDelay):
                
                T_subTrack = T_start + self.delaysCommon
                
                
                X_subTrack = func_X(T_subTrack)
                Y_subTrack = func_Y(T_subTrack)
                Z_subTrack = func_Z(T_subTrack)
                
                startIndex += next((i for i,x in enumerate(TimeArray[startIndex:] - TimeArray[startIndex]) if x >= self.timeWindow), None)

                
                T_start = TimeArray[startIndex]
                
                SquaredDisp_X.append(squaredDisp(X_subTrack))
                SquaredDisp_Y.append(squaredDisp(Y_subTrack))
                SquaredDisp_Z.append(squaredDisp(Z_subTrack))
                
                counter += 1
                
            print('no:of Subtracks from Track {}: {}'.format(ii, counter))
        
        return SquaredDisp_X, SquaredDisp_Y, SquaredDisp_Z
    
    
    def computeMSD(self):
        
        SquaredDisp_X, SquaredDisp_Y, SquaredDisp_Z = self.computeSqDisp()
        
        
        nSubTracks = len(SquaredDisp_X)
        
        print('Total number of subtracks : {}'.format(nSubTracks))
        
        MSD_matrix_X = np.zeros((nSubTracks,self.tLen))
        MSD_matrix_Y = np.zeros((nSubTracks,self.tLen))
        MSD_matrix_Z = np.zeros((nSubTracks,self.tLen))
        
#        print(np.shape(MSD_matrix_X))
        
       
        for ii in range(nSubTracks):
            
#            print(np.shape(SquaredDisp_X[ii]))
            MSD_matrix_X[ii,:len(SquaredDisp_X[ii])] =  SquaredDisp_X[ii]
            MSD_matrix_Y[ii,:len(SquaredDisp_Y[ii])] =  SquaredDisp_Y[ii]
            MSD_matrix_Z[ii,:len(SquaredDisp_Z[ii])] =  SquaredDisp_Z[ii]
            
#        subTracks = np.array(range(nSubTracks))
        
#        print(len(subTracks))
#        print(len(self.delaysCommon))
#        print(np.shape(MSD_matrix_X))
#        plt.figure(4)
#        ax = plt.contourf(self.delaysCommon,subTracks,MSD_matrix_X==0)
#        plt.colorbar(ax)
#        plt.axis('image')
#        plt.show()
            
        Weights = np.sum(MSD_matrix_X!=0, axis = 0)
        
        
#        plt.figure(4)
#        plt.plot(self.delaysCommon, Weights, color = 'r')
#        plt.show()
#        
#        print(np.shape(Weights))
#        print(np.shape(MSD_matrix_X))
        
        MSD_matrix_X[MSD_matrix_X==0] = np.nan
        MSD_matrix_Y[MSD_matrix_Y==0] = np.nan
        MSD_matrix_Z[MSD_matrix_Z==0] = np.nan
        
        
        
        
        
        
        self.MSD_X = np.nanmean(MSD_matrix_X,axis = 0)
        self.MSD_Y = np.nanmean(MSD_matrix_Y,axis = 0)
        self.MSD_Z = np.nanmean(MSD_matrix_Z,axis = 0)
        
        self.MSD_XY = np.nanmean(MSD_matrix_X + MSD_matrix_Y, axis = 0)
        
        self.stdev_X = np.nanstd(MSD_matrix_X,axis = 0)
        self.stdev_Y = np.nanstd(MSD_matrix_Y,axis = 0)
        self.stdev_Z = np.nanstd(MSD_matrix_Z,axis = 0)
        
        self.stdev_XY = np.nanstd(MSD_matrix_X + MSD_matrix_Y, axis = 0)
        
        
            
            
        
                
#------------------------------------------------------------------------------
        # Plotting functions
#------------------------------------------------------------------------------
    def plotMSD(self, figname = 1, saveFolder = None, orgName = None, plot_analytical = False):
        
        
#        f, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True)
        
        plt.figure(figsize = (8,6))
        ax1 = plt.gca()
        PlotUtils.errorfill(self.delaysCommon, self.MSD_X, self.stdev_X, ax = ax1, color = 'r', label = 'Horizontal (X)')

        ax3 = plt.gca()
        
        PlotUtils.errorfill(self.delaysCommon, self.MSD_Z, self.stdev_Z, ax = ax3, color = 'b', label = 'Vertical (Z)')
#        
        if(plot_analytical==True and self.testFlag==True):
            ax3.plot(self.delaysCommon, self.DiffCoeff*self.delaysCommon, color = 'k', marker = 'o')


        ax3.set_xlabel('Time (s)')
        ax3.set_ylabel('MSD')
        plt.xlim(0,np.max(self.delaysCommon))
        plt.legend(loc='upper left')
        if(saveFolder is not None and orgName is not None):
            plt.savefig(os.path.join(saveFolder,orgName+'_MSD_Linear.svg'),bbox_inches='tight',dpi=150)

        plt.show()
        
        
        # 
        time = np.linspace(10,100,100)
        A = 1
        B = 0.01
        X_linear = B*time
        X_quad = A*(time**2)
        
        
        plt.figure(figsize = (4,3))
        ax1 = plt.gca()
        ax1.plot(self.delaysCommon, self.MSD_X, color = 'r', label = 'Horizontal (X)')
        ax1.set_yscale('log')
        ax1.set_xscale('log')
        
        ax3 = plt.gca()
        
        ax3.plot(self.delaysCommon, self.MSD_Z, color = 'b', label = 'Vertical')#        ax3.set_yscale('log')
        ax3.plot(time,X_linear,'k-')
        ax3.plot(time,X_quad, 'k--')
        ax3.set_yscale('log')
        ax3.set_xscale('log')
        ax3.set_xlabel('Time (s)')
        ax3.set_ylabel('MSD')
        
        plt.xlim(0,np.max(self.delaysCommon))
#        plt.legend(loc='upper left')
        if(saveFolder is not None and orgName is not None):

            plt.savefig(os.path.join(saveFolder,orgName+'_MSD_Log.svg'),bbox_inches='tight',dpi=150)

        plt.show()






        
        
        
        
  
        
        
        

# def main():
    
  

    
    
    
    
    
    
# if __name__ == '__main__':
#     main()    
        

        
        
