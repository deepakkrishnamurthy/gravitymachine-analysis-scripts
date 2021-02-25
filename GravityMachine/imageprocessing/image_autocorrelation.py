'''
Script to calculate the auto-correlation/cross-correlation between pairs of images

-Deepak
'''
import numpy as np
from scipy import signal
import os
import cv2

path = '/Users/deepak/Dropbox/GravityMachine/FocusTracking/two-hole aperture/pyro 10x 1.1 mm, 1um dz, 201 points_2020-01-23 17-51-4.544342/cropped'


FileList = os.listdir(path)

print(FileList)

FileList.sort()

count = 0
nImages = 1
for file in FileList:
    
    

    
    print(file)
    
    image = cv2.imread(os.path.join(path, file),0)
    
    print(np.shape(image))
    
    
    cv2.imshow('Image', image)
    cv2.waitKey(1)
    
    
    auto_corr = signal.correlate(image, image, mode = 'same')
    
    
    
    if(count>nImages):
        break
    
    count+=1
    
    
cv2.imshow('auto-corr', auto_corr)
cv2.waitKey(1)