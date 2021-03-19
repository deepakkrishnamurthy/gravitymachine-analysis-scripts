import numpy as np
import cv2

# -*- coding: utf-8 -*-
"""
Created on Mon May  7 19:44:40 2018

@author: Francois and Deepak
"""

import numpy as np
import cv2
from scipy.ndimage.filters import laplace

#color is a vector HSV whose size is 3


def default_lower_HSV(color):
    c=[0,100,100]
    c[0]=np.max([color[0]-10,0])
    c[1]=np.max([color[1]-40,0])
    c[2]=np.max([color[2]-40,0])
    return np.array(c,dtype="uint8")

def default_upper_HSV(color):
    c=[0,255,255]
    c[0]=np.min([color[0]+10,178])
    c[1]=np.min([color[1]+40,255])
    c[2]=np.min([color[2]+40,255])
    return np.array(c,dtype="uint8")

def threshold_image(image_BGR,LOWER,UPPER):
    image_HSV=cv2.cvtColor(image_BGR,cv2.COLOR_BGR2HSV)
    imgMask = cv2.inRange(image_HSV, LOWER, UPPER)  #The tracked object will be in white
    # imgMask = cv2.erode(imgMask, None, iterations=2) # Do a series of erosions and dilations on the thresholded image to reduce smaller blobs
    # imgMask = cv2.dilate(imgMask, None, iterations=2)
    
    return imgMask

def threshold_image_gray(image_gray, LOWER, UPPER):
    imgMask = np.array((image_gray >= LOWER) & (image_gray <= UPPER), dtype='uint8')
    
    # imgMask = cv2.inRange(cv2.UMat(image_gray), LOWER, UPPER)  #The tracked object will be in white
    # imgMask = cv2.erode(imgMask, None, iterations=2) # Do a series of erosions and dilations on the thresholded image to reduce smaller blobs
    # imgMask = cv2.dilate(imgMask, None, iterations=2)
    
    return imgMask



def bgr2gray(image_BGR):
    return cv2.cvtColor(image_BGR,cv2.COLOR_BGR2GRAY)

def crop(image,center,imSize): #center is the vector [x,y]
    imH,imW,*rest=image.shape  #image.shape:[nb of row -->height,nb of column --> Width]
    xmin = max(10,center[0] - int(imSize))
    xmax = min(imW-10,center[0] + int(imSize))
    ymin = max(10,center[1] - int(imSize))
    ymax = min(imH-10,center[1] + int(imSize))
    return np.array([[xmin,ymin],[xmax,ymax]]),np.array(image[ymin:ymax,xmin:xmax])

def get_bbox(cnt):
    return cv2.boundingRect(cnt)


def find_centroid_enhanced(image,last_centroid):
    #find contour takes image with 8 bit int and only one channel
    #find contour looks for white object on a black back ground
    # This looks for all contours in the thresholded image and then finds the centroid that maximizes a tracking metric
    # Tracking metric : current centroid area/(1 + dist_to_prev_centroid**2)
    contours = cv2.findContours(image, cv2.RETR_TREE,cv2.CHAIN_APPROX_SIMPLE)[-2]
    centroid=False
    isCentroidFound=False
    if len(contours)>0:
        all_centroid=[]
        dist=[]
        for cnt in contours:
            M = cv2.moments(cnt)
            if M['m00']!=0:
                cx = int(M['m10']/M['m00'])
                cy = int(M['m01']/M['m00'])
                centroid=np.array([cx,cy])
                isCentroidFound=True
                all_centroid.append(centroid)
                dist.append([cv2.contourArea(cnt)/(1+(centroid-last_centroid)**2)])

    if isCentroidFound:
        ind=dist.index(max(dist))
        centroid=all_centroid[ind]

    return isCentroidFound,centroid

def find_centroid_enhanced_Rect(image,last_centroid):
    #find contour takes image with 8 bit int and only one channel
    #find contour looks for white object on a black back ground
    # This looks for all contours in the thresholded image and then finds the centroid that maximizes a tracking metric
    # Tracking metric : current centroid area/(1 + dist_to_prev_centroid**2)
    contours = cv2.findContours(image, cv2.RETR_TREE,cv2.CHAIN_APPROX_SIMPLE)[-2]
    centroid=False
    isCentroidFound=False
    rect = False
    if len(contours)>0:
        all_centroid=[]
        dist=[]
        for cnt in contours:
            M = cv2.moments(cnt)
            if M['m00']!=0:
                cx = int(M['m10']/M['m00'])
                cy = int(M['m01']/M['m00'])
                centroid=np.array([cx,cy])
                isCentroidFound=True
                all_centroid.append(centroid)
                dist.append([cv2.contourArea(cnt)/(1+(centroid-last_centroid)**2)])

    if isCentroidFound:
        ind=dist.index(max(dist))
        centroid=all_centroid[ind]
        cnt = contours[ind]
        xmin,ymin,width,height = cv2.boundingRect(cnt)
        xmin = max(0,xmin)
        ymin = max(0,ymin)
        width = min(width, imW - int(cx))
        height = min(height, imH - int(cy))
        rect = (xmin, ymin, width, height)


    return isCentroidFound,centroid, rect

def find_centroid_basic(image):
    #find contour takes image with 8 bit int and only one channel
    #find contour looks for white object on a black back ground
    # This finds the centroid with the maximum area in the current frame
    contours = cv2.findContours(image, cv2.RETR_TREE,cv2.CHAIN_APPROX_SIMPLE)[-2]
    centroid=False
    isCentroidFound=False
    if len(contours)>0:
        cnt = max(contours, key=cv2.contourArea)
        M = cv2.moments(cnt)
        if M['m00']!=0:
            cx = int(M['m10']/M['m00'])
            cy = int(M['m01']/M['m00'])
            centroid=np.array([cx,cy])
            isCentroidFound=True
    return isCentroidFound,centroid

def find_centroid_basic_Rect(image):
    #find contour takes image with 8 bit int and only one channel
    #find contour looks for white object on a black back ground
    # This finds the centroid with the maximum area in the current frame and alsio the bounding rectangle. - DK 2018_12_12
    imH,imW = image.shape
    contours = cv2.findContours(image, cv2.RETR_TREE,cv2.CHAIN_APPROX_SIMPLE)[-2]
    centroid=False
    isCentroidFound=False
    rect = False
    if len(contours)>0:
        # Find contour with max area
        cnt = max(contours, key=cv2.contourArea)
        M = cv2.moments(cnt)

        if M['m00']!=0:
            # Centroid coordinates
            cx = int(M['m10']/M['m00'])
            cy = int(M['m01']/M['m00'])
            centroid=np.array([cx,cy])
            isCentroidFound=True

             # Find the bounding rectangle
            xmin,ymin,width,height = cv2.boundingRect(cnt)
            xmin = max(0,xmin)
            ymin = max(0,ymin)
            width = min(width, imW - xmin)
            height = min(height, imH - ymin)
            rect = (xmin, ymin, width, height)

    return isCentroidFound,centroid, rect

def get_image_center_width(image):
    ImShape=image.shape
    ImH,ImW=ImShape[0],ImShape[1]
    return np.array([ImW*0.5,ImH*0.5]),ImW


def YTracking_Objective_Function(image, color):
    #variance methode
    if(image.size is not 0):
        if(color):
            image = bgr2gray(image)
        mean,std=cv2.meanStdDev(image)
        return std[0][0]**2
    else:
        return 0



def colorThreshold(image, threshLow, threshHigh):
    
    thresh_low = np.array(threshLow,dtype='uint8')
    thresh_high = np.array(threshHigh,dtype='uint8')
    
    img_type = None
    try:
        imW, imH, numChannels = np.shape(image)
        if(numChannels==3):
            img_type = 'color'
    except:
        imW, imH = np.shape(image)
        img_type = 'gs'


    if(img_type == 'color'):
        im_th = threshold_image(image, thresh_low, thresh_high)
    else:
        im_th = threshold_image_gray(image, thresh_low[2], thresh_high[2])
    
    kernel3 = np.ones((3,3),np.uint8)
    kernel5 = np.ones((5,5),np.uint8)
    kernel7 = np.ones((7,7),np.uint8)
    kernel15 = np.ones((15,15),np.uint8)
    
    
    im_th = cv2.morphologyEx(im_th, cv2.MORPH_CLOSE, kernel5)
    
    im_th = cv2.morphologyEx(im_th, cv2.MORPH_CLOSE, kernel7)
    
    im_th = cv2.morphologyEx(im_th, cv2.MORPH_OPEN, kernel7)
    
    im_th = cv2.morphologyEx(im_th, cv2.MORPH_CLOSE, kernel15)

    
#    plt.matshow (im_th , cmap=cm.Greys_r )
#    plt.show()
    
    # Select the largest contour as the final mask
    cnts = cv2.findContours(im_th,cv2.RETR_TREE,cv2.CHAIN_APPROX_SIMPLE)[-2]
    
    print('No:of contours:{}'.format(len(cnts)))
    
    if(len(cnts)>0):
        largestContour = max(cnts, key=cv2.contourArea)
    
        M = cv2.moments(largestContour)
        
        # Approximate the contour to 1% of its arcLength
        epsilon = 0.01*cv2.arcLength(largestContour,True)
        approxContour = cv2.approxPolyDP(largestContour,epsilon,True)
                
        
#        orgCentroid = np.array((int(M["m10"] / M["m00"]), int(M["m01"] / M["m00"])))
        
#        ellipse = cv2.fitEllipse(approxContour)
            
            
            # in mm
    
        return approxContour
    
    else:
        return None