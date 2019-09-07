import numpy as np
import cv2

def colorThreshold(self,image):
    
    thresh_low = np.array(self.threshLow,dtype='uint8')
    thresh_high = np.array(self.threshHigh,dtype='uint8')
    

    
    im_th = cv2.inRange(image, thresh_low, thresh_high)
    
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