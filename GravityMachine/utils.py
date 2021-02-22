# General utility functions
import numpy as np

def velocity_central_diff(time, data):
        
    assert len(time) == len(data), 'Time and data lengths not equal'
    res = all(i < j for i, j in zip(time, time[1:]))
    assert res == True, "Time should be strictly increasing" 
    
    data_len = len(time)
    velocity = np.zeros(data_len)
    
    for i in range(data_len):
        if(i == 0):
            # Forward difference at the start points
            velocity[i] = (data[i+1] - data[i])/(time[i+1] - time[i])
        
        elif(i == data_len-1):
            # Backward difference at the end points
            velocity[i] = (data[i] - data[i-1])/(time[i] - time[i-1])
    
        else:
            # Central difference for all other points
            velocity[i] = (data[i+1] - data[i-1])/(time[i+1] - time[i-1])
            
    return velocity