
units = {'Time':'(s)', 'X':'(mm)','Y':'(mm)','Z':'(mm)'}
imgFormat = ['.png', '.svg']
# Map between the header names in the CSV file and the internal variable names
# X_objStage  Y_objStage  Z_objStage  Theta_stage X_image Z_image
VARIABLE_MAPPING = {'Time':'Time', 'X':'X_objStage','Y':'Y_objStage','Z':'Z_objStage','Image name':'DF1', 'X_image':'X_image', 'Z_image':'Z_image'}

VARIABLE_MAPPING_LEGACY = {'Time':'Time', 'X_obj':'Xobj','Y_obj':'Yobj','Z_obj':'ZobjWheel','Image name':'Image name', 'X_image':'Xobj_image', 'Z_image':'Zobj'}
