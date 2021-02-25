import openpiv.tools
import openpiv.process
import openpiv.scaling
import openpiv.validation
import openpiv.filters
import numpy as np

#frame_a  = openpiv.tools.imread( 'exp1_001_a.bmp' )
#frame_b  = openpiv.tools.imread( 'exp1_001_b.bmp' )
rootFolder = '/Users/deepak/Dropbox/GravityMachine/DataAnalysis/TestData/seacucumber9_PIV/'

frame_a = openpiv.tools.imread(rootFolder + 'IMG_31126.tif')

frame_b = openpiv.tools.imread(rootFolder + 'IMG_31127.tif')

u, v, sig2noise = openpiv.process.extended_search_area_piv( frame_a.astype(np.int32), frame_b.astype(np.int32), window_size=24, overlap=12, dt=0.02, search_area_size=64, sig2noise_method='peak2peak' )

x, y = openpiv.process.get_coordinates( image_size=frame_a.shape, window_size=24, overlap=12 )

u, v, mask = openpiv.validation.sig2noise_val( u, v, sig2noise, threshold = 1.1 )

u, v, mask = openpiv.validation.global_val( u, v, (-1000, 2000), (-1000, 1000) )

u, v = openpiv.filters.replace_outliers( u, v, method='localmean', max_iter=10, kernel_size=2)

x, y, u, v = openpiv.scaling.uniform(x, y, u, v, scaling_factor = 96.52 )

openpiv.tools.save(x, y, u, v, mask, 'exp1_001.txt' )

openpiv.tools.display_vector_field('exp1_001.txt', scale=100, width=0.0025)
