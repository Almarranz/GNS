#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 10:38:06 2022

@author: alvaromartinez
"""

import numpy as np
from scipy import ndimage
from astropy.io import fits
import glob

field = '6'
band = 'H'

common_path = '/Users/alvaromartinez/Desktop/PhD/HAWK/GNS_2/data/GNS/2021/' + band + '/' + field +'/ims/'
pruebas = '/Users/alvaromartinez/Desktop/PhD/HAWK/GNS_2/pruebas/'
bpm_name = 'bpm.fits'


debug = 0
bpm = fits.getdata(common_path + bpm_name)

#second struct doesnt group data by the corners, first one does

#struct=np.ones((3,3))
struct=[[0, 1, 0],
 [1, 1, 1],
 [0, 1, 0]] 

# %

mask=np.ones(shape=(bpm.shape))
# makes data from FITS usable in scipy.ndimage
bpm1=bpm.byteswap().newbyteorder()
#labels the strcutures in bpm, and print the number of structures found
id_regions, num_ids =ndimage.label(bpm1, structure=struct)
print(num_ids)

 #sums the values of the pixels in each structres individually
id_sizes =np.array(ndimage.sum(bpm1,id_regions, range(0,num_ids+1)))
print(id_sizes)
#choose regions smaller than 'reg' pixels and makes them zeros
reg =6
area_mask = (id_sizes < reg)
bpm1[area_mask[id_regions]] = 0
#makes the mask white with blcks holes
where_0 = np.where(bpm1 == 0)
where_1 = np.where(bpm1 == 1)
bpm1[where_0] = 1
bpm1[where_1] = 0


mask[:,0:11] = 0
mask[0:10,:] = 0
mask[:,4073:4096] = 0
mask[:,2044:2052] = 0
mask[4087:4096,:] = 0
mask[2043:2053,:] = 0

fits.writeto(common_path+ 'mask_%sreg_py.fits'%(reg-1),bpm1*mask,overwrite=True)












