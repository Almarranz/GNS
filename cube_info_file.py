#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul  3 17:50:17 2022

@author: alvaromartinez
"""

# =============================================================================
# Creates the cube_info.txt file that will be used by selecgood.pro 
# =============================================================================

import numpy as np
from scipy import ndimage
from astropy.io import fits
import glob
import sys
from astropy.io import fits
import os

field = '6'
band = 'H'

common_path = '/Users/alvaromartinez/Desktop/PhD/HAWK/GNS_2/data/GNS/2021/' + band + '/' + field +'/ims/'
pruebas = '/Users/alvaromartinez/Desktop/PhD/HAWK/GNS_2/pruebas/'
cubes = '/Users/alvaromartinez/Desktop/PhD/HAWK/GNS_2/data/GNS/2021/' + band + '/' + field +'/cubes/'

with open(common_path + 'cube_info.txt','w') as f:
    f.write('')
# Select the cube with bad frames and put them in the 'cube' array
# Number the bad frames with cero in slices. 

with open(common_path + 'bad_cubes_%s.txt'%(field),'w') as f_bad:
    f_bad.write('')

cube = [11,24,31,34,36,41,45,58]    
normal =[1,2,3,4,5,6,7,8]
slices=[[0,2,3,4,5,6,7,8],[0,2,3,4,5,6,7,8],[0,0,0,0,0,6,7,8],[0,0,0,4,5,6,7,8],
        [0,2,3,4,5,6,7,8],[0,0,0,4,5,6,7,8],[0,0,3,4,5,6,7,8],
        [0,2,3,4,5,6,7,8]]
print(len(cube),len(slices))
if len(cube) != len(slices):
    sys.exit('Check the slices!')
for i in range(1,71):
    if i in cube:
        print(i)
        with open(common_path + 'cube_info.txt','a') as f:
            sl = str(slices[cube.index(i)]).replace('[', '').replace(']', '')
            # f.write('%s \n'%(slices[cube.index(i)]))  
            f.write('%s \n'%(sl))  
        
        with open(common_path + 'bad_cubes_%s.txt'%(field),'a') as f_bad:
            f_bad.write(str(i) + '\n')
            
    else:
        with open(common_path + 'cube_info.txt','a') as f:
            nor =str(normal).replace('[', '').replace(']', '')
            f.write(str(nor) + '\n') 












