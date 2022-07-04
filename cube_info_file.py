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

with open(pruebas + 'cube_info.txt','w') as f:
    f.write('')
    
cube = [58]    
normal =[1,2,3,4,5,6,7,8]
slices=[[0,2,3,4,5,6,7,8]]
for i in range(1,71):
    if i in cube:
        print(i)
        with open(pruebas + 'cube_info.txt','a') as f:
            sl = str(slices[cube.index(i)]).replace('[', '').replace(']', '')
            # f.write('%s \n'%(slices[cube.index(i)]))  
            f.write('%s \n'%(sl))  
    else:
        with open(pruebas + 'cube_info.txt','a') as f:
            nor =str(normal).replace('[', '').replace(']', '')
            f.write(str(nor) + '\n') 












