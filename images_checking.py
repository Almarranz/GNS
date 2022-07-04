#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul  3 16:26:43 2022

@author: alvaromartinez
"""

# =============================================================================
# Load all the dejitter images in the same cube for ocular inspection
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
# chips =[1,2,3,4]
chips =[1]
n_offsets = 70
cube = np.zeros((n_offsets*8,4096,4096))
for chip in chips:
    f = len(glob.glob(common_path+ 'chip%s_cube*.gz'%(chip)))
    # cube = np.zeros((f*8,2048,2048))
    
    i=0
    with open(pruebas +'check_chip_all.txt', 'w' ) as texto:
        texto.write('')
        
    for files in sorted(glob.glob(common_path+ 'chip%s_cube*.gz'%(chip)),key=os.path.getmtime):
        print(os.path.basename(files))
        data = fits.getdata(files)
        for s in range((data.shape[0])):
            if chip ==1:
                cube[i,0:2048,0:2048] = data[s,:,:]
                i += 1
                with open(pruebas +'check_chip_all.txt', 'a' ) as texto:
                    texto.write('%s %s \n'%(i,os.path.basename(files)))
                print(i,os.path.basename(files))
            elif chip ==2:
                 cube[i,2048:4096,2048:4096] = data[s,:,:]
                 i += 1
            elif chip ==3:
                 cube[i,2048:4096,0:2048] = data[s,:,:]  
                 i += 1
            elif chip ==4:
                 cube[i,0:2048,2048:4096] = data[s,:,:]   
                 i += 1
    fits.writeto(pruebas + 'check_all_chip.fits',cube, overwrite=True)             
    
    
   # %% 
# dark_c1 = dark[0:2047,0:2047]
# dark_c2 = dark[2048:4095,0:2047]
# dark_c3 = dark[2048:4095,2048:4095]
# dark_c4 = dark[0:2047,2048:4095]
   
cube_test =cube[i,2048:4096,2048:4096]
    
    
    
    
    
    