#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 10:14:07 2022

@author: alvaromartinez
"""

# =============================================================================
# A gasgano kind of tool to classifided files in flat dark or science
# =============================================================================
# Change to pull
import numpy as np
import os
import glob
from astropy.io import fits

data = '/Users/alvaromartinez/Desktop/Phd/HAWKI/GNS_2/raw/field6/'
# %%
field = 'GC H F6'
with open(data + 'dark.sof','w') as f:
    f.write('')
with open(data + 'flat.sof','w') as ff:
    ff.write('')
with open(data + 'sky.sof','w') as fsk:
    fsk.write('')
with open(data + 'science.sof', 'w') as fsc:
    fsc.write(data + 'science.sof')
count = 0


for file in glob.glob(data + 'HAWKI*.fits'):

    # name = os.path.basename(file)
    name = file
    header = fits.open(file)[0].header
    print(header['OBJECT'])
    if header['OBJECT'] == 'DARK':
        count +=1
        # print(name)
        with open(data + 'dark.sof','a') as f:
            f.write(name +' DARK \n')
            f.close
    elif header['OBJECT'] == 'FLAT':
        count +=1
        # print(name)
        with open(data + 'flat.sof','a') as ff:
            ff.write(name +' FLAT_TWILIGHT \n')
            ff.close()
    elif header['OBJECT'] == 'GC_H_Sky':
        count +=1
        with open(data + 'sky.sof', 'a') as fsk:
            fsk.write(name + ' SKY \n')
            fsk.close()
    elif header['OBJECT'] == field:
        count +=1
        with open(data + 'science.sof', 'a') as fsc:
            fsc.write(name + ' SCIENCE \n')
            fsc.close()
    else:
        print(len(name)*'*')
        print(name,header['OBJECT'])
        print(len(name)*'*')
print(count, len(file))

            
            
            
            
            
            