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
import shutil
# %%
field_n =6
date = '22-05-28'#date of the flats
data = '/Users/alvaromartinez/Desktop/Phd/HAWK/GNS_2/H/Field/%s/'%(field_n)
sky = '/Users/alvaromartinez/Desktop/Phd/HAWK/GNS_2/H/Sky/%s/'%(field_n)
flat = '/Users/alvaromartinez/Desktop/Phd/HAWK/GNS_2/H/Flat/%s/'%(date)
dark = '/Users/alvaromartinez/Desktop/Phd/HAWK/GNS_2/Dark/%s/'%(field_n)
field = 'GC H F%s'%(field_n)
# %%
with open(dark + 'dark.sof','w') as f:
    f.write('')
with open(flat + 'flat.sof','w') as ff:
    ff.write('')
with open(sky + 'list.txt','w') as fsk:
    fsk.write('')
with open(data + 'list.txt', 'w') as fsc:
    fsc.write('')
count = 0


for file in glob.glob(data + 'HAWKI*.fits'):

    # name = os.path.basename(file)
    name = file
    header = fits.open(file)[0].header
    print(header['OBJECT'])
    if header['OBJECT'] == 'DARK' and 'WinDarks' not in header['ORIGFILE'] :
        count +=1
        # print(name)
        with open(dark + 'dark.sof','a') as f:
            f.write(dark+os.path.basename(name) +' DARK \n')
            f.close
        os.replace(name, dark+os.path.basename(name))
    elif header['OBJECT'] == 'FLAT':
        count +=1
        # print(name)
        with open(flat + 'flat.sof','a') as ff:
            ff.write(flat+os.path.basename(name) +' FLAT_TWILIGHT \n')
            ff.close()
        os.replace(name, flat+os.path.basename(name))
    elif header['OBJECT'] == 'GC_H_Sky':
        count +=1
        with open(sky + 'list.txt', 'a') as fsk:
            fsk.write(os.path.basename(name) + '\n')
            fsk.close()
        os.replace(name, sky+os.path.basename(name))
    elif header['OBJECT'] == field:
        count +=1
        with open(data + 'list.txt', 'a') as fsc:
            fsc.write(os.path.basename(name) + '\n')
            fsc.close()
    else:
        print(len(name)*'*')
        print(name,header['OBJECT'],header['ORIGFILE'])
        print(len(name)*'*')

with open(flat + 'flat.sof','a') as ff:
    ff.write(dark+'darkcomb.fits MASTER_DARK')
    ff.close()

            
            
            
            
            