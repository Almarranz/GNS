#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 10:14:07 2022

@author: alvaromartinez
"""

# =============================================================================
# A gasgano kind-of-tool for classifing files in flat, dark or science
# =============================================================================

# =============================================================================
# WARNING!!!!
# For some reason the dark that we have to use for these data are no loger
# atached to the science data. In order to get the right darks you have to donwload
# them directly from the archive (http://archive.eso.org/eso/eso_archive_main.html)
# specify the corret date and the instrument (HAWKI), check the 'calibration'
# box and fill the field TPL ID with 'HAWKI_img_cal_WinDarks'
# =============================================================================
import numpy as np
import os
import glob
from astropy.io import fits
import shutil
# %%
field_n =6
date =  field_n #date of the flats
data = '/Users/alvaromartinez/Desktop/Phd/HAWK/GNS_2/H/Field/%s/'%(field_n)
sky = '/Users/alvaromartinez/Desktop/Phd/HAWK/GNS_2/H/Sky/%s/'%(field_n)
flat = '/Users/alvaromartinez/Desktop/Phd/HAWK/GNS_2/H/Flat/%s/'%(date)
dark = '/Users/alvaromartinez/Desktop/Phd/HAWK/GNS_2/Dark/%s/'%(field_n)
common_path = '/Users/alvaromartinez/Desktop/Phd/HAWK/GNS_2/Data/GNS/2021/H/6/ims/'
field = 'GC H F%s'%(field_n)
# %%
with open(dark + 'list.sof','w') as f:
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
    if header['OBJECT'] == 'DARK' and 'WinDarks' in header['ORIGFILE'] :
        count +=1
        print(len(name)*'*')
        print(name,header['OBJECT'],header['ORIGFILE'])
        print(len(name)*'*')
        with open(dark + 'list.txt','a') as f:
            f.write(os.path.basename(name) + '\n')
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
        print(len(name)*'-')
        print(name,header['OBJECT'],header['ORIGFILE'])
        print(len(name)*'-')

with open(flat + 'flat.sof','a') as ff:
    ff.write(common_path+'dark_ext.fits MASTER_DARK')
    ff.close()

            
            
            
            
            