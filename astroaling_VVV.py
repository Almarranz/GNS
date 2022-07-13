#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 12 12:08:07 2021

@author: amartinez
"""

#!/usr/bin/env python
# coding: utf-8
# In[1]:


import numpy as np
import matplotlib.pyplot as plt
import astroalign as aa
from astropy.io import fits
import sys
import time
# In[12]:


band='J'
field=6
chip =1 

VVV='/Users/alvaromartinez/Desktop/PhD/HAWK/GNS_2/VVV/Fields/%s/'%(band)
tmpdir = VVV+'/Fields/%s/tmp/'%(band)
GNS_stars= '/Users/alvaromartinez/Desktop/Phd/HAWK/GNS_2/Data/GNS/2021/H/%s/data/'%(field)




# In[3]:


VVV_all=fits.open(VVV+'Field%s.fits.gz'%(field))[0].data

xsize_vvv,ysize_vvv=VVV_all.shape[0],VVV_all.shape[1]
VVV_c1=VVV_all[0:round(xsize_vvv/2)-1,0:round(ysize_vvv/2)-1]



# In[5]
dic_gns={}
dic_f={}
# max_pc = 100
# min_pc = 0
for i in range(1,5):#loads GNS stars in two separate list: coordinates and fluxes
    stars_lst = np.loadtxt(GNS_stars+'stars_%s.txt'%(i))
    print('Stars in orginal list %s'%(len(stars_lst)))
    # flux_max = np.percentile(stars_lst[:,2],max_pc)
    # flux_min = np.percentile(stars_lst[:,2],min_pc)
    # ok_stars=np.where((stars_lst[:,2]<flux_max) & (stars_lst[:,2]>flux_min))
    # stars_good = stars_lst[ok_stars]
    # dic_gns['gns_c%s'%(i)]=stars_lst[ok_stars][:,0:2]
    # dic_f['f_c%s'%(i)]=stars_lst[ok_stars][:,2]
    
    
    dic_gns['gns_c%s'%(i)]=stars_lst[:,0:2]
    dic_f['f_c%s'%(i)]=stars_lst[:,2]
    
    


vvv_all=np.loadtxt(VVV+'Field%s_stars.txt'%(field))
mag_max = np.percentile(vvv_all[:,4],10)
mag_min = np.percentile(vvv_all[:,4],0)
ok_vvv = np.where((vvv_all[:,4]>=mag_min)&(vvv_all[:,4]<=mag_max))
vvv_s=np.loadtxt(VVV+'Field%s_stars.txt'%(field),usecols=(0,1))
vvv_s = vvv_s[ok_vvv]
scale=0.34/0.106
# In[6]:
dic_vvv={}#crops VVV image in four regions (chips like) for astroalign
#crop list in chips
chip1=np.where((vvv_s[:,0]<round(ysize_vvv/2)+100) & (vvv_s[:,1]<round(xsize_vvv/2)+100))
dic_vvv['vvv_c1']=vvv_s[chip1]*scale
chip2=np.where((vvv_s[:,0]>round(ysize_vvv/2)-100) & (vvv_s[:,1]<round(xsize_vvv/2)+100))
dic_vvv['vvv_c2']=vvv_s[chip2]*scale
chip3=np.where((vvv_s[:,0]>round(ysize_vvv/2)-100) & (vvv_s[:,1]>round(xsize_vvv/2)-100))
dic_vvv['vvv_c3']=vvv_s[chip3]*scale
chip4=np.where((vvv_s[:,0]<round(ysize_vvv/2)+100) & (vvv_s[:,1]>round(xsize_vvv/2)-100))
dic_vvv['vvv_c4']=vvv_s[chip4]*scale



#%%
color = ['g','y','orange','r']

for c in range(1,5):
    fig, ax = plt.subplots(1,1,figsize=(10,10))
    ax.set_title('Chip %s color %s'%(c,color[c-1]))
    ax.scatter(vvv_s[:,0]*scale, vvv_s[:,1]*scale,alpha=0.05 )
    ax.scatter(dic_vvv['vvv_c%s'%(c)][:,0], dic_vvv['vvv_c%s'%(c)][:,1], s=10,color =color[c-1])

# sys.exit()
# In[7]:
max_offset =2.5
for i in range(chip,5):
    check_x,check_y=max_offset+1,max_offset+1
    while abs(check_x) >max_offset or abs(check_y)>max_offset  :# only when tranformation gets better than 1 chip desplacement the coordinates are stored
        print('starting aa')
        print('Stars in transformed list %s'%(len(dic_gns['gns_c%s'%(i)])))
        tic = time.perf_counter()
        m,(_,_)= aa.find_transform(dic_gns['gns_c%s'%(i)],dic_vvv['vvv_c%s'%(i)],max_control_points=250)
        print('For chip%s'%(i)+'\n'+"Translation: (x, y) = (%.2f, %.2f)"%(m.translation[0],m.translation[1]))
        print("Rotation: %.3f degrees"%(m.rotation * 180.0 / np.pi))
        print("Scale factor: %.4f"%(m.scale))
        toc = time.perf_counter()
        print('it took %.2f seconds'%(toc-tic))
        test_gns = aa.matrix_transform(dic_gns['gns_c%s'%(i)], m.params)
        
        
        print(20*'#'+'\n'+'CHECKING'+'\n'+20*'#')
        tic = time.perf_counter()
        t,(_,_)= aa.find_transform(test_gns,dic_vvv['vvv_c%s'%(i)],max_control_points=250)
        print("Translation: (x, y) = (%.2f, %.2f)"%(t.translation[0],t.translation[1]))
        print("Rotation: %.3f degrees"%(t.rotation * 180.0 / np.pi))
        print("Scale factor: %.4f"%(t.scale))
        toc = time.perf_counter()
        print('checking took %.2f seconds'%(toc-tic))
        check_x= t.translation[0]
        check_y= t.translation[1]
    # sys.exit('line 117')
    print('___TRANSFORMED GNS chip %s (field %s)___'%(i,field))   
    print('Stars in transformed list %s'%(len(dic_gns['gns_c%s'%(i)])))
    dic_gns['gns_c%s'%(i)] = aa.matrix_transform(dic_gns['gns_c%s'%(i)], m.params)
    dic_gns['gns_c%s'%(i)]=np.c_[dic_gns['gns_c%s'%(i)],dic_f['f_c%s'%(i)]]
    # np.savetxt(pruebas + 'aa_stars_%s.txt'%(i),dic_gns['gns_c%s'%(i)])
    np.savetxt(GNS_stars+ 'aa_stars_%s.txt'%(i),dic_gns['gns_c%s'%(i)],fmt ='%.6f %.6f %.6f')

print('DONE')
# sys.exit()












