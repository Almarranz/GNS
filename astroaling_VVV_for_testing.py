#!/usr/bin/env python
# coding: utf-8
# In[1]:


import numpy as np
import matplotlib.pyplot as plt
import astroalign as aa
from astropy.io.fits import getheader
from astropy.io import fits
from scipy.spatial import distance
import sys

# In[12]:


band='H'
field=9
chip=1
VVV='/Users/amartinez/Desktop/PhD/HAWK/GNS_2/VVV/Fields/%s/'%(band)
file ='/Users/amartinez/Desktop/PhD/HAWK/GNS_2/data/GNS/2021/%s/%s/ims/'%(band,field)
GNS_stars= '/Users/amartinez/Desktop/PhD/HAWK/GNS_2/data/GNS/2021/%s/%s/data/'%(band,field)
scripts='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/scripts/'
pruebas='/Users/amartinez/Desktop/PhD/HAWK/GNS_2/pruebas/'
stream = open(scripts+'polywarp.py')
read_file = stream.read()
exec(read_file)


# In[3]:


VVV_all=fits.open(VVV+'Field%s.fits'%(field))[0].data
l_exp=fits.open(file+'lnx_jitter_%s.fits'%(chip))[0].data
x,y=VVV_all.shape[0],VVV_all.shape[1]
VVV_c1=VVV_all[0:round(x/2)-1,0:round(y/2)-1]


# In[4]:


# p, (pos_img, pos_img_rot) = aa.find_transform(l_exp,VVV_c1,max_control_points=200
#                                              )
# print("\nTranslation: (x, y) = ({:.2f}, {:.2f})".format(*p.translation))
# print("Rotation: %.3f degrees"%(p.rotation * 180.0 / np.pi))
# print("Rotation: {:.3f} degrees".format(p.rotation * 180.0 / np.pi))
# print("\nScale factor: {:.4f}".format(p.scale))
            


# In[5]
dic_gns={}
dic_f={}
for i in range(1,3):
    dic_gns['gns_c%s'%(i)]=np.loadtxt(GNS_stars+'stars_%s.txt'%(i),usecols=(0,1))
    dic_f['f_c%s'%(i)]=np.loadtxt(GNS_stars+'stars_%s.txt'%(i),usecols=(2))
# gns_s=np.loadtxt(GNS_stars+'stars_%s.txt'%(chip),usecols=(0,1))
# gns_s2=np.loadtxt(GNS_stars+'stars_%s.txt'%(2),usecols=(0,1))

vvv_s=np.loadtxt(VVV+'Field%s_stars.txt'%(field),usecols=(0,1))
scale=0.34/0.106

# In[6]:
dic_vvv={}
#crop list in chips
chip1=np.where((vvv_s[:,0]<round(1453/2)+100) & (vvv_s[:,1]<round(1453/2)+100))
dic_vvv['vvv_c1']=vvv_s[chip1]*scale
chip2=np.where((vvv_s[:,0]>round(1453/2)-100) & (vvv_s[:,1]<round(1453/2)+100))
dic_vvv['vvv_c2']=vvv_s[chip2]*scale
chip3=np.where((vvv_s[:,0]>round(1453/2)-100) & (vvv_s[:,1]>round(1453/2)-100))
dic_vvv['vvv_c3']=vvv_s[chip3]*scale
chip4=np.where((vvv_s[:,0]<round(1453/2)+100) & (vvv_s[:,1]>round(1453/2)-100))
dic_vvv['vvv_c4']=vvv_s[chip4]*scale



#%%
# m,(_,_)= aa.find_transform(gns_s2,vvv_c2,max_control_points=200)
# print('For chip%s'%(2)+'\n'+"Translation: (x, y) = (%.2f, %.2f)"%(m.translation[0],m.translation[1]))
# print("Rotation: %.3f degrees"%(m.rotation * 180.0 / np.pi))
# print("Rotation: %.3f degrees"%(m.rotation * 180.0 / np.pi))
# print("\nScale factor: %.4f"%(m.scale))

# In[7]:

for i in range(1,3):
    check_x,check_y=2,2
    while abs(check_x) >1 or abs(check_y)>1  :
        m,(_,_)= aa.find_transform(dic_gns['gns_c%s'%(i)],dic_vvv['vvv_c%s'%(i)],max_control_points=450)
        print('For chip%s'%(i)+'\n'+"Translation: (x, y) = (%.2f, %.2f)"%(m.translation[0],m.translation[1]))
        print("Rotation: %.3f degrees"%(m.rotation * 180.0 / np.pi))
        print("Scale factor: %.4f"%(m.scale))
        
        test_gns = aa.matrix_transform(dic_gns['gns_c%s'%(i)], m.params)
        
        print(20*'#'+'\n'+'CHECKING'+'\n'+20*'#')
        t,(_,_)= aa.find_transform(test_gns,dic_vvv['vvv_c%s'%(i)],max_control_points=450)
        print("Translation: (x, y) = (%.2f, %.2f)"%(t.translation[0],t.translation[1]))
        print("Rotation: %.3f degrees"%(t.rotation * 180.0 / np.pi))
        print("Scale factor: %.4f"%(t.scale))
        check_x= t.translation[0]
        check_y= t.translation[1]
        
    print('___NOW TRANSFORMING GNS___')   
    dic_gns['gns_c%s'%(i)] = aa.matrix_transform(dic_gns['gns_c%s'%(i)], m.params)
    dic_gns['gns_c%s'%(i)]=np.c_[dic_gns['gns_c%s'%(i)],dic_f['f_c%s'%(i)]]
    np.savetxt(pruebas + 'aa_stars_%s.txt'%(i),dic_gns['gns_c%s'%(i)])
print('Saliendo')
sys.exit()
    # In[8]:

# for i in range(1,3):
#     dic_gns['gns_c%s'%(i)] = aa.matrix_transform(dic_gns['gns_c%s'%(i)], m.params)

#     t,(_,_)= aa.find_transform(dic_gns['gns_c%s'%(i)],dic_vvv['vvv_c%s'%(i)],max_control_points=450)
#     print("Translation: (x, y) = (%.2f, %.2f)"%(t.translation[0],t.translation[1]))
#     print("Rotation: %.3f degrees"%(t.rotation * 180.0 / np.pi))
#     print("Rotation: %.3f degrees"%(t.rotation * 180.0 / np.pi))
#     print("\nScale factor: %.4f"%(t.scale))


# In[9]:






distancia=1
q=1#choose the chip(quadrant) do you want to align
for loop in range(5):
    diff=[]
    for i in range(len(dic_vvv['vvv_c%s'%(q)])):
    #compara las distancia entre los puntos y guarda las menores que a, si hay mas de dos puntos con distancias menores que a, guarda la más perqueña
        dist=distance.cdist(dic_vvv['vvv_c%s'%(q)][i:i+1,0:2],dic_gns['gns_c%s'%(q)][:,0:2], 'euclidean')
        d=np.where(dist<distancia)
        if len(d[1])>0:
            diff.append((dic_vvv['vvv_c%s'%(q)][i],dic_gns['gns_c%s'%(q)][d[1][np.argmin(dist[d])]]))
    print('Common stars %s'%(len(diff)))           
    x1=[]
    y1=[]
    x2=[]
    y2=[]
    for i in range(len(diff)):
        x1.append(diff[i][0][0])
        y1.append(diff[i][0][1])
        x2.append(diff[i][1][0])
        y2.append(diff[i][1][1])

    x1=np.array(x1)
    y1=np.array(y1)
    x2=np.array(x2)
    y2=np.array(y2)

    Kx=[]
    Ky=[]
    degree=1
    Kx,Ky=polywarp(x1,y1,x2,y2,degree=degree)
    #print(Kx[0,0])
    xi=np.zeros(len(dic_gns['gns_c%s'%(q)]))
    yi=np.zeros(len(dic_gns['gns_c%s'%(q)]))
    x=[]
    y=[]
    x=dic_gns['gns_c%s'%(q)][:,0]
    y=dic_gns['gns_c%s'%(q)][:,1]
    x=np.array(x)
    y=np.array(y)
    #for k in range(degree+1):
     #   for m in range(degree+1):
      #      xi=xi+Kx[k,m]*(x**k)*(y**m)
       #     yi=yi+Ky[k,m]*(x**k)*(y**m)
    for k in range(degree+1):
        for m in range(degree+1):
            xi=xi+Kx[k,m]*x**k*y**m
            yi=yi+Ky[k,m]*x**k*y**m
    dic_gns['gns_c%s'%(q)][:,0]=xi
    dic_gns['gns_c%s'%(q)][:,1]=yi
    '''
    diff=[]
    im_a=1
    for i in range(dic_stars['listE_im'+str(im_a)].shape[0]): #compara las distancia entre los puntos y guarda las menores que a, si hay mas de dos puntos con distancias menores que a, guarda la más perqueña
                dist=distance.cdist(dic_stars['listE_im'+str(im_a)][i:i+1,0:2],im2[:,0:2], 'euclidean')
                d=np.where(dist<distancia)
                if len(d[1])>0:
                    #diff.append((dic_stars['stars_im'+str(im_a)][i],im2[d[1][np.argmin(dist[d])]]))
                    diff.append((im2[d[1][np.argmin(dist[d])]]))
    dic_listas['stars_1'+str(l)]=diff
    '''
#     np.savetxt(tmp+'cube_im'+str(l)+'_chip'+str(chip)+'_aligned_py.txt',im2,fmt='%.5f')
#     print('Saved aligned list of im%s'%(l))


# In[ ]:

dic_gns['gns_c%s'%(i)]=np.c_[dic_gns['gns_c%s'%(i)],dic_f['f_c%s'%(i)]]



# In[14]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




