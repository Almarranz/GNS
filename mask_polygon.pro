PRO mask_polygon, image, ds9reg, value = value


; PURPOSE:
;		Mask polygon regions in an image using ds9-format polygon files.
;
; INPUT:
;		image - 2D image array to be masked
;
;		ds9reg - string, name of the ds9 region file containing the polygon 
;			information. This file can be generated with ds9 of course.  
;
;           Or, any text file that contains lines starting with 'polygon('
;           (no space or any other characters in front) and end with ')' 
;			can be used.  Lines not starting with 'polygon(' will be igored.  
;			The format should be:
;				polygon(x1,y1,x2,y2,x3,y3,x4,y4...,xn,yn)
;			for polygons with n points.  The polygon loop should not be 
;			closed, i.e., xn not equal to x1, same for y.  There should
;			be no space at all in the line.  xn and yn can be integers or
;			floats.  They should be IMAGE COORDINATES.  RA/Dec will not be
;			accepted.  The points CAN be outside the image.  Multiple polygons
;			can be masked at the same time in one polygon file.
;		value - value to fill in the masked regions.  Default is zero, if
;			not provided.  All polygons will be masked with the same value.
;
;
;


field = '9'
band = 'H'
ds9reg='/Users/amartinez/Desktop/PhD/HAWK/GNS_2/H/regions/mask_c4.reg'
;~ mask_path = '/Users/amartinez/Desktop/PhD/HAWK/GNS_2/data/GNS/2021/' + band + '/' + field +'/ims/'
;~ image=mask_path+ 'mask.fits

pruebas= '/Users/amartinez/Desktop/PhD/HAWK/GNS_2/pruebas/'
image=readfits(pruebas +'mask.fits')

IF n_elements(value) EQ 0 THEN value = 0.0
imsize = size(image,/dimen)




openr, ds9, ds9reg, /get_lun
string=''



WHILE NOT EOF(ds9) DO BEGIN   ; main loop
readf, ds9, string
print, string
IF strmid(string, 0, 8) NE 'polygon(' THEN goto, skip   

; extract coordinate components
coors = strmid(string,8,strlen(string)-9)
coors = float(strsplit(coors,',',/extract))
npoint = n_elements(coors)/2
ind = indgen(npoint)
x = coors[ind*2]-1
y = coors[ind*2+1]-1
merge_vector, x, x[0]  ; to close the loop
merge_vector, y, y[0]

minx = max([min(x),0])
maxx = min([max(x),imsize[0]-1])
miny = max([min(y),0])
maxy = min([max(y),imsize[1]-1])
print,minx,maxx,miny,maxy

subimage = image[minx:maxx, miny:maxy]
print,size(subimage)

subsize = size(subimage,/dimen)
print,subsize

xind = lindgen(subsize[0],subsize[1]) mod subsize[0]
yind = lindgen(subsize[0],subsize[1]) / subsize[0]
x = float(x-minx)
y = float(y-miny)

;xind=177.-minx
;yind=303.-miny


; calculate sum of angles
;theta = fltarr(subsize[0],subsize[1])
theta=0.0
FOR i=0,npoint-1 DO BEGIN
   theta1 = atan(y[i]-yind, x[i]-xind)
   theta2 = atan(y[i+1]-yind,x[i+1]-xind)

   dtheta = theta2 - theta1
   A = where(dtheta GT !pi)
   IF total(A) NE -1 THEN dtheta[A] = dtheta[A] - 2*!pi 
   B = where(dtheta LT -!pi)
   IF total(B) NE -1 THEN dtheta[B] = dtheta[B] + 2*!pi
   
   theta = theta + dtheta
;   print,theta1,theta2,dtheta,theta
ENDFOR



A = where(abs(theta) GT !pi) 
IF total(A) NE -1 THEN subimage[A] = value
image[minx:maxx, miny:maxy] = subimage

skip:
ENDWHILE  ; end of main loop
free_lun, ds9

writefits, pruebas+'new_mask.fits',image

END
