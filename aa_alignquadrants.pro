PRO AA_ALIGNQUADRANTS,chip

; PURPOSE: Loads the list for GNS stars already transfom with rotation and translation (python astroalign).
; Then search for more stars and refine the
; alignment.
; 
; This script loops over the four quadrants for each pointing.
; The pointing ('field') must be set by hand.

; ------------EDIT HERE---------
; size of VVV reference image
; This is the printed output  the end of prep_refim.pro
xsize_ref = 1453
ysize_ref = 1453
;~ ysize_ref = 655

;~ chip = 1      
chip_nr = strn(chip)
band = 'H'
field_nr = 9

; ----------------CAN BE EDITED, but USUALLY NOT NECESSARY ----------
 

; To create small aligned longexposure images
; So that we can cut away the unnecessary zeros around the images.
xs = 2400
ys = 2400
;~ ys = 1200
; Approximate offsets of fields within aligned HAWK-I frame
; These offsets can be seen inside the 
; images lnx_jitter_?_aligned.fits.gz
; that are being produced by this script
x_off = [0,2200,2200,0]
y_off = [0,0,2200,2200]
;~ y_off = [50,50,900,900]

; size of dejittered HAWK-I long exposures
xsize_quad = 2700
ysize_quad = 2700
;~ ysize_quad = 1500
scale=0.34/0.106

; --------------------------------------------------------------------
VVV='/Users/amartinez/Desktop/PhD/HAWK/GNS_2/VVV/'
im_path = '/Users/amartinez/Desktop/PhD/HAWK/GNS_2/data/GNS/2021/'+band+'/' + strn(field_nr) + '/ims/'
data_path = '/Users/amartinez/Desktop/PhD/HAWK/GNS_2/data/GNS/2021/'+band+'/' + strn(field_nr) + '/data/'
tmp_path = '/Users/amartinez/Desktop/PhD/HAWK/GNS_2/data/GNS/2021/'+band+'/' + strn(field_nr) + '/tmp/'
ref_file =  VVV +'/Fields/H/Field' + strn(field_nr) + '_stars.txt'
pruebas='/Users/amartinez/Desktop/PhD/HAWK/GNS_2/pruebas/'

;~ im_path = '/home/data/GNS/2015/'+band+'/' + strn(field_nr) + '/ims/'
;~ data_path = '/home/data/GNS/2015/'+band+'/' + strn(field_nr) + '/data/'
;~ tmp_path = '/home/data/GNS/2018/'+band+'/' + strn(field_nr) + '/tmp/'
;~ ref_file = '/home/data/VVV/ForHAWKI/J/Fields/Field' + strn(field_nr) + '_stars.txt'



; Read list of reference stars
readcol, ref_file, x_ref, y_ref, a_ref, d_ref, m_ref, sx_ref, sy_ref, sm_ref, c, Format = 'A,A,A,A,A,A,A,A,A'


x_ref=float(x_ref)
y_ref=float(y_ref)
a_ref=float(a_ref)
d_ref=float(d_ref)
m_ref=float(m_ref)
sx_ref=float(sx_ref)
sy_ref=float(sy_ref)
sm_ref=float(sm_ref)
c=float(c)


; scale ref pixel scale
; -----------------------------
x_ref_scaled = x_ref * scale
y_ref_scaled = y_ref * scale

xsize_final = round(xsize_ref * scale)
ysize_final = round(ysize_ref * scale)

;stop
; now read 'original' HAWK-I stars. We need the original list for Kx, Ky 
 readcol, data_path + 'stars_' + chip_nr + '.txt', x, y, f, Format='A,A,A'
 x=float(x)
 y=float(y)
 f=float(f)
; now read transformed HAWK-I stars list
 readcol, data_path + 'aa_stars_' + chip_nr + '.txt', xi, yi, fi, Format='A,A,A'
 ;~ readcol, data_path + 'stars_' + chip_nr + '.txt', x, y, f, Format='A,A,A'
 xi=float(xi)
 yi=float(yi)
 fi=float(fi)
 
 dmax = 1.0
 compare_lists, x_ref_scaled, y_ref_scaled, xi, yi, x1c, y1c, x2c, y2c, MAX_DISTANCE=dmax, SUBSCRIPTS_1=subc1, SUBSCRIPTS_2 = subc2, SUB1 = sub1, SUB2 = sub2
 nc = n_elements(subc1)
 print, 'Found ' + strn(nc) + ' common stars.'

 ; iterative degree 1 alignment
 ; ------------------------------

  lim_it = 1
  count=0
  comm=[]
  it=0
  lim_it=1 ;cosider convergece when the number of common stars  'lim_it' times.
	 
  while count lt lim_it do begin
  it=it+1
 ;~ for it = 1, 3 do begin
  degree = 1
  polywarp, x_ref_scaled[subc1], y_ref_scaled[subc1], x[subc2], y[subc2], degree, Kx, Ky
  print, Kx
  print, Ky
  xi = Kx[0,0] + Kx[0,1]*x + Kx[1,0]*y + Kx[1,1]*x*y
  yi = Ky[0,0] + Ky[0,1]*x + Ky[1,0]*y + Ky[1,1]*x*y
  compare_lists, x_ref_scaled, y_ref_scaled, xi, yi, x1c, y1c, x2c, y2c, MAX_DISTANCE=dmax, SUBSCRIPTS_1=subc1, SUBSCRIPTS_2 = subc2, SUB1 = sub1, SUB2 = sub2
  nc = n_elements(subc1)
  print, 'Iteration ' + strn(it)
  print, 'Found ' + strn(nc) + ' common stars.'
  comm=[comm,nc]
    if (n_elements(comm) gt 2) then begin
	   if comm[-2] ge comm[-1] then begin
	   count=count+1
	   endif else begin
	   count=0
	   endelse
	endif
  endwhile
;~ endfor
;STOP
 ; iterative degree 2 alignment
 ; ------------------------------

 print, 'Now Degree 2 alignment.'
 lim_it = 1
  count=0
  comm=[]
  it=0
  
	 
  while count lt lim_it do begin
  it=it+1
 ;~ for it = 1, 3 do begin
  degree = 2
  polywarp, x_ref_scaled[subc1], y_ref_scaled[subc1], x[subc2], y[subc2], degree, Kx, Ky
  print, Kx
  print, Ky
  xi = Kx[0,0] + Kx[0,1]*x + Kx[1,0]*y + Kx[1,1]*x*y + Kx[0,2]*x^2 + Kx[1,2]*x^2*y + Kx[2,2]*x^2*y^2 + Kx[2,0]*y^2 + Kx[2,1]*y^2*x 
  yi = Ky[0,0] + Ky[0,1]*x + Ky[1,0]*y + Ky[1,1]*x*y + Ky[0,2]*x^2 + Ky[1,2]*x^2*y + Ky[2,2]*x^2*y^2 + Ky[2,0]*y^2 + Ky[2,1]*y^2*x
  compare_lists, x_ref_scaled, y_ref_scaled, xi, yi, x1c, y1c, x2c, y2c, MAX_DISTANCE=dmax, SUBSCRIPTS_1=subc1, SUBSCRIPTS_2 = subc2, SUB1 = sub1, SUB2 = sub2
  nc = n_elements(subc1)
  print, 'Iteration ' + strn(it)
  print, 'Found ' + strn(nc) + ' common stars.'
  comm=[comm,nc]
    if (n_elements(comm) gt 2) then begin
	   if comm[-2] ge comm[-1] then begin
	   count=count+1
	   endif else begin
	   count=0
	   endelse
	endif
  endwhile
  ;~ endfor

 ; determine transformation parameters for image and save them
 ; for later use
 ; ------------------------------------------------------------
 polywarp,  x[subc2], y[subc2], x_ref_scaled[subc1], y_ref_scaled[subc1], degree, Kx, Ky
 
  x_fin = Kx[0,0] + Kx[0,1]*x + Kx[1,0]*y + Kx[1,1]*x*y + Kx[0,2]*x^2 + Kx[1,2]*x^2*y + Kx[2,2]*x^2*y^2 + Kx[2,0]*y^2 + Kx[2,1]*y^2*x 
  y_fin = Ky[0,0] + Ky[0,1]*x + Ky[1,0]*y + Ky[1,1]*x*y + Ky[0,2]*x^2 + Ky[1,2]*x^2*y + Ky[2,2]*x^2*y^2 + Ky[2,0]*y^2 + Ky[2,1]*y^2*x 
 
 SAVE, Kx, Ky, FILENAME=data_path + 'AlignPars_chip' + chip_nr

 ; transform image and mask
 ; mask pixels that have coverage < 0.9 to avoid 
 ; later problems with division by small numbers
 ; ---------------------------------------------

 im = readfits(im_path + 'lnx_jitter_'+chip_nr + '_wt.fits')
 ;~ im = readfits(im_path + 'lnx_jitter_'+chip_nr + '_wt.fits.gz')
 transim = POLY_2D(im,Kx,Ky,2,xsize_final,ysize_final,CUBIC=-0.5,MISSING=0)
; transim = transim * transmask
 writefits, im_path + 'lnx_jitter_'+chip_nr+ '_aligned_wt.fits', transim
 ;~ writefits, im_path + 'lnx_jitter_'+chip_nr+ '_aligned_wt.fits.gz', transim, /COMPRESS

 im = readfits(im_path + 'lnx_jitter_'+chip_nr + '.fits')
 ;~ im = readfits(im_path + 'lnx_jitter_'+chip_nr + '.fits.gz')
 transim = POLY_2D(im,Kx,Ky,2,xsize_final,ysize_final,CUBIC=-0.5,MISSING=0)
; transim = transim * transmask
 writefits, im_path + 'lnx_jitter_'+chip_nr+ '_aligned.fits', transim
 ;~ writefits, im_path + 'lnx_jitter_'+chip_nr+ '_aligned.fits.gz', transim, /COMPRESS

 transim_im = transim

 im = readfits(im_path + 'lnx_jitter_'+chip_nr + '_sig.fits')
 ;~ im = readfits(im_path + 'lnx_jitter_'+chip_nr + '_sig.fits.gz')
 transim = POLY_2D(im,Kx,Ky,2,xsize_final,ysize_final,CUBIC=-0.5,MISSING=0 )
; transim = transim * transmask
 writefits, im_path + 'lnx_jitter_'+chip_nr+ '_aligned_sig.fits', transim
 ;~ writefits, im_path + 'lnx_jitter_'+chip_nr+ '_aligned_sig.fits.gz', transim, /COMPRESS
 
 transim_noise = transim

 ; to check alignment with VVV, create
 ; a transformed image inside the VVV frame of reference
 polywarp,   x[subc2], y[subc2], x_ref[subc1], y_ref[subc1], degree, Kx, Ky
 im = readfits(im_path + 'lnx_jitter_'+chip_nr + '.fits')
 ;~ im = readfits(im_path + 'lnx_jitter_'+chip_nr + '.fits.gz')
 transim = POLY_2D(im,Kx,Ky,2,xsize_ref,ysize_ref,CUBIC=-0.5,MISSING=0)
 writefits, im_path + 'lnx_jitter_'+chip_nr+ '_VVV.fits', transim
 ;~ writefits, im_path + 'lnx_jitter_'+chip_nr+ '_VVV.fits.gz', transim, /COMPRESS

 im = readfits(im_path + 'lnx_jitter_'+chip_nr + '_wt.fits')
 ;~ im = readfits(im_path + 'lnx_jitter_'+chip_nr + '_wt.fits.gz')
 transim = POLY_2D(im,Kx,Ky,2,xsize_ref,ysize_ref,CUBIC=-0.5,MISSING=0)
 writefits, im_path + 'lnx_jitter_'+chip_nr+ '_VVV_wt.fits', transim
 ;~ writefits, im_path + 'lnx_jitter_'+chip_nr+ '_VVV_wt.fits.gz', transim, /COMPRESS
 

; Now cut out the aligned HAWK-I images from the large frames to get
; rid of the zeros
; --------------------------------------------------------------------

  xlo = x_off[chip-1]
  ylo = y_off[chip-1]
  xhi = x_off[chip-1] + xs - 1
  yhi = y_off[chip-1] + ys - 1

lnx = transim_im[xlo:xhi,ylo:yhi]
lnx_noise = transim_noise[xlo:xhi,ylo:yhi]
 
 
 
writefits, data_path + 'lnx_aligned_' + chip_nr + '.fits', lnx
writefits, data_path + 'lnx_aligned_' + chip_nr + '_sig.fits', lnx_noise
;~ writefits, data_path + 'lnx_aligned_' + chip_nr + '.fits.gz', lnx, /COMPRESS
;~ writefits, data_path + 'lnx_aligned_' + chip_nr + '_sig.fits.gz', lnx_noise, /COMPRESS


;To have an idea of the displacement that we have, we take the initial list and the corrected one to see the difference in positions.


; median global offset
x_dif = median(x_fin - x)
y_dif = median(y_fin - y) 
 
; origin of arrows in displaced image
x0_new = x_dif + x
y0_new = y_dif + y
 
; max size of image
x_dis = max(x_fin-x0_new)
y_dis = max(y_fin-y0_new)
 

print,  strn(x_dis*0.106) + 'arcsec    ' + strn(y_dis*0.106) + 'arcsec'

dx = x_dis*0.106
dy = y_dis*0.106




;stop

 
forprint, TEXTOUT= data_path + 'displ_arcsec_distortion_chip' + chip_nr + '.txt', dx, dy, /NOCOMMENT 

;Drawing the detector


set_plot,'PS', /interpolate

device, XOFFSET=0, YOFFSET=0, $
   FILENAME= data_path + 'Alignment_resid_chip' + chip_nr +'.eps', XSIZE=24., YSIZE=11., $
   /portrait, /color, BITS_PER_PIXEL=8, encapsulated=1
!P.MULTI=[0,1,0]
!P.CHARSIZE=1.8
!X.THICK=4
!Y.THICK=4
!P.THICK=4.0
!P.CHARTHICK=4

!P.COLOR=0

;cgoplot, x1c, y1c, Color='orange', PSYM=3

xf = (-xi[subc2]+x1c)*60
yf = (-yi[subc2]+y1c)*60
; subtract global shift
xf = xf - median(xf)
yf = yf - median(yf)
x0 = xi[subc2]
y0 = yi[subc2]
xf = xf+xi[subc2]
yf = yf+yi[subc2]

ran = RANDOMU(Seed, n_elements(xf))

num = sort(ran)
num = num[0:330]
a = xi[subc2[num]]
b = yi[subc2[num]]

cgplot, a, b, Color='black', XRANGE = [0,xsize_quad], YRANGE = [0,ysize_quad], XTITLE='X-Axis [pixels]', YTITLE='Y-Axis [pixels]', XSTYLE=1, PSYM=3, YSTYLE=1
;cgplot, a, b, Color='black', XRANGE = [2200,4600], YRANGE = [900,2000], XTITLE='[X]', YTITLE='[Y]', XSTYLE=1, PSYM=3, YSTYLE=1


ARROW, x0[num], y0[num], xf[num], yf[num], /DATA, COLOR=4, HSIZE=100


device, /close




END
