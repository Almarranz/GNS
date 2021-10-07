PRO DEEPASTRO

;~ tmpdir = './Fields/tmp/'
;~ dir = './Fields/'
VVV='/Users/amartinez/Desktop/PhD/HAWK/GNS_2/VVV/'
tmpdir = VVV+'/Fields/H/tmp/'
;~ tmpdir = 'tmp/'
;~ dir = './'
ZP = 23.171
innam = 'Field9'

im = readfits(VVV+'/Fields/H/'+ innam + '.fits',header)
;~ im = readfits(dir + innam + '.fits.gz',header)
EXTAST, header, astr
ZP = SXPAR(header,'PHOTZP')
print, 'Zero Point: ' + strn(ZP)
sz = size(im)
n1 = sz[1]
n2 = sz[2]

;data = im[122:1452,*]
;mmm, data, skymod, sigma , skew, HIGHBAD = 2.0e4
;print, 'Sky background and sigma: ' + strn(skymod) + ', ' + strn(sigma)
noise = im
;noise[*,*] = sigma
noise[*,*] = 3.0 ; very roughly estimated by eye

psf = readfits(VVV+'/Fields/H/' + innam + '_psf.fits')
psf = psf/total(psf)   ; normalize PSF
; Settings for StarFinder you are most likely to 
; want to play with
; these parameters apply to the final StarFinder run
; not to PSF extraction
; ------------------------------------------------

sf_thresh = [3.,3.]
sf_back_box = 20
deblend = 1
deblost = 0
min_correlation = 0.7


; General settings for StarFinder
; --------------------------

correl_mag = 4
niter = 2
compbg = 1
rel_thresh = 1
guide_x = ""
guide_y = ""


; StarFinder run
;####################

  starfinder, im, psf, X_BAD=x_bad, Y_BAD = y_bad, $
        BACKGROUND = background, BACK_BOX = sf_back_box, $
        sf_thresh, REL_THRESHOLD = rel_thresh, /PRE_SMOOTH, $
        NOISE_STD = noise, min_correlation, $
        CORREL_MAG = correl_mag, INTERP_TYPE = 'I', $
	ESTIMATE_BG = compbg, DEBLEND = deblend, DEBLOST = deblost, $
        N_ITER = niter, SILENT=0, $
	GUIDE_X = guide_x, GUIDE_Y = guide_y, $
	SV_SIGMA_R = 0.0085, SV_SIGMA_A= 0.0050, $
      	x, y, f, sx, sy, sf, c, STARS = stars, $
        LOGFILE = logfilename, /CUBIC

  writefits, VVV+'/Fields/H/' + innam + '_stars.fits', stars
  writefits, VVV+'/Fields/H/' + innam + '_bg.fits', background
  writefits, VVV+'/Fields/H/' + innam + '_resid.fits', im-stars-background

  ;~ writefits, dir + innam + '_stars.fits', stars, /COMPRESS
  ;~ writefits, dir + innam + '_bg.fits', background, /COMPRESS
  ;~ writefits, dir + innam + '_resid.fits', im-stars-background, /COMPRESS

  sigma_gauss = fwhm(psf)/2.355 ; compute standard deviation of Gaussian corresponding to PSF
  dat = ptr_new({X_size: 10, Y_size: 10, Sigma_x: sigma_gauss, Sigma_y: sigma_gauss, Angle: 0.0})
  im = image_model(x,y,f,n1,n2,'gaussian', dat)
  writefits, VVV+'/Fields/H/' + innam + '_map.fits', im


  ; Calibrate and save list
  ; select stars in region with more than covfrac coverage
  nstars = n_elements(f)
  openw, outp, VVV+'/Fields/H/' + innam + '_stars.txt', /get_lun
  m = ZP - 2.5*alog10(f)
  sm = 2.5/alog(10.) * sf/f
  XY2AD, x, y, astr, a, d
  for s = 0, nstars-1 do begin
   xi = round(x[s]) & yi = round(y[s])
   printf, outp, format='(9f13.6)', x[s], y[s], a[s], d[s], m[s], sx[s], sy[s], sm[s], c[s]
  endfor
  free_lun, outp



  print, 'Finished.'

END
