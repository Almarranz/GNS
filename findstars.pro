; Auxiliary routine taken almost literally from StarFinder's
; take a close look at selected reference stars
; with the option to reject individual sources
; XPSF_EXTRACT_CONFIRM


PRO confirm_stars, image, wnum, display_opt, x_in, y_in, $
			  psfsize, x_out, y_out

	on_error, 2
	if  n_elements(x_in) eq 0 or n_elements(y_in) eq 0  then  return
   	sub_arrays, image, x_in, y_in, psfsize, stack
   	nstars = n_elements(x_in)
   	for  n = 0L, nstars - 1  do begin
   	   xn = -1  &  yn = -1
   	   opt = default_display_opt(stack[*,*,n])
   	   opt.reverse = display_opt.reverse
   	   opt.stretch = display_opt.stretch
   	   opt.color_table = display_opt.color_table
   	   display_image, stack[*,*,n], wnum, OPTIONS = opt
   	   msg = dialog_message('Confirm this star?', /QUESTION)
   	   if  strlowcase(msg) eq 'no'  then begin
   	      x_in[n] = -1  &  y_in[n] = -1
   	   endif
   	endfor
   	w = where(x_in ge 0 and y_in ge 0, n_confirm)
   	if  n_confirm ne 0  then begin
   	   x_out = x_in[w]  &  y_out = y_in[w]
   	endif
	display_image, image, wnum, OPTIONS = display_opt
	return
end

; ==============================

PRO FINDSTARS, field_nr, chip_nr

field_nr = 6
;~ chip_nr = 1
band = 'H'

data_path = '/Users/alvaromartinez/Desktop/PhD/HAWK/GNS_2/data/GNS/2021/'+band+'/' + strn(field_nr) + '/ims/'
out_path = '/Users/alvaromartinez/Desktop/PhD/HAWK/GNS_2/data/GNS/2021/'+band+'/' + strn(field_nr) +  '/data/'
tmp_path = '/Users/alvaromartinez/Desktop/PhD/HAWK/GNS_2/data/GNS/2021/'+band+'/' + strn(field_nr) +  '/tmp/'
;~ data_path = '/data/GNS/2015/'+band+'/' + strn(field_nr) + '/ims/'
;~ out_path = '/data/GNS/2015/'+band+'/' + strn(field_nr) + '/data/'
;~ tmp_path = '/data/GNS/2015/'+band+'/' + strn(field_nr) + '/tmp/'


; Find stars in images of each Chip
; that can serve for image alignment.


; General StarFinder settings
; ---------------------------
psf_size = 60.
maskrad = 20
threshold = [5.]
back_box = 20.
deblend = 0
min_correlation = 0.7
weigh_med = 0
unweighted = 1
upper_lev = 3.0e4
n_fwhm_back = 9.0
n_fwhm_fit = 2.0
n_fwhm_match = 1.0
mag_fac = 2L
n_width = 3.0
norm_max = 1
correl_mag = 2.0
niter = 2
compbg = 0
rel_thresh = 1
guide_x = ""
guide_y = ""

;loop chips
; ----------------

for i_chip = chip_nr, chip_nr do begin
  chip_nr = strn(i_chip)
  im = readfits(data_path + 'lnx_jitter_' + strn(i_chip)+ '.fits')
  ;~ im = readfits(data_path + 'lnx_jitter_' + strn(i_chip)+ '.fits.gz')
  sz = size(im)
  n1 = sz[1]
  n2 = sz[2]
  noise = readfits(data_path + 'lnx_jitter_' + strn(i_chip)+ '_sig.fits')
  ;~ noise = readfits(data_path + 'lnx_jitter_' + strn(i_chip)+ '_sig.fits.gz')

 ; choose PSF reference stars
 ; --------------------------

  device, decomposed = 0
  boxsize = 31. ; width of box within which maximum is searched after a click
  disp_opt = default_display_opt(im)
  disp_opt.stretch = 'logarithm'
  disp_opt.range = max(im)*[1.e-3,1.0]
  disp_opt.large = 1
  display_image, im, wnum, OPTIONS = disp_opt, MODIFY_OPT = modify_opt
  click_on_max, im, /MARK, BOXSIZE = boxsize, x_0, y_0
  measure_centroid, im, x_0, y_0, boxsize
  ques = dialog_message('Selection OK?', /QUESTION)
  if strlowcase(ques) eq 'no' then begin
   confirm_stars, im, wnum, disp_opt, x_0, y_0, psf_size, x_psf, y_psf
  endif else begin
   x_psf = x_0
   y_psf = y_0
  endelse
  psf_extract,x_psf, y_psf, x_secondary, y_secondary, im, $
      psf_size, psf, psf_fwhm, back, $
      N_FWHM_BACK = n_fwhm_back, N_FWHM_FIT = n_fwhm_fit, $
      INTERP_TYPE = 'I', UPPER_LEVEL = upper_lev, $
      N_FWHM_MATCH = n_fwhm_match, N_WIDTH = n_width, $
      MAG_FAC = mag_fac, UNWEIGHTED = (weigh_med eq 0) and 1B, $
      NORM_MAX = (norm_max eq 1) and 1B, /CUBIC

   writefits, 'tmppsf.fits', psf
;  STOP

   MASK_PSF, psf, maskrad, PSF_MASKED=psf_masked, PSF_OFFSET=psf_offset, WINGS=wings
   psf = psf_masked
   psf = psf/total(psf)  ; normalization of PSF
   writefits, data_path + 'psf_chip' + chip_nr + '.fits', psf
;   STOP

   forprint, TEXTOUT= out_path + 'psf_fwhm_chip' + chip_nr + '.txt', psf_fwhm*0.106, /NOCOMMENT


 ;    Extract PSF for field and re-run StarFinder
 ; --------------------------------------------

   threshold = [3.,3.]
   back_box = 0.
   deblend = 0
   starfinder, im, psf, X_BAD=x_bad, Y_BAD = y_bad, $
        BACKGROUND = background, BACK_BOX = back_box, $
        threshold, REL_THRESHOLD = rel_thresh, /PRE_SMOOTH, $
        NOISE_STD = noise, min_correlation, $
        CORREL_MAG = correl_mag, INTERP_TYPE = 'I', $
	ESTIMATE_BG = compbg, DEBLEND = deblend, N_ITER = niter, SILENT=0, $
	GUIDE_X = guide_x, GUIDE_Y = guide_y, $
	SV_SIGMA_R = 0.0085, SV_SIGMA_A= 0.0050, $
      	x, y, f, sx, sy, sf, c, STARS = stars, $
        LOGFILE = logfilename, /CUBIC
   forprint, TEXTOUT= out_path + 'stars' + '_' + chip_nr + '.txt', /NOCOMMENT, x, y, f

; OPTIONAL OUTPUT FOR DEBUGGING
; -----------------------------
;   writefits, tmp_path + 'im' + '_' + chip_nr + '.fits' , im
;   writefits, data_path + 'resid' + '_' + chip_nr + '.fits' , im - stars, /COMPRESS
   ;~ writefits, data_path + 'stars' + '_' + chip_nr + '.fits' , stars
   ;~ writefits, data_path + 'stars' + '_' + chip_nr + '.fits' , stars, /COMPRESS
;   dat = ptr_new({X_size: 20, Y_size: 20, Sigma_x: 1.5, Sigma_y: 1.5, Angle: 0.0})
;   map = image_model(x,y,f,n1,n2,'gaussian', dat)
 ;  writefits, data_path + 'map' + '_' + chip_nr + '.fits', map, /COMPRESS

endfor

 wdelete, wnum

END

