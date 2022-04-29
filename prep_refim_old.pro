PRO PREP_REFIM_old

field = '9'
band = 'H'

; Determine FoV
; We assume 60"  jitter window for all observations, even
; though the 2015 observations were done with a 30" jitter
;
; The size of the fina FoV depends on the windowing in the
; y-direciton. It is 768 pixels for DIR = 1.26s and 512 pixels
; for DIT = 0.85s
; EDIT THIS PARAMETER IF NECESSARY
x_size = 2048.
y_size = 2048.
; y_size = 512.
jitterbox = 60.

; Necessary size of extracted VVV image:
; Detector FoV + jitterbox/2. + 10. arcseconds safety margin
; Note that here it is assumed that jitterbox > 15" - gap between the chips !

pxscl_hawki = 0.106 ; pixel scale of HAWK-I
pxscl_vvv = 0.34  ; pixel scale of Omegacam 
x_vvv = round((2*x_size * pxscl_hawki + jitterbox/2. + 10.)/pxscl_vvv)
y_vvv = round((2*y_size * pxscl_hawki + jitterbox/2. + 10.)/pxscl_vvv)

; Use header of first exposure to obtain telescope pointing at start
; of the observations.

list = 'list.txt'
raw_path = '/Users/amartinez/Desktop/PhD/HAWK/GNS_2/'+ band+'/Field/'+field+'/'
pruebas= '/Users/amartinez/Desktop/PhD/HAWK/GNS_2/pruebas/'
;~ raw_path= pruebas
;~ raw_path = '/home/data/raw/2015/' + band +'/Field/' + field + '/'
readcol, raw_path + list, names, FORMAT='A'
im = readfits(raw_path + names[0], header)

; get offset from original pointing )IN RADIANS!!!!!!!!
x_off = strsplit(header[613],' ', /extract)
y_off = strsplit(header[614],' ', /extract) 
print,x_off,y_off
stop
;~ x_off = strsplit(header[425],' ', /extract)
;~ y_off = strsplit(header[426],' ', /extract) 
x_off = float(x_off[5])
y_off = float(y_off[5])

;Convert to degrees
x_off = x_off/3600.
y_off = y_off/3600.
crval1 = SXPAR(header, 'RA')
crval2 = SXPAR(header, 'DEC')
alpha = crval1 - x_off
delta = crval2 - y_off


vvvim = readfits('GC_VVV_J.fits',vvvhdr) ; load b333 from VVV
;~ vvvim = readfits('GC_VVV_J.fits.gz',vvvhdr) ; load b333 from VVV
EXTAST, vvvhdr, astr
AD2XY, alpha, delta, astr, x_cen, y_cen


x0 = round(x_cen - x_vvv/2.)
x1 = x0 + x_vvv
; Some fields are at eastern edge of b333
; this needs to be taken care of
x0_expand = 0
if (x0 lt 0) then begin
  x0_expand = abs(x0)
  x0 = 0
  tmpim = fltarr(x_vvv,y_vvv)
endif
y0 = round(y_cen - y_vvv/2.)
y1 = y0 + y_vvv


hextract, vvvim, vvvhdr, refim, refhdr, x0, x1, y0, y1
; Some fields are at eastern edge of b333
; this needs to be taken care of.
; Careful: The following Will screw up header astrometry.
if (x0_expand gt 0) then begin
  tmpim[x0_expand:x_vvv-1,*] = refim[*,*]
endif

SXADDPAR, refhdr, 'BITPIX', 32
refim = float(refim)
writefits, pruebas + strn(field) + '_J.fits', refim, refhdr
;~ writefits, 'Fields/Field' + strn(field) + '_J.fits.gz', refim, refhdr, /COMPRESS

;stop
END
