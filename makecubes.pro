PRO MAKECUBES, field_path, ims, field

field = '9'
band = 'H'

field_path = '/Users/amartinez/Desktop/PhD/HAWK/GNS_2/data/GNS/2021/' + band + '/' + field +'/'
;~ field_path = '/data/GNS/2015/' + band + '/' + field + '/'
ims = 'ims/'

indir = field_path + 'cleaned/'
outdir = field_path + 'ims/'

;n_frames = 8 ; number of frames in cube, including the last one (long exposure)

n_offset = 68



indir_mask = '/Users/amartinez/Desktop/PhD/HAWK/GNS_2/data/GNS/2021/' + band + '/' + field +'/' + '/ims/'
;~ indir_mask = '/data/GNS/2015/' + band + '/' + field + '/ims/'
maskname = 'mask.fits'



; Change format of image cubes
; -----------------------------

; images

for j = 1, n_offset do begin
  
 cube = readfits(indir + 'cube'+ strn(j) + '.fits', header)
 sz=size(cube)
 n_frames=sz[3]; I have eliminated some frames, so now not all cubes have the same number of frames.
 ;~ cube = readfits(indir + 'cube'+ strn(j) + '.fits.gz', header)

 ;Chip 1 
 quad = fltarr(2048,2048,n_frames)  
 quad[*,*,*] = cube[0:2047,0:2047,*]
 writefits, outdir + 'chip1_cube' + strn(j) + '.fits', quad, header
 ;~ writefits, outdir + 'chip1_cube' + strn(j) + '.fits', quad, header,/COMPRESS
 
  
 ;Chip 2 
 quad = fltarr(2048,2048,n_frames)  
 quad[*,*,*] = cube[2048:4095,0:2047,*]
 writefits, outdir + 'chip2_cube' + strn(j) + '.fits', quad, header
 ;~ writefits, outdir + 'chip2_cube' + strn(j) + '.fits', quad, header, /COMPRESS
 
 
  
 ;Chip 3 
 quad = fltarr(2048,2048,n_frames)  
 quad[*,*,*] = cube[2048:4095,2048:4095,*]
 writefits, outdir + 'chip3_cube' + strn(j) + '.fits', quad, header
 ;~ writefits, outdir + 'chip3_cube' + strn(j) + '.fits', quad, header, /COMPRESS
 
 
  
 ;Chip 4 
 quad = fltarr(2048,2048,n_frames)  
 quad[*,*,*] = cube[0:2047,2048:4095,*]
 writefits, outdir + 'chip4_cube' + strn(j) + '.fits', quad, header
 ;~ writefits, outdir + 'chip4_cube' + strn(j) + '.fits', quad, header, /COMPRESS   


 print, 'Created chip cubes for offset' + strn(j)

endfor



print, 'Chip cubes created.'

; mask

mask = readfits(indir_mask + maskname)

;Chip1
mask_chip = mask[0:2047,0:2047]
writefits, outdir + 'mask_chip1.fits', mask_chip
;~ writefits, outdir + 'mask_chip1.fits', mask_chip, /COMPRESS


;Chip2
mask_chip = mask[2048:4095,0:2047]
writefits, outdir + 'mask_chip2.fits', mask_chip
;~ writefits, outdir + 'mask_chip2.fits', mask_chip, /COMPRESS


;Chip3
mask_chip = mask[2048:4095,2048:4095]
writefits, outdir + 'mask_chip3.fits', mask_chip
;~ writefits, outdir + 'mask_chip3.fits', mask_chip, /COMPRESS


;Chip4
mask_chip = mask[0:2047,2048:4095]
writefits, outdir + 'mask_chip4.fits', mask_chip
;~ writefits, outdir + 'mask_chip4.fits', mask_chip, /COMPRESS


END
