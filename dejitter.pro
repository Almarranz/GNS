pro dejitter, field, chip

;field = 10
;chip = 4

for chip = chip, chip do begin

band = 'H'
indir = '/Users/alvaromartinez/Desktop/PhD/HAWK/GNS_2/data/GNS/2021/'+band+'/' + strn(field) + '/ims/'
outdir = '/Users/alvaromartinez/Desktop/PhD/HAWK/GNS_2/data/GNS/2021/'+band+'/' + strn(field) + '/cubes/'
pruebas= '/Users/alvaromartinez/Desktop/PhD/HAWK/GNS_2/pruebas/'

;~ pruebas= '/Users/amartinez/Desktop/PhD/HAWK/GNS_2/pruebas/'

;~ outdir=pruebas
;~ indir = '/data/GNS/2015/'+band+'/' + strn(field) + '/ims/'
;~ outdir = '/data/GNS/2015/'+band+'/' + strn(field) + '/cubes/'

; choose final image size large enough to accommodate jittered
; images for 2048 x 768 pixels with 60" jitter box
; in x: (2048 * 0.106 + 60.)/0.106 = 2614.04
; in y: (768 * 0.106 + 60.)/0.106 = 1334.04
; Later on, the images will be chopped up in so-called sub-cubes.
; The sub-cubes are 600 x 600 pixels and are shifted by 300 pixels
; To make this fit nicely, we now round the above numbers to
; 2700 and 1500

sizex_new = 2700
sizey_new = 2700

; width in x and y 
; of small region used 
; for fine alignment
wx = 1000
wy = 1000

;elements_cube = 8 no all cubes have the same numer of elements now.

; Set the jitter box width 
; equal to or larger than the jitter box used 
; during the observations.
; This value is needed to make sure that the dejittered images
; are made large enough so that no data fall outside the frame.
;
; The jitter_offset paraemter makes sure that an exposure
; with zero offset falls into the middle of the frame after alignment
;~ jitter_offset = round((60./2.)/0.106) ; maximal offset of a jittered exposure from initial pointing
;jitter offset = 330 ; This value has been chosen empirically (systematics in FITS Header cumoffset?)
jitter_offset = round((2700-2048)/2) ; This value has been chosen empirically (systematics in FITS Header cumoffset?)

n_offsets = 70
;~ n_offsets = 68


openw, out1, outdir + 'list_chip'+strn(chip)+'.txt', /get_lun
openw, out2, outdir + 'masklist_chip'+strn(chip)+'.txt', /get_lun
openw, out3, 'jitter_offsets'+strn(chip)+'.txt', /get_lun

new_ref = fltarr(sizex_new,sizey_new)
;~ new_ref_cube = fltarr(sizex_new,sizey_new,elements_cube)
new_mask_ref = fltarr(sizex_new,sizey_new)
lnx_cube = fltarr(sizex_new,sizey_new,n_offsets)
mask_cube = fltarr(sizex_new,sizey_new,n_offsets)


for i = 1, n_offsets do begin


  print, 'Pointing : ' + strn(i)   

  ;~ printf, out1, 'cube_jitter_' + strn(chip) + '_' + strn(i) + '.fits'
  ;~ printf, out2, 'mask_jitter_' + strn(chip) + '_' + strn(i) + '.fits'
  printf, out1, 'cube_jitter_' + strn(chip) + '_' + strn(i) + '.fits.gz'
  printf, out2, 'mask_jitter_' + strn(chip) + '_' + strn(i) + '.fits.gz'
  
  ;~ cube = readfits(indir + 'chip' + strn(chip) + '_cube' + strn(i) + '.fits', header)
  ;~ mask = readfits(indir + 'mask_chip' + strn(chip) + '.fits')
  cube = readfits(indir + 'chip' + strn(chip) + '_cube' + strn(i) + '.fits.gz', header)
  mask = readfits(indir + 'mask_chip' + strn(chip) + '.fits.gz')
  
  ;################################################################
  sz=size(cube)
  elements_cube = sz[3] ; many cubes have fewer number of slices
  
  
  new_ref_cube = fltarr(sizex_new,sizey_new,elements_cube-1) ; smaller than input cube because we omit the long exposure (last frame)
	
  ;################################################################
  ; Read in cumulative offset from initial pointing in pixels
  x_off = strsplit(header[619],'HIERARCH ESO SEQ CUMOFFSETX = ', ESCAPE = '/', /extract)
  y_off = strsplit(header[620],'HIERARCH ESO SEQ CUMOFFSETY = ', ESCAPE = '/', /extract) 
  
  
  x_off_header = fix(x_off[0])
  y_off_header = fix(y_off[0])
  
    
  print, 'Offsets from Fits Header: ' + strn(x_off_header) + ', ' + strn(y_off_header)
    
STOP

  pos_x = jitter_offset - x_off_header
  pos_y = jitter_offset - y_off_header

  if i eq 1 then begin
   
   ref_header = header

   n_size = size(mask)
   size_old_x = n_size[1]
   size_old_y = n_size[2]
 
   new_ref_cube[pos_x:pos_x + size_old_x -1 ,pos_y:pos_y + size_old_y - 1,*] = cube[*,*,0:elements_cube-2]
   new_mask_ref[pos_x:pos_x + size_old_x -1 ,pos_y:pos_y + size_old_y - 1] = mask
  
   ;~ writefits, outdir + 'cube_jitter_' + strn(chip) + '_' + strn(i) + '.fits', new_ref_cube, header
   ;~ writefits, outdir + 'mask_jitter_' + strn(chip) + '_' + strn(i) + '.fits', new_mask_ref
   print, 'cube_jitter_' + strn(chip) + '_' + strn(i) + '.fits'
   help, new_ref_cube
   writefits, outdir + 'cube_jitter_' + strn(chip) + '_' + strn(i) + '.fits', new_ref_cube, header, /COMPRESS
   writefits, outdir + 'mask_jitter_' + strn(chip) + '_' + strn(i) + '.fits', new_mask_ref, /COMPRESS

   ; store long exposures in reference image    
   new_ref[pos_x:pos_x + size_old_x -1 ,pos_y:pos_y + size_old_y - 1,*] = cube[*,*,elements_cube-1];last exposure is 8-1, not 8. The cube has 8 slices, from 0 to 7.
   ;~ new_ref[pos_x:pos_x + size_old_x -1 ,pos_y:pos_y + size_old_y - 1,*] = cube[*,*,elements_cube]
   lnx_cube[*,*,0] = new_ref[*,*]
   mask_cube[*,*,0] = new_mask_ref[*,*]
 
   ; extract a small region from near image centre
   ; that is used for final fine alignment
   small_ref = fltarr(wx,wy)
   a = sizex_new/2-wx/2
   b = sizey_new/2-wy/2
   small_ref = new_ref[a:a+wx, b:b+wy]
   ; set borders of small image to 1
   small_ref[0:2,*]=1
   small_ref[wx-2:wx,*]=1
   small_ref[*,0:2]=1
   small_ref[*,wy-2:wy]=1

   printf, out3, x_off_header, y_off_header, FORMAT='(2F8.2)'
   
;   writefits, outdir + 'aaa_ref.fits', small_ref

  endif else begin

   new = fltarr(sizex_new,sizey_new)
   new_cube=fltarr(sizex_new,sizey_new,elements_cube-1); smaller than input cube because we omit the long exposure (last frame)
   new_mask = fltarr(sizex_new,sizey_new)
   lnx_tmp = fltarr(sizex_new,sizey_new)
    
   new[pos_x:pos_x + size_old_x -1 ,pos_y:pos_y + size_old_y - 1] = cube[*,*,elements_cube-1]
   ;~ new[pos_x:pos_x + size_old_x -1 ,pos_y:pos_y + size_old_y - 1] = cube[*,*,elements_cube]
   ;new_mask[pos_x:pos_x + size_old_x -1 ,pos_y:pos_y + size_old_y - 1] = mask ;I have called this line later down. 
  
   ; extract a small region from near image centre
   ; that is used for final fine alignment
   small_new_ref = fltarr(wx,wy)
   a = sizex_new/2-wx/2
   b = sizey_new/2-wy/2
   small_new_ref = new[a:a+wx, b:b+wy]
   ; set borders of small image to 1  
   small_new_ref[0:2,*]=1
   small_new_ref[wx-2:wx,*]=1
   small_new_ref[*,0:2]=1
   small_new_ref[*,wy-2:wy]=1
;   writefits, outdir + 'aaaa_ref.fits', small_new_ref
 
;   OFFSET  =  alignoffset(small_ref, small_new_ref, 1) 
;   x_off = round(OFFSET[0])
;   y_off = round(OFFSET[1])


    correl_optimize, small_ref, small_new_ref, x_off, y_off, MAGNIFICATION=4, /NUMPIX ; /NUMPIX is ESSENTIAL
    print, 'Offsets from correlation: ' + strn(x_off) + ', ' + strn(y_off)

    ; offsets in header and in correlation point into differeent directions
    ; not tested....
    printf, out3, (x_off_header - x_off), (y_off_header - y_off), FORMAT='(2F8.2)'
   
   ; This is to catch grave problems with correl_optimize
   ; These problems will hopefully not happen any more with /NUMPIX.
   if abs(x_off) gt 100 or abs(y_off) gt 100 then begin
    x_off = 0
    y_off = 0
    print,'###################################'
    print, 'x_off and y_off CORRECTED to zero'
    print,'###################################'
   endif
 
   new_cube[pos_x + x_off:pos_x + size_old_x -1 + x_off ,pos_y + y_off:pos_y + size_old_y + y_off- 1,*] = cube[*,*,0:elements_cube-2]
   lnx_tmp[pos_x + x_off:pos_x + size_old_x -1 + x_off ,pos_y + y_off:pos_y + size_old_y + y_off- 1] = cube[*,*,elements_cube-1]
   ;~ lnx_tmp[pos_x + x_off:pos_x + size_old_x -1 + x_off ,pos_y + y_off:pos_y + size_old_y + y_off- 1] = cube[*,*,elements_cube]
   new_mask[pos_x + x_off:pos_x + size_old_x -1 + x_off ,pos_y + y_off:pos_y + size_old_y + y_off- 1] = mask; here we are applying the shifting also to the mask
   lnx_cube[*,*,i-1] = lnx_tmp[*,*]
   mask_cube[*,*,i-1] = new_mask[*,*]   

   print, 'cube_jitter_' + strn(chip) + '_' + strn(i) + '.fits'
   help, new_cube
   writefits, outdir + 'cube_jitter_' + strn(chip) + '_' + strn(i) + '.fits', new_cube, /COMPRESS
   writefits, outdir + 'mask_jitter_' + strn(chip) + '_' + strn(i) + '.fits', new_mask, /COMPRESS
  

endelse
; STOP

endfor

free_lun, out1, out2, out3

;~ writefits, outdir + 'lnx_jitter_cube_' + strn(chip) + '.fits', lnx_cube, ref_header
;~ writefits, outdir + 'lnx_jitter_mask_' + strn(chip) + '.fits', mask_cube
writefits, outdir + 'lnx_jitter_cube_' + strn(chip) + '.fits', lnx_cube, ref_header, /COMPRESS
writefits, outdir + 'lnx_jitter_mask_' + strn(chip) + '.fits', mask_cube, /COMPRESS


endfor
;stop
 

end