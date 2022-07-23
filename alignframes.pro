PRO ALIGNFRAMES, field_nr, chip_nr


; PURPOSE: Uses the image distortion solutions found by 
; ALIGNQUADRANTS to align all the individual frames.
; NOTE: This script also removes the first exposures from each cube
;       becasue they are generally contaminated by vertical stripes 
;       particularly on chip 3

;----- EDIT HERE-----
;field_nr = 10
;chip_nr = 1

 ; ------THESE PAREMETERS MAY NEED TO BE ADJUSTED------
; see alignquadrants.pro
; Image size is the size of the cut-out LNX images for each chip
; They are located in the data directory!

x_out = 2700
y_out = 2700
;~ x_out = 2460
;~ y_out = 2460

; offsets of fields within transformed frame (see alignquadrants.pro)
; used to cut out the chips from the large algined field
;~ x_off = [0,0,0,0]
;~ y_off = [0,0,0,0]
x_off = [80,2100,2100,80]
y_off = [80,80,2100,2100]
; size of VVV reference image
xsize_ref = 1453
ysize_ref = 1453

; -------------------------------------------------------


      
      cube_path = '/Users/alvaromartinez/Desktop/PhD/HAWK/GNS_2/data/GNS/2021/H/' + strn(field_nr) + '/cubes/'
      data_path = '/Users/alvaromartinez/Desktop/PhD/HAWK/GNS_2/data/GNS/2021/H/' + strn(field_nr) + '/data/'
      out_path =  '/Users/alvaromartinez/Desktop/PhD/HAWK/GNS_2/data/GNS/2021/H/' + strn(field_nr) + '/aligned/'
      ;~ cube_path = '/home/data/GNS/2015/'+band+'/' + strn(field_nr) + '/cubes/'
      ;~ data_path = '/home/data/GNS/2015/'+band+'/' + strn(field_nr) + '/data/'
      ;~ out_path =  '/home/data/GNS/2015/'+band+'/' + strn(field_nr) + '/aligned/'
      pruebas= '/Users/amartinez/Desktop/PhD/HAWK/GNS_2/pruebas/'



scale = 0.34/0.106
; HAWKI --> HAWKI
;scale = 1.0



x_trans = round(xsize_ref * scale)
y_trans = round(ysize_ref * scale)



 ; Loop over four quadrants
 ; ------------------------

 for chip = chip_nr, chip_nr do begin

  chip_nr = strn(chip)

  ; Read image distortion parameters computed by:q!
  ; alignquadrants.pro
  RESTORE, data_path + 'AlignPars_chip' + strn(chip_nr)

  readcol, cube_path + 'list_chip' + chip_nr + '.txt', cube_names, FORMAT='A
  readcol, cube_path + 'masklist_chip' + chip_nr + '.txt', mask_names, FORMAT='A'
  n_cubes = n_elements(cube_names)
  

  xlo = x_off[chip-1]
  ylo = y_off[chip-1]
  xhi = x_off[chip-1] + x_out - 1
  yhi = y_off[chip-1] + y_out - 1
  
  print,xlo,xhi
  for ic = 0, n_cubes - 1 do begin

 
   masks = readfits(cube_path + mask_names[ic])
   cube = readfits(cube_path + cube_names[ic])

   sz = size(cube)
   print,sz
   n3 = sz[3]
   
   ; Note: Remove first frame in cube
   ;~ outcube = fltarr(x_out,y_out,n3-1)
   ;~ outmasks = fltarr(x_out,y_out,n3-1)
   
   ;Note: keeps the first frame:
   outcube = fltarr(x_out,y_out,n3)
   outmasks = fltarr(x_out,y_out,n3)
   
   
   trans_mask=fltarr(x_trans+200,y_trans+200)
  ;Because the masks are always the same, I put this here instead of the "for" cycle.
  ; 
   trans_mask[100:x_trans+99,0:y_trans-1] = POLY_2D(masks,Kx,Ky,2,x_trans,y_trans,CUBIC=-0.5,MISSING=0)
   mask=trans_mask
   ;~ mask = POLY_2D(masks,Kx,Ky,2,x_trans-200,y_trans-200,CUBIC=-0.5,MISSING=0)
   
   bad = where(mask lt 1)
   mask[bad] = 0
   
   for j = 0, n3 -1 do begin ; Keep the first frame in cube
   ;~ for j = 1, n3 -1 do begin ; Remove first frame in cube
    trans_im=fltarr(x_trans+200,y_trans+200)
    im = cube[*,*,j]
    ;~ writefits, pruebas + 'im_test_chip'+chip_nr+'_'+strn(ic)+'.fits',im;·································
    trans_im[100:x_trans+99,0:y_trans-1] = POLY_2D(im,Kx,Ky,2,x_trans,y_trans,CUBIC=-0.5,MISSING=0)
    im=trans_im
    ;~ im = POLY_2D(im,Kx,Ky,2,x_trans-200,y_trans-200,CUBIC=-0.5,MISSING=0)
    
    im = im * mask
    ;~ writefits, pruebas + 'imKK_test_chip'+chip_nr+'_'+strn([ic]+1)+'.fits',im;··································
    ;~ print, 'imKK_test_chip'+chip_nr+'_'+strn([ic]+1)+'.fits'
    ;~ print,Kx[0,0]
    ;~ stop
    
    outcube[*,*,j-1] = im[xlo:xhi,ylo:yhi]
    outmasks[*,*,j-1] = mask[xlo:xhi,ylo:yhi]
    print, j
   endfor
   
;    writefits, out_path + mask_names[ic], outmasks
;    writefits, out_path + cube_names[ic], outcube
   print, mask_names[ic]
   print,'Done Cube', ic +1
   writefits, out_path + mask_names[ic], outmasks, /COMPRESS
   writefits, out_path + cube_names[ic], outcube, /COMPRESS
 

 endfor
 ; save names of cubes and masks
 forprint, /NOCOMMENT, TEXTOUT=out_path + 'list_chip' + chip_nr + '.txt', cube_names
 forprint, /NOCOMMENT, TEXTOUT=out_path + 'masklist_chip' + chip_nr + '.txt', mask_names


endfor

END
