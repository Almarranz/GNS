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

x_out = 2400
y_out = 1200

; offsets of fields within transformed frame (see alignquadrants.pro)
; used to cut out the chips from the large algined field
x_off = [0,2200,2200,0]
y_off = [50,50,900,900]

; size of VVV reference image
xsize_ref = 1453
ysize_ref = 655

; -------------------------------------------------------


      band = 'H'
      cube_path = '/home/data/GNS/2015/'+band+'/' + strn(field_nr) + '/cubes/'
      data_path = '/home/data/GNS/2015/'+band+'/' + strn(field_nr) + '/data/'
      out_path =  '/home/data/GNS/2015/'+band+'/' + strn(field_nr) + '/aligned/'
      


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

  readcol, cube_path + 'list_chip' + chip_nr + '.txt', cube_names, FORMAT='A'
  readcol, cube_path + 'masklist_chip' + chip_nr + '.txt', mask_names, FORMAT='A'
  n_cubes = n_elements(cube_names)

  xlo = x_off[chip-1]
  ylo = y_off[chip-1]
  xhi = x_off[chip-1] + x_out - 1
  yhi = y_off[chip-1] + y_out - 1

  for ic = 0, n_cubes - 1 do begin

 
   masks = readfits(cube_path + mask_names[ic])
   cube = readfits(cube_path + cube_names[ic])

   sz = size(cube)
   n3 = sz[3]
   ; Note: Remove first frame in cube
   outcube = fltarr(x_out,y_out,n3-1)
   outmasks = fltarr(x_out,y_out,n3-1)
  
  ;Because the masks are always the same, I put this here instead of the "for" cycle.
  ; 
   mask = POLY_2D(masks,Kx,Ky,2,x_trans,y_trans,CUBIC=-0.5,MISSING=0)
   bad = where(mask lt 1)
   mask[bad] = 0
   
   for j = 1, n3 -1 do begin ; Remove first frame in cube
    im = cube[*,*,j]

    im = POLY_2D(im,Kx,Ky,2,x_trans,y_trans,CUBIC=-0.5,MISSING=0)
    im = im * mask
    outcube[*,*,j-1] = im[xlo:xhi,ylo:yhi]
    outmasks[*,*,j-1] = mask[xlo:xhi,ylo:yhi]
    print, j
   endfor

   writefits, out_path + mask_names[ic], outmasks, /COMPRESS
   writefits, out_path + cube_names[ic], outcube, /COMPRESS
 

 endfor
 ; save names of cubes and masks
 forprint, /NOCOMMENT, TEXTOUT=out_path + 'list_chip' + chip_nr + '.txt', cube_names
 forprint, /NOCOMMENT, TEXTOUT=out_path + 'masklist_chip' + chip_nr + '.txt', mask_names


endfor

END
