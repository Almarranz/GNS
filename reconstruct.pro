pro reconstruct, field_nr, chip_nr


; PURPOSE: Piece the holographically reduced sub-images together
;          Basically, this means reversing the work of sub-cubes.
;          We have to be careful, however, because holography can 
;          induce positional offsets. Also, there is differential
;          tip-tilt in the original data. Its effect will be reduced
;          with this reconstruction, too.
;
; METHOD: Use stars in the reconstructed sub-images and determine
;         their relative overall offsets between overlapping subimages
;         to correct for shifts in x and y.
;         All bands will be aligned with respect to H, which is the
;         band with the best compromise between number of sources 
;         and saturation/crowding.
;         

field = strn(field_nr)
chip = 'chip' + strn(chip_nr)

; Basic parameters
; size of final image
xaxis = 2400L
yaxis = 1200L  
sub_size_x0 = 600
sub_size_y0 = 600

rebfac = 2   ; rebin factor
sub_size_x0 = sub_size_x0 * rebfac
sub_size_y0 = sub_size_y0 * rebfac
border = 25 * rebfac  ; borders of input images will be masked
                      ; to suppress edge effects of holography

tolerance = 2 ; Correponds to dmax in search for 
              ; common stars between a subimmages
              ;  for fine-alignment of subimages
              ; with rebfac = 2, tolerance = 2
              ; corresponds to half the resolution limit

; set bands = 1 for simultaneous reconstruction of JHK
; otherwise, only H will be reconstructed
bands = 0
bands = 1

; Numbers of first and last sub-images (usually not to be changed, see
; runholo.pro)
i_x0 = 0
i_x1 = 6

i_y0 = 0
i_y1 = 2

; number of sub-images
;nsub = 9
nsub = 3

; input and output directories
indirk = '/data/GNS/2015/Ks/'
indirh = '/data/GNS/2015/H/'
indirj = '/data/GNS/2015/J/'
 
outdirh = indirh
outdirj = indirj
outdirk = indirk
 
pathh = indirh + field + '/cubeims_sub/chip' + strn(chip_nr) + '/'
pathj = indirj + field + '/cubeims/chip' + strn(chip_nr) + '/'
pathk = indirk + field + '/cubeims/chip' + strn(chip_nr) + '/'

; General StarFinder settings
estim_bg = 0
back_box = 0
correl_mag = 2.0
deblend = 0
deblost = 0
niter = 2
rel_thresh = 1
GUIDE_X = ''
GUIDE_Y = ''
; ***** TEST EFFECT OF FOLLOWING PARAMETERS HERE *****
; It may make sense to go down to Threshold = [5.]
; Also, Threshold = [5.,5.] (two iterations) may make sense
; ... at the cost of computing time.
Threshold = [20.]
 ; Lower values may make sense
 ; down to min_corrrelation = 0.7 (StarFinder default)
 ; lower correlation values give larger number of stars
 ; Higher values give higher quality positions
 ; and sepped up computing time
min_correlation = 0.9
;Rainerś preliminary assessment on 29/1/2021
; ==> With threshold = [10.] and min_correlation = 0.9 the standard
; deviations of the computed offsets appear to be <=0.0035 pixels 
; i.e. smaller than 0.2 mas. Good enough for proper motions.
; With threshold = [5.] and min_correlation = 0.9  the code is notably
; slower, but the alignment uncertainty does not appear to improve.
; It apepars to be best to use high SNR stars.
; I get very good results for Threshold = [100.] and min_correlation =
; 0.9.
; WARNINGL I have ONLY checked the uncertainties, it can be that the
; mean shifts show additional systematic effects.


; variables that save all relative offsets
x_off_tot = []
y_off_tot = []

z_tmp = []

thissigj =  fltarr(xaxis*rebfac,yaxis*rebfac)
thisexpj =  fltarr(xaxis*rebfac,yaxis*rebfac)
thisimj =  fltarr(xaxis*rebfac,yaxis*rebfac)

tmp_subfield = fltarr(xaxis*rebfac,yaxis*rebfac)
tmp_subfield_exp = fltarr(xaxis*rebfac,yaxis*rebfac)
tmp_subfield_noise = fltarr(xaxis*rebfac,yaxis*rebfac)
;tmp_subfield_wt = fltarr(xaxis*rebfac,yaxis*rebfac) 
     
 
; Define simple Gaussian PSF for rapid source detection
; with StarFinder
psf = psf_gaussian(NPIXEL=20,FWHM=2*rebfac,/NORMALIZE,/DOUBLE)
;writefits, 'psf.fits', psf

; STEP 1: ALIGN J NAD K SUB-IMAGES WITH THEIR H COUNTERPARTS
; ------------------------------------------------------------

for i_x = i_x0, i_x1 do begin
  for i_y = i_y0, i_y1 do begin 

     name = 'holo_' + strn(i_x) + '_' +  strn(i_y)

     ;Mask borders of each input sub-image
     ; because borders are affected by  ringing (caused by division in 
     ; Fourier space in holography code)

;     mask_borders = fltarr(sub_size_x0, sub_size_y0)
;     mask_borders[*,*] = 0
;     mask_borders[border:sub_size_x0-border-1, border:sub_size_y0-border-1] = 1


; 1) Read input subimages
; mask the borders
; -----------------------

     h = readfits(pathh + name + '.fits.gz'); * mask_borders
     h_exp = readfits(pathh + name + '_wt.fits.gz');  * mask_borders    
     hsub = readfits(pathh + name + '_s1.fits.gz'); * mask_borders
     for i_s = 2, nsub do begin
       tmpim = readfits(pathh + name + '_s' + strn(i_s) + '.fits.gz'); * mask_borders
       hsub = [[[hsub]],[[tmpim]]]
     endfor
     
     h_noise = readfits(pathh + name + '_sigma.fits.gz'); * mask_borders
     h_noisesub = readfits(pathh + name + '_sigma_s1.fits.gz'); * mask_borders
     for i_s = 2, nsub do begin
       tmpim = readfits(pathh + name + '_sigma_s' + strn(i_s) + '.fits.gz'); * mask_borders
       h_noisesub = [[[h_noisesub]],[[tmpim]]]
     endfor


     ; Read J and Ks input if bands = 1
     ; Borders will only be masked after alignment with H-band images
     ; ------------------------------------------------------------
     if bands eq 1 then begin
  
       j = readfits(pathj + name + '.fits.gz')
       j_exp = readfits(pathj + name + '_wt.fits.gz')
       jsub = readfits(pathj + name + '_s1.fits.gz')
       for i_s = 2, nsub do begin
         tmpim = readfits(pathj + name + '_s' + strn(i_s) + '.fits.gz')
         jsub = [[[jsub]],[[tmpim]]]
       endfor
     
       j_noise = readfits(pathj + name + '_sigma.fits.gz')
       j_noisesub = readfits(pathj + name + '_sigma_s1.fits.gz')
       for i_s = 2, nsub do begin
         tmpim = readfits(pathj + name + '_sigma_s' + strn(i_s) + '.fits.gz')
         j_noisesub = [[[j_noisesub]],[[tmpim]]]
       endfor

       k = readfits(pathk + name + '.fits.gz')
       k_exp = readfits(pathk + name + '_wt.fits.gz')
       ksub = readfits(pathk + name + '_s1.fits.gz')
       for i_s = 2, nsub do begin
         tmpim = readfits(pathk + name + '_s' + strn(i_s) + '.fits.gz')
         ksub = [[[ksub]],[[tmpim]]]
       endfor
     
       k_noise = readfits(pathk + name + '_sigma.fits.gz')
       k_noisesub = readfits(pathk + name + '_sigma_s1.fits.gz')
       for i_s = 2, nsub do begin
         tmpim = readfits(pathk + name + '_sigma_s' + strn(i_s) + '.fits.gz')
         k_noisesub = [[[k_noisesub]],[[tmpim]]]
       endfor
     
     endif    
     
  ; 2) Align with H-band via identifying and matching common stars
  ; ---------------------------------------------------------------
       starfinder, h, psf, X_BAD=x_bad, Y_BAD = y_bad, $
        BACKGROUND = background, BACK_BOX = back_box, $
        threshold, REL_THRESHOLD = rel_thresh, /PRE_SMOOTH, $
        NOISE_STD = h_noise, min_correlation, $
        CORREL_MAG = correl_mag, INTERP_TYPE = 'I', $
        ESTIMATE_BG = estim_bg, DEBLEND = deblend, N_ITER = niter, SILENT=0, $
        GUIDE_X = guide_x, GUIDE_Y = guide_y, $
        SV_SIGMA_R = 0.0085, SV_SIGMA_A= 0.0050, $
        xh, yh, fh, sx, sy, sf, c, STARS = stars, $
        LOGFILE = logfilename, /CUBIC
   
     if bands eq 1 then begin        
       starfinder, j, psf, X_BAD=x_bad, Y_BAD = y_bad, $
        BACKGROUND = background, BACK_BOX = back_box, $
        threshold, REL_THRESHOLD = rel_thresh, /PRE_SMOOTH, $
        NOISE_STD = j_noise, min_correlation, $
        CORREL_MAG = correl_mag, INTERP_TYPE = 'I', $
        ESTIMATE_BG = estim_bg, DEBLEND = deblend, N_ITER = niter, SILENT=0, $
        GUIDE_X = guide_x, GUIDE_Y = guide_y, $
        SV_SIGMA_R = 0.0085, SV_SIGMA_A= 0.0050, $
        xj, yj, fj, sx, sy, sf, c, STARS = stars, $
        LOGFILE = logfilename, /CUBIC;, /NO_SLANT
        
        starfinder, k, psf, X_BAD=x_bad, Y_BAD = y_bad, $
        BACKGROUND = background, BACK_BOX = back_box, $
        threshold, REL_THRESHOLD = rel_thresh, /PRE_SMOOTH, $
        NOISE_STD = k_noise, min_correlation, $
        CORREL_MAG = correl_mag, INTERP_TYPE = 'I', $
        ESTIMATE_BG = estim_bg, DEBLEND = deblend, N_ITER = niter, SILENT=0, $
        GUIDE_X = guide_x, GUIDE_Y = guide_y, $
        SV_SIGMA_R = 0.0085, SV_SIGMA_A= 0.0050, $
        xk, yk, fk, sx, sy, sf, c, STARS = stars, $
        LOGFILE = logfilename, /CUBIC;, /NO_SLANT
        
        
        ; Compare JHKs lists to find common stars
        ; Correct relative shifts in X
        ; and Y vs H-band positions
        ; ---------------------------------------

        compare_lists, xh, yh, xj, yj, xhjc, yhjc, xjc, yjc, MAX_DISTANCE=tolerance, SUBSCRIPTS_1=subch, $
         SUBSCRIPTS_2 = subcj, SUB1 = onlyh, SUB2 = onlyj

        compare_lists, xh, yh, xk, yk, xhkc, yhkc, xkc, ykc, MAX_DISTANCE=tolerance, SUBSCRIPTS_1=subchk,$
         SUBSCRIPTS_2 = subck, SUB1 = onlyhk, SUB2 = onlyk
     
        xj_off = median(xhjc - xjc)
        yj_off = median(yhjc - yjc)
        xj = xj + xj_off
        yj = yj + yj_off

        xk_off = median(xhkc - xkc) 
        yk_off = median(yhkc - ykc)
        xk = xk + xk_off
        yk = yk + yk_off  
      
        j = image_shift(j, xj_off, yj_off); * mask_borders
        for i_s = 0, nsub-1 do begin
         tmpim = jsub[*,*,i_s]
         tmpim = image_shift(tmpim, xj_off, yj_off); * mask_borders
         jsub[*,*,i_s] = tmpim
        endfor
       
        k = image_shift(k, xk_off, yk_off); * mask_borders
        for i_s = 0, nsub-1 do begin
         tmpim = ksub[*,*,i_s]
         tmpim = image_shift(tmpim, xj_off, yj_off); * mask_borders
         ksub[*,*,i_s] = tmpim
        endfor

        j_noise = image_shift(j_noise, xj_off, yj_off); * mask_borders
        for i_s = 0, nsub-1 do begin
         tmpim = j_noisesub[*,*,i_s]
         tmpim = image_shift(tmpim, xj_off, yj_off); * mask_borders
         j_noisesub[*,*,i_s] = tmpim
        endfor
        
        k_noise = image_shift(k_noise, xk_off, yk_off); * mask_borders
        for i_s = 0, nsub-1 do begin
         tmpim = k_noisesub[*,*,i_s]
         tmpim = image_shift(tmpim, xk_off, yk_off); * mask_borders
         k_noisesub[*,*,i_s] = tmpim
        endfor                
        
        j_exp = image_shift(j_exp, xj_off, yj_off); * mask_borders
        k_exp = image_shift(k_exp, xk_off, yk_off); * mask_borders
         
        writefits, pathj +  'offH_holo_' + strn(i_x) + '_' +  strn(i_y) + '.fits.gz', j
        writefits, pathj +  'offH_holo_' + strn(i_x) + '_' +  strn(i_y) + '_wt.fits.gz', j_exp
        for i_s = 1, nsub do writefits, pathj +  'offH_holo_' + strn(i_x) + '_' +  strn(i_y) + '_s' + strn(i_s) + '.fits.gz', jsub[*,*,i_s-1]

        writefits, pathk +  'offH_holo_' + strn(i_x) + '_' +  strn(i_y) + '.fits.gz', k
        writefits, pathk +  'offH_holo_' + strn(i_x) + '_' +  strn(i_y) + '_wt.fits.gz', k_exp
        for i_s = 1, nsub do writefits, pathk +  'offH_holo_' + strn(i_x) + '_' +  strn(i_y) + '_s' + strn(i_s) + '.fits.gz', ksub[*,*,i_s-1]

        writefits, pathj + 'offH_holo_' + strn(i_x) + '_' +  strn(i_y) + '_sigma.fits.gz' , j_noise
        for i_s = 1, nsub do writefits, pathj +  'offH_holo_' + strn(i_x) + '_' +  strn(i_y) + '_sigma_s' + strn(i_s) + '.fits.gz', j_noisesub[*,*,i_s-1]

        writefits, pathk + 'offH_holo_' + strn(i_x) + '_' +  strn(i_y) + '_sigma.fits.gz' , k_noise
        for i_s = 1, nsub do writefits, pathk +  'offH_holo_' + strn(i_x) + '_' +  strn(i_y) + '_sigma_s' + strn(i_s) + '.fits.gz', k_noisesub[*,*,i_s-1]
        
                                       
     endif        
        
;3) Place subimages in correct place in the reconstructed image
;   This is similar to undoing the " chopping"  up of the field in
;   subcubes,pro, but taking into account offsets due to 
;   atmospheric efects (tip-tilt) and holography.
; -------------------------------------------------------------


  ; The following computation may  only be correct only true for 
  ; certain ratios of sub_shift/sub_size
  ; depends on settings in subcubes.pro
  ; Here, the assumption is that each half of a sub-image
  ; overlaps with one half of its neighbour (see subcibes.pro)

     x_sub_shift = sub_size_x0/2
     y_sub_shift = sub_size_y0/2
     nx = xaxis/x_sub_shift - 1 
     ny = yaxis/y_sub_shift - 1
       
     ; offset of sub-image (subcube) 
     ; within the original (and also reconstructed) image
     ; as defined in subcubes.pro
     x0 =  i_x*x_sub_shift
     x1 = x0 + sub_size_x0 - 1 
     y0 =  i_y*y_sub_shift 
     y1 = y0 + sub_size_y0 - 1             


     ; Now we have to carefully remove any differential shifts
     ; between the subimages that may have been caused by the 
     ; application of the speckle holography algorithm
     ; -------------------------------------------------------

     ; add offset of sub-image to positions
     ; of detected stars within this subimage (xh, yh)
     xh = xh + x0
     yh = yh + y0            
    
     ; make list of [x,y] positions of  all stars detected in all subimages         
     z_tmp = [[xh], [yh]] 
      
           
      if (i_x eq i_x0 and i_y eq i_y0) then begin  ; initialize variable for 0 offset (lower left corner)          
          
        z_acu = z_tmp  ; accumulated list of *corrected* positions for entire H band image
        x_app = [x0]   ; list of all sub-image offsets along x for current x0
        y_app = [y0]   ; list of all sub-image offsets along y for current x0
        x_ant = x0     ; saves offset of current sub-image
                
      endif else begin
        
      ; z_tmp[*,0], z_tmp[*,1], are [x,y] positions
      ; in current sub-image
      ; z_acu[*,0], z_acu[*,1] are [x,y] all already
      ; *corrected* positions in reconstructed image 
      compare_lists, z_tmp[*,0], z_tmp[*,1], z_acu[*,0], z_acu[*,1], xc_tmp, yc_tmp, xc_acu, yc_acu,$
           MAX_DISTANCE=tolerance, SUBSCRIPTS_1=subchk, SUBSCRIPTS_2 = subck, SUB1 = onlyhk, SUB2 = onlyk
        
        
        ; align star positions of this sub-image with full list
        off_x = median(xc_acu-xc_tmp)
        off_y = median(yc_acu-yc_tmp)
        nc = n_elements(xc_tmp)
        doff_x = stddev(xc_acu-xc_tmp)/sqrt(nc)
        doff_y = stddev(yc_acu-yc_tmp)/sqrt(nc)
        print, 'Median offset in x: ' + strn(off_x) + '+- ' + strn(doff_x)
        print, 'Median offset in y: ' + strn(off_y) + '+- ' + strn(doff_y)

        z_corrected = [[z_tmp[*,0] + off_x], [z_tmp[*,1] + off_y]]
        z_acu = [z_acu, z_corrected]  ; add the corrected positions to the end of the full list 
        
        if (x0 eq x_ant) then begin
          
          x_app = [x_app, off_x + x0]
          y_app = [y_app, off_y + y0]
          
          x_ant = x0
          
        endif else begin
          
          x_off_tot = [[x_off_tot], [x_app]] ; list of all x-offsets
          y_off_tot = [[y_off_tot], [y_app]] ; list of al y-offsets
          x_app = []
          y_app = []
          
          x_app = [x0+off_x]
          y_app = [y0+off_y]

          x_ant = x0
          
        endelse
        
        if ((i_x eq i_x1) and (i_y eq i_y1)) then begin
  
          x_off_tot = [[x_off_tot], [x_app]]
          y_off_tot = [[y_off_tot], [y_app]]  
          x_app = []
          y_app = []

        endif        
    
      endelse
      
            
      ; Delete unnecessary rows in the accumulated positions
;      ; Is this really necessary? I have commented this part.
      ; WIll make the code a bit slower, but less complex.
;      if i_x gt i_x0  then begin
;        z_acu = transpose(z_acu)
;        index = where(z_acu[0,*] gt x_off_tot[0,(i_x-i_x0)-1], count)
;        z_acu = removerows(z_acu, index)
;        z_acu = transpose(z_acu)
;      endif 

    endfor
   endfor

  ; Reconstruct images
  ; --------------------

  name = 'holo_'
  suffix = ''
  build, pathh, field, chip, i_x0, i_x1, i_y0, i_y1, x_off_tot, y_off_tot, xaxis, yaxis, rebfac, name, border, suffix
  for i_s = 1, nsub do begin
    suffix = '_s' + strn(i_s)
    build, pathh, field, chip, i_x0, i_x1, i_y0, i_y1, x_off_tot, y_off_tot, xaxis, yaxis, rebfac, name, border, suffix
  endfor

  if bands eq 1 then begin
    name = 'offH_holo_'

    suffix = ''
    build, pathj, field, chip, i_x0, i_x1, i_y0, i_y1, x_off_tot, y_off_tot, xaxis, yaxis, rebfac, name, border, suffix
    for i_s = 1, nsub do begin
      suffix = '_s' + strn(i_s)
      build, pathj, field, chip, i_x0, i_x1, i_y0, i_y1, x_off_tot, y_off_tot, xaxis, yaxis, rebfac, name, border, suffix
    endfor

    suffix = ''
    build, pathk, field, chip, i_x0, i_x1, i_y0, i_y1, x_off_tot, y_off_tot, xaxis, yaxis, rebfac, name, border, suffix
    for i_s = 1, nsub do begin
      suffix = '_s' + strn(i_s)
      build, pathk, field, chip, i_x0, i_x1, i_y0, i_y1, x_off_tot, y_off_tot, xaxis, yaxis, rebfac, name, border, suffix
    endfor
  endif


end
