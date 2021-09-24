PRO RUNHOLO, field_nr, chip

; THESE PARAMETERS CAN BE EDITED
; USUALLY, THEY STAY UNCHANGED

      iter = ''  ; if '', then results of SF on long exposure is used - default
                  ; if '2', then  results of SF on holo rebin 1 is used


      band = 'H'
      ZP = 26.3 ; rough mean ZP for all HAWK-I detectors
                ; estimated from ESO QC data base
                ; USE ZP FOR CORRECT FILTER!
      DIT = 1.26   ; DIT [s] of observat1ions
      mag_min = 16  ; faintest mag for reference stars
      mag_max = 10  ; brightest mag for reference stars

      psf_path = '/home/data/GNS/2015/'+band+'/' + strn(field_nr) + '/data/'
      in_path = '/home/data/GNS/2015/'+band+'/' + strn(field_nr) + '/subcubes/'
      out_path = '/home/data/GNS/2015/'+band+'/' + strn(field_nr) + '/cubeims/'
      outpsf_path = '/home/data/GNS/2015/'+band+'/' + strn(field_nr) + '/psfs/'
      data_path = '/home/data/GNS/2015/'+band+'/' + strn(field_nr) + '/data/'
      tmp_path = '/home/data/GNS/2015/'+band+'/' + strn(field_nr) + '/tmp/'
      tmp_path_s = tmp_path +  'tmp'+strn(chip)+'/'

      x_cube = 600 ; xaxis length of sub-cube
      y_cube = 600 ; yaxis length of sub-cube
      x_large = 2400  ; xaxis length of large cube
      y_large = 1200  ; yaxis length of large cube
      
      SIGDEV = 3.0 ; standard deviations used in RESISTANT_MEAN

; Use SSA PSF as a reference for PSF size and seeing
 ; PSF should be roughly the same for all chips
 ; use PSF from long exposure, NOT holography (NOT holo2)
  ssa_psf = readfits(psf_path + 'psf_chip'+strn(chip)+'_holo.fits')
  psf_fwhm = fwhm(ssa_psf)

; Parameters for reference star selection
; ---------------------------------------
delta_r = 2. * psf_fwhm
delta_mag = 2.5 ; 2.5 mag correspond to a factor of 10
n_ref_max = 20  ; max number of reference stars to use in second iteration

; Parameters for holography
; --------------------------

debug = 1
rebfac = 2L
nsub = 9         ; number of sub-images (HOLO_MOSAIC_SUBIM.PRO) or jackknife images (HOLO_MOSAIC.PRO)
nsigma = [1.,2.] ; noise thesholds for PSF 
satlevel = 4.0e4 ; saturation threshold of data
unweighted = 0
psfnoise = 0
smoothmask = 0 ; to suppress edge effects from mask
holo_iter = n_elements(nsigma) - 1
normrad = fwhm(ssa_psf)
starlist = tmp_path_s + 'stars.txt'
psf_size = rebfac * (x_cube < y_cube)
airy = psf_gaussian(NPIXEL=psf_size,FWHM=2*rebfac,/NORMALIZE,/DOUBLE)
bord = 0  ; width of cosine shaped transition region from 0 to 3 sigma noise threshold
out_iter = 10    ; output intermediate result after processing out_iter images
psfout = 0       ; psf output?
holoout = 0
subpix = 1       ; sub-pixel alignment of PSFs, Y/N
psfavg = 0      ; set > 0 for PSF mean superposition, else median
clip = 0
minsupp = 0.9
maxsupp = 1.1
weightframes = 0

maskrad = round(5*psf_fwhm)
boxhw = round(7*psf_fwhm)   ; half width of size of PSF image for holography
psf_border = 2*psf_fwhm

;print, maskrad, boxhw
;determine number of pixels in circular aperture within maskrad/2
; needed for reference star selection
dummy = fltarr(2*boxhw+1,2*boxhw+1)
dummy[*,*] = 1
circmask = circ_mask(dummy,boxhw,boxhw,maskrad)
n_inner = total(circmask)

; large circular mask to avoid reference sources near saturated stars
dummy = fltarr(4*boxhw+1,4*boxhw+1)
dummy[*,*] = 1
circbig = circ_mask(dummy,2*boxhw,2*boxhw,2*maskrad)


x_sub_shift = x_cube/2
y_sub_shift = y_cube/2
nx = x_large/x_sub_shift - 1
ny = y_large/y_sub_shift - 1


 chip_nr = strn(chip)
 ; read positions and fluxes of stars detected on SSA image of this chip
 readcol, data_path + 'stars' + '_' + chip_nr + '_holo' + iter + '.txt', x_chip, y_chip, f_chip, correl

 for i_x = 0, nx -1 do begin
  for i_y = 0, ny -1 do begin
; for i_x = 1, 1 do begin
;  for i_y = 1, 1 do begin

   ; Compute longexposure image (for debugging) and support region
   filenam = '_' + strn(i_x) + '_' + strn(i_y)
   readcol, in_path  + 'chip' + chip_nr + '/masklist' + filenam + '.txt', mnames, FORMAT='A'
   readcol, in_path  + 'chip' + chip_nr + '/list' + filenam + '.txt', cnames, FORMAT='A'
   lxp = fltarr(x_cube,y_cube)
   support = fltarr(x_cube,y_cube)
   nm = n_elements(mnames)
   for j = 0, nm-1 do begin
    maskcube = readfits(in_path + 'chip' +  chip_nr + '/' + mnames[j])
    support = support + total(maskcube,3)
    cube = readfits(in_path + 'chip' +  chip_nr + '/' + cnames[j])
    lxp = lxp + total(cube,3)
   endfor
   accept = where(support gt 0)
   lxp[accept] = lxp[accept]/support[accept]
   support = support/max(support)
   writefits, tmp_path_s + 'lxp.fits', lxp
   writefits, tmp_path_s + 'support.fits', support

   ; Select stars in this sub-field
   xlo = i_x * x_sub_shift
   ylo = i_y * y_sub_shift   
   x_sub = x_chip - xlo
   y_sub = y_chip - ylo
   subind = where(x_sub ge 0 and x_sub le (x_cube-1) and y_sub ge 0 and y_sub le (y_cube-1), n_sub) 
   x_sub = x_sub[subind]
   y_sub = y_sub[subind]
   f_sub = f_chip[subind]
   ; save list of stars to pass it on to holo_mosaic
   forprint, TEXTOUT=starlist, /NOCOMMENT, x_sub, y_sub, f_sub
   dat = ptr_new({X_size: 20, Y_size: 20, Sigma_x: 2., Sigma_y: 2., Angle: 0.0})
   substars = image_model(x_sub,y_sub,f_sub,x_cube,y_cube,'gaussian', dat)
   writefits, tmp_path_s + 'sources.fits', substars

  ; Select PSF  reference stars among stars in selected brightness range
   ; Use a triple criterion: 
   ; (1) Minimum and maximum brightness.
   ; (2) No brighter star within delta_r.
   ; (3) Any star within delta_r  must be at least delta_mag fainter.
 ; ------------------------------------------------------------------- 

   m_sub = ZP - 2.5 * alog10(f_sub/DIT)

   ord = sort(m_sub)
   x_sub = x_sub[ord]
   y_sub = y_sub[ord]
   m_sub = m_sub[ord]

   ; select by magnitude 
   ref_ind = where(m_sub le mag_min and m_sub ge mag_max,n_ref) 
   x_psf = x_sub[ref_ind]
   y_psf = y_sub[ref_ind]
   m_psf = m_sub[ref_ind]

  ; PSF stars must be isolated
   isolated_stars, x_sub, y_sub, m_sub, x_psf, y_psf, m_psf, delta_mag, delta_r, ind_iso
   x_psf =x_psf[ind_iso]
   y_psf =y_psf[ind_iso]
   m_psf =m_psf[ind_iso]
   isolated_stars, x_sub, y_sub, m_sub, x_psf, y_psf, m_psf, 0, maskrad, ind_iso
   x_psf =x_psf[ind_iso]
   y_psf =y_psf[ind_iso]
   m_psf =m_psf[ind_iso]
   n_ref = n_elements(m_psf)
   print, 'Found '+ strn(n_ref) + ' isolated reference stars.'

  ; exclude potentially saturated reference stars and those close to saturated sources
   sat_pixels = where(lxp gt satlevel, complement=not_saturated,n_saturated)
   if (n_saturated gt 0) then begin
     sat_mask = lxp
     sat_mask[not_saturated] = 0
     sat_mask[sat_pixels] = 1
     sat_mask = CONVOLVE(sat_mask,circbig)
     goodpix = where(sat_mask lt 1,complement=maskpix)
     sat_mask[maskpix] = 0
     sat_mask[goodpix] = 1
     writefits, tmp_path_s + 'saturation_mask.fits', sat_mask
     accept = []
     for s = 0, n_ref-1 do begin
       xx = round(x_psf[s])
       yy = round(y_psf[s])
       xx = 0 > xx & xx = (x_cube-1) < xx
       yy = 0 > yy & yy = (y_cube-1) < yy
       if (sat_mask[xx,yy] gt 0) then accept = [accept,s]
     endfor
     x_psf = x_psf[accept]
     y_psf = y_psf[accept]
     m_psf = m_psf[accept]
   endif
   n_ref = n_elements(m_psf)
   print, 'Found '+ strn(n_ref) + ' isolated and unsaturated reference stars.'


   ; Sort from brightest to faintest: Important!
   ; ------------------------------------------
   ord = sort(m_psf)
   x_psf = x_psf[ord]
   y_psf = y_psf[ord]
   m_psf = m_psf[ord]
   n_ref = n_elements(m_psf)
   f_psf = 10^(0.4*(ZP-m_psf))
   print, 'Found '+ strn(n_ref) + ' reference stars.'
   dat = ptr_new({X_size: 20, Y_size: 20, Sigma_x: 2., Sigma_y: 2., Angle: 0.0})
   refstars = image_model(x_psf,y_psf,f_psf,x_cube,y_cube,'gaussian', dat)
   writefits, tmp_path_s + 'refstars.fits', refstars

  ; Run holography
  ; -------------------- 
    filenam = '_' + strn(i_x) + '_' + strn(i_y)
    indir =  in_path  + 'chip' + chip_nr + '/'
    inlist = in_path  + 'chip' + chip_nr + '/list' + filenam + '.txt'
    maskdir = in_path  + 'chip' + chip_nr + '/'
    masklist = maskdir + 'masklist' + filenam + '.txt'

    outnam = out_path + 'chip' + chip_nr + '/' + 'holo_' + strn(i_x) + '_' +  strn(i_y)
    refsources = fltarr(3,n_ref)
    refsources[0,*] = x_psf
    refsources[1,*] = y_psf
    refsources[2,*] = f_psf
    mrad = maskrad
    bhw = boxhw
    nrad = normrad
    n_mask_secondary = 0
    out_psf_path_s =  outpsf_path+ 'chip'  + chip_nr + '/'

    holo_iter = n_elements(nsigma) - 1
    ; using the saturaiton mask will not
    ; necessarily improve first estimation
    ; of PSF
    ; ESTIM_BG usually not necessary and costs execution time
    HOLO_MOSAIC, indir, inlist, 'holo',  mrad, nrad, nsigma, rebfac, refsources, starlist, maskdir=maskdir, masklist=masklist, DEBUG = debug, iter=holo_iter, AIRY=airy, OUT_ITER=out_iter, PSFOUT = psfout, UNWEIGHTED = unweighted, SUBPIX=subpix, tmpdir=tmp_path_s, BOXHW=bhw, NSUB=nsub, MINSUPP=minsupp, MAXSUPP = maxsupp, REBITER=0, RAWOUT = 0, N_MASK_SECONDARY=n_mask_secondary, PSF_FRAC = 1.0, CORRECT_SKY = 0, SMOOTHMASK = smoothmask, N_REF_MAX = n_ref_max, PR = pr, SATLEVEL=satlevel, PSF_BORDER=psf_border, ESTIM_BG = 0, SAT_MASK = 0; sat_mask

;    HOLO_MOSAIC_SUBIM, indir, inlist, 'holo',  mrad, nrad, nsigma, rebfac, refsources, starlist, maskdir=maskdir, masklist=masklist, DEBUG = debug, iter=holo_iter, AIRY=airy, OUT_ITER=out_iter, PSFOUT = psfout, UNWEIGHTED = unweighted, SUBPIX=subpix, tmpdir=tmp_path_s, BOXHW=bhw, NSUB=nsub, MINSUPP=minsupp, MAXSUPP = maxsupp, REBITER=0, RAWOUT = 0, N_MASK_SECONDARY=n_mask_secondary, PSF_FRAC = 1.0, CORRECT_SKY = 0, SMOOTHMASK = smoothmask, N_REF_MAX = n_ref_max, PR = pr, SATLEVEL=satlevel, PSF_BORDER=psf_border, ESTIM_BG = 0, SAT_MASK = 0; sat_mask

 
   ; Make mosaic and noise map for deep image 
   ; 1) Deep image
   readcol, tmp_path_s + 'holo_ims.txt', holonames, FORMAT='A'
   readcol, tmp_path_s + 'weights.txt', expnames, FORMAT='A'
   n_holo = n_elements(holonames)
   cube = fltarr(x_cube*rebfac,y_cube*rebfac,n_holo)
   cube_wt = fltarr(x_cube*rebfac,y_cube*rebfac,n_holo)
   for j = 0, n_holo-1 do begin
    im = readfits(tmp_path_s + holonames[j]  + '.fits.gz')
    wt = readfits(tmp_path_s + expnames[j] + '.fits.gz')
    bad = where(wt lt 1,complement=good)
    wt[bad] = 0
    wt[good] = 1
    cube_wt[*,*,j] = wt
    cube[*,*,j] = im * wt
   endfor
   mosaic = fltarr(x_cube*rebfac,y_cube*rebfac)
   sigma = fltarr(x_cube*rebfac,y_cube*rebfac)
   weights = fltarr(x_cube*rebfac,y_cube*rebfac)
   for x = 0, x_cube*rebfac -1 do begin
    for y = 0, y_cube*rebfac - 1 do begin
      vals = cube[x,y,*]
      wt = cube_wt[x,y,*]
      good = where(wt gt 0,count)
      if count gt 1 then begin
        RESISTANT_Mean,vals[good],SIGDEV,m,s,Num_Rej
        mosaic[x,y] = m
        sigma[x,y] = s
        weights[x,y] = total(wt)
      endif else begin
        weights[x,y] = 0
        sigma[x,y] = 0
        mosaic[x,y] = 0
      endelse
    endfor
   endfor
   writefits, outnam + '.fits', mosaic, /COMPRESS
   writefits, outnam + '_sigma.fits', sigma, /COMPRESS
   writefits, outnam + '_wt.fits', weights, /COMPRESS


   ; Make mosaic and noise map for sub-maps
   ; 2) Submaps
   ; use same weight maps for submaps than for overall image
   cube = fltarr(x_cube*rebfac,y_cube*rebfac,n_holo)
   for j_sub = 1, nsub do begin
    for j = 0, n_holo-1 do begin           
     im = readfits(tmp_path_s + holonames[j] + '_s' + strn(j_sub)  + '.fits.gz')
;     im = readfits(tmp_path_s + holonames[j] + '_s' + strn(j_sub)  + '.fits')
     wt = cube_wt[*,*,j]
     bad = where(wt lt 1,complement=good)
     wt[bad] = 0
     wt[good]= 1
     cube[*,*,j] = im * wt
    endfor
    mosaic = fltarr(x_cube*rebfac,y_cube*rebfac)
    sigma = fltarr(x_cube*rebfac,y_cube*rebfac)
    for x = 0, x_cube*rebfac -1 do begin
     for y = 0, y_cube*rebfac - 1 do begin
       vals = cube[x,y,*]
       wt = cube_wt[x,y,*]
       good = where(wt gt 0,count)
       if count gt 1 then begin
         RESISTANT_Mean,vals[good],SIGDEV,m,s,Num_Rej
         mosaic[x,y] = m
         sigma[x,y] = s
       endif else begin
         sigma[x,y] = 0
         mosaic[x,y] = 0
       endelse
     endfor
    endfor
    writefits, outnam + '_s' + strn(j_sub) +'.fits', mosaic, /COMPRESS
    writefits, outnam + '_sigma_s' + strn(j_sub) +'.fits', sigma, /COMPRESS
  endfor

  
   ; Delete all temporary files
   print, 'Finished sub-field ' + strn(i_x) + ', ' + strn(i_y)

   spawn, 'rm ' + tmp_path_s + '*.fits.gz' 
;   spawn, 'rm ' + tmp_path_s + '*.fits' 

  endfor
 endfor

 print, 'Finished chip ' + chip_nr

END
