pro merge_sublists, field_nr, chip_nr


; PURPOSE: Piece the lists of the holographically reduced sub-images together
;          use input from align_lists.pro
;
; METHOD: All sub-fields will be aligned progressively,
;         starting with the left lower sub-field and aligning
;         new lists relative to the currently existing
;         full list
; 
;         Finally, all multiple measurements will be
;         averaged, with the standard deviation computed as the "PSF"
;         uncertainty.
;         The PSF uncertainty will 
;         be the median PSF uncertainty of stars within the same
;         magnitude range.      
;
;         Stars near to the borders of sub-images will be deselected 
;         to avoid bias due to edge effects.

Band = 'H'
field = strn(field_nr)
chip = 'chip' + strn(chip_nr)
t_exp = 3.32
if Band eq 'J' then ZP = 26.3
if Band eq 'H' then ZP = 26.3
if Band eq 'Ks' then ZP = 25.6

; Set minimum brightness for stars
; that are used for relative calibration
; of the sub-fields
if Band eq 'J' then mag_cal = 18
if Band eq 'H' then mag_cal = 16
if Band eq 'Ks' then mag_cal = 14

; input and output directories
basedir = '/Users/alvaromartinez/Desktop/Phd/HAWK/GNS_2/Data/GNS/2021/'
plotdir = basedir + 'H/6/photo/chip4/plots/'
indir = '/Users/alvaromartinez/Desktop/Phd/HAWK/GNS_2/Data/GNS/2021/'+ Band + '/' + field
outdir = indir
photdir = indir + '/photo/chip' + strn(chip_nr) + '/lists/'
imdir = indir + '/cubeims/chip' + strn(chip_nr) + '/'

; Basic parameters
; size of final image
xaxis = 2700L
yaxis = 2700L  
sub_size_x0 = 1350
sub_size_y0 = 1350
edge = 20


rebfac = 2   ; rebin factor
edge =  edge * rebfac
xaxis = xaxis * rebfac
yaxis = yaxis * rebfac
sub_size_x0 = sub_size_x0 * rebfac
sub_size_y0 = sub_size_y0 * rebfac
xlo = edge-1
xhi = sub_size_x0 - edge
ylo = edge
yhi = sub_size_y0 - edge
; Here, the assumption is that each half of a sub-image
; overlaps with one half of its neighbour (see subcubes.pro)
x_sub_shift = sub_size_x0/2
y_sub_shift = sub_size_y0/2
nx = xaxis/x_sub_shift - 1 
ny = yaxis/y_sub_shift - 1
       
n_multi = 4 ; maximum number of measurements for a given source (four overlapping sub-fields)
; tolerance = 2.0 
tolerance = 10.0 ; Correponds to dmax in search for 
              ; common stars between a subimmages
              ;  for fine-alignment of subimages
              ; with rebfac = 2, tolerance = 2
              ; corresponds to half the resolution limit
              ; For S/N = 10 position uncertainty accoridng to statistics
              ; shoudl be FWHM/10 = 0.4 pixels
              ; With tolerance = 2 we get some multiple detections at ks



; (1) determine the max size of the arrays containing the measurements
; ----------------------------------------------------------------------
nmax = 0L
for i_x = 0, nx-1 do begin
  for i_y = 0, ny-1 do begin 

     name = 'list_' + strn(i_x) + '_' +  strn(i_y)
     ; for J and Ks use lists pre-aligned with H
     ; that are produced by align_fields.pro
     if (Band ne 'H') then name = 'offH_' + name    
     
     readcol, photdir + name + '.txt', x, y, f, sx, sy, sf, c, /SILENT
     nmax = nmax + n_elements(f)
     print, 'nmax, i_x, i_y',nmax, i_x, i_y
  endfor
endfor
print, 'Max number of sources: ' + strn(nmax)
n_all = fltarr(nmax)
x_all = fltarr(n_multi, nmax)
y_all = fltarr(n_multi, nmax)
f_all = fltarr(n_multi, nmax)
sx_all = fltarr(n_multi, nmax)
sy_all = fltarr(n_multi, nmax)
sf_all = fltarr(n_multi, nmax)
c_all = fltarr(n_multi, nmax)

; (2) Read the data, determine astro offsets,
; correct them, and store all measurements in arrays
; -------------------------------------------------

n_detected = 0L ; number of stars already in combined list
f_offsets = fltarr(nx,ny) ; used just for testing
if (Band eq 'H') then begin
  x_offsets = fltarr(nx,ny)
  y_offsets = fltarr(nx,ny)
  f_offsets = fltarr(nx,ny)
  f_offsets[0,0] = 1.0
;  RESTORE, photdir + 'f_cal_sub.SAV'
endif else begin
  RESTORE, 'x_offsets.SAV'
  RESTORE, 'y_offsets.SAV'
endelse

for i_x = 0, nx-1 do begin
  for i_y = 0, ny-1 do begin 

     name = 'list_' + strn(i_x) + '_' +  strn(i_y)
     ; for J and Ks use lists pre-aligned with H
     ; that are produced by align_fields.pro
     if (Band ne 'H') then name = 'offH_' + name    

     readcol, photdir + name + '.txt', x, y, f, sx, sy, sf, c, /SILENT

     ; Deselect stars near edge of image     
     in_frame = where((x gt xlo) and (x lt xhi) and (y gt ylo)  and (y lt yhi))
     x = x[in_frame]
     y = y[in_frame]
     f = f[in_frame]
     sx = sx[in_frame]
     sy = sy[in_frame]
     sf = sf[in_frame]
     c = c[in_frame]

    ; Apply pre-offsets for this sub-field
    ; ----------------------------------------
     ; offset of sub-image
     ; within the original image
     ; as defined in subcubes.pro
     xoff_0 =  i_x * x_sub_shift/2
     yoff_0 =  i_y * y_sub_shift /2
     x = x + xoff_0
     y = y + yoff_0
    print,'xoff_0,y_off0',xoff_0,yoff_0
    ; compare list for this sub-field with current list of 
    ; all detected stars
    ; -----------------------------------------------------
    if (n_detected gt 0) then begin
      x1 = xl
      y1 = yl
    endif else begin
      x1 = [0.0]
      y1 = [0.0]
    endelse
    x2 = x
    y2 = y

    compare_lists, x1, y1, x2, y2, x1c, y1c, x2c, y2c, MAX_DISTANCE=tolerance, SUBSCRIPTS_1=subc1, SUBSCRIPTS_2 = subc2, SUB1 = sub1, SUB2 = sub2

    nc = n_elements(subc2)
    if (i_x eq 0) and (i_y eq 0) then nc = 0L
    print, 'Number of common stars: ' + string(nc)
    if (sub1[0] gt (-1)) then n_new = n_elements(sub2)
    print, 'Number of new stars: ' + string(n_new)

    ; Now we have to carefully remove any differential shifts
    ; between the subimages that may have been caused by the 
    ; application of the speckle holography algorithm. 
    ; Ks and J have already been aligned with H. 
    ; Therefore we compute the offsets only for H.
    ;
    ; We also correct for flux offsets.
    ; -------------------------------------------------------
    
    
    if (i_x gt 0) or (i_y gt 0) then begin
;    if (i_x gt 0)  then begin
    
      if (nc lt 1) then begin
        print, 'No common stars, cannot compute offsets.'
        STOP
      endif

      print
      f_cal = f[subc2]
      mags = ZP-2.5*alog10(f_cal)
      cal_ind = where(mags lt mag_cal,n_cal)
      print, 'Field: ' +  strn(i_x) + ', ' + strn(i_y) + ':'
      print, 'Number of calibration stars:  ' + strn(n_cal)

      if (Band eq 'H') then begin
        x_off = median(x1c[cal_ind] - x2c[cal_ind])
        y_off = median(y1c[cal_ind] - y2c[cal_ind])
        dx_off = stddev(x1c[cal_ind] - x2c[cal_ind])/sqrt(n_cal)
        dy_off = stddev(y1c[cal_ind] - y2c[cal_ind])/sqrt(n_cal)
        x_offsets[i_x,i_y] = x_off
        y_offsets[i_x,i_y] = y_off
        print, 'Offsets in x, y:'
        print, strn(x_off) + ', ' + strn(y_off)
        print, 'Uncertainty of offsets in x, y:'
        print, strn(dx_off) + ', ' + strn(dy_off)
      endif else begin
        x_off =  x_offsets[i_x,i_y]
        y_off =  y_offsets[i_x,i_y]
      endelse
  
 ;     f_off = median(fl[subc1[cal_ind]]/f[subc2[cal_ind]])
 ;     df_off = stddev(fl[subc1[cal_ind]]/f[subc2[cal_ind]])/sqrt(n_cal)
 ;     f_offsets[i_x,i_y] = f_off
;      cghistoplot, (fl[subc1[cal_ind]]/f[subc2[cal_ind]])
;      f = f * f_off ; Scaling the flux does not improve things,
;      rather to the contrary.
;      It appears that flux variations are very small <1%, but corrections
;      will pull it all into a particular direction. There can be the
;      odd outlier with f_off = 0.96 

      x = x + x_off
      y = y + y_off
      f_off = f_offsets[i_x,i_y]
;     Uncomment if you wish to photometrically cross-calibrate the holography subfields
;      f = f*f_off
      print, 'Flux calibration factor:' + strn(f_off)

    endif

  ; save values of new common stars in common list
  ; i.e. create additional entries (measurements) for
  ; stars already in the common list
  ; ------------------------------------------------

    for jj = 0L, (nc-1) do begin
      n_all[subc1[jj]] = n_all[subc1[jj]] + 1.0
      ind = n_all[subc1[jj]] - 1
      x_all[ind,subc1[jj]] = x[subc2[jj]]
      y_all[ind,subc1[jj]] = y[subc2[jj]]
      f_all[ind,subc1[jj]] = f[subc2[jj]]
      sx_all[ind,subc1[jj]] = sx[subc2[jj]]
      sy_all[ind,subc1[jj]] = sy[subc2[jj]]
      sf_all[ind,subc1[jj]] = sf[subc2[jj]]
      c_all[ind,subc1[jj]] = c[subc2[jj]]
    endfor

  ; save values of new stars in common list
  ; i.e. create entries for new stars
  ; -----------------------------------------

   for jj = 0L, (n_new-1) do begin
      n_all[(n_detected+jj)] = 1.0
      x_all[0,(n_detected+jj)] = x[sub2[jj]]
      y_all[0,(n_detected+jj)] = y[sub2[jj]]
      f_all[0,(n_detected+jj)] = f[sub2[jj]]
      sx_all[0,(n_detected+jj)] = sx[sub2[jj]]
      sy_all[0,(n_detected+jj)] = sy[sub2[jj]]
      sf_all[0,(n_detected+jj)] = sf[sub2[jj]]
      c_all[0,(n_detected+jj)] = c[sub2[jj]]
    endfor

    ; Now create a list with the mean positions
    ; of already detected stars
    n_detected = n_detected + n_new ; current total number of detected stars
    nl = n_all[0:n_detected-1]
    xl = fltarr(n_detected)
    yl = fltarr(n_detected)
    fl = fltarr(n_detected)
    for i_star = 0, n_detected-1 do begin    
      vals = f_all[*,i_star]
      good = where(vals gt 0)
      fl[i_star] = mean(vals[good])
      vals = x_all[*,i_star]
      xl[i_star] = mean(vals[good])
      vals = y_all[*,i_star]
      yl[i_star] = mean(vals[good])
    endfor

  endfor
endfor

; save offsets
; -------------
writefits, photdir + 'flux_offsets_simple.fits', f_offsets

if (Band eq 'H') then begin
  SAVE, x_offsets, FILENAME='x_offsets.SAV'
  SAVE, y_offsets, FILENAME='y_offsets.SAV'
endif

        
;3) Compute mean values for all detected stars and measurements
; --------------------------------------------------------------

print, 'Detected ' + strn(n_detected) + ' stars.'
nl = n_all[0:n_detected-1]
xl = fltarr(n_detected)
yl = fltarr(n_detected)
fl = fltarr(n_detected)
cl = fltarr(n_detected)
sxl = fltarr(n_detected) ; random uncertainty
syl = fltarr(n_detected)
sfl = fltarr(n_detected)
dxl = fltarr(n_detected) ; PSF uncertainty
dyl = fltarr(n_detected)
dfl = fltarr(n_detected)

for i_star = 0, n_detected-1 do begin    

  vals = f_all[*,i_star]
  good = where(vals gt 0,n_good)

  vals = x_all[*,i_star]
  xl[i_star] = mean(vals[good])
  if (n_good gt 1) then dxl[i_star] = stddev(vals[good])/sqrt(n_good)
  vals = y_all[*,i_star]
  yl[i_star] = mean(vals[good])
  if (n_good gt 1) then dyl[i_star] = stddev(vals[good])/sqrt(n_good)
  vals = f_all[*,i_star]
  fl[i_star] = mean(vals[good])
  if (n_good gt 1) then dfl[i_star] = stddev(vals[good])/sqrt(n_good)
  vals = sx_all[*,i_star]
  sxl[i_star] = mean(vals[good])
  vals = sy_all[*,i_star]
  syl[i_star] = mean(vals[good])
  vals = sf_all[*,i_star]
  sfl[i_star] = mean(vals[good])
  vals = c_all[*,i_star]
  cl[i_star] = mean(vals[good])
endfor

; 4) Determine the PSF uncertainty
; The basic assumptino is that the PSF uncertainty is not a 
; function of magnitude not of position of the sub-field.
; It can be easiest measured by using the brightest stars
; -------------------------------------------------------------------


; Assign PSF uncertainties to stars without multiple observations
; Plot systematic and statistical uncerrtainties
; ----------------------------------------------------------------

m = ZP - 2.5*alog10(fl/t_exp)
multi_obs = where(nl gt 1,complement=single_obs)
psf_cal_ind = where(nl gt 1 and m lt mag_cal)

; median relative PSF error
psferr_f = median(dfl[psf_cal_ind]/fl[psf_cal_ind])
print, 'Median PSF flux error: ' + strn(psferr_f)
psferr_x = median(dxl[psf_cal_ind])
print, 'Median PSF x position error: ' + strn(psferr_x)
psferr_y = median(dyl[psf_cal_ind])
print, 'Median PSF y position error: ' + strn(psferr_x)

; assign individual PSF uncertainties
; for plots
dfl[single_obs] = fl[single_obs] * psferr_f
dxl[single_obs]  = median(dxl[psf_cal_ind])
dyl[single_obs] = median(dyl[psf_cal_ind])


dm = 2.5/alog(10) *  dfl/fl
sm = 2.5/alog(10) *  sfl/fl

; Make some control plots
; -========================

blue = cgColor('Dodger Blue')

; X-uncertainties
; --------------
set_plot,'PS',/interpolate
device, XOFFSET=0, YOFFSET=0, $
      FILENAME = plotdir + 'sx_' +  chip + '_' + band + '_' + '.eps', XSIZE=15., YSIZE=15., $
      /portrait, /color,BITS_PER_PIXEL=8, encapsulated=1
!P.MULTI=[1,1,0]
!P.CHARSIZE=1.2
!P.THICK=2.0
!P.CHARTHICK=2.0
cgplot, m, sxl,  PSYM=3, Color = 'black', XTITLE='[' + band + ']', YTITLE='STatistical uncertainty in x [pixels]',$
      THICK = 0.5, XRANGE = [10,22], YRANGE = [0,0.5]
device, /close 

set_plot,'PS',/interpolate
device, XOFFSET=0, YOFFSET=0, $
      FILENAME = plotdir + 'dx_' +  chip + '_' + band + '_' + '.eps', XSIZE=15., YSIZE=15., $
      /portrait, /color,BITS_PER_PIXEL=8, encapsulated=1
!P.MULTI=[1,1,0]
!P.CHARSIZE=1.2
!P.THICK=2.0
!P.CHARTHICK=2.0
cgplot, m, dxl,  PSYM=3, Color = 'black', XTITLE='[' + band + ']', YTITLE='PSF uncertainty in x [pixels]',$
      THICK = 0.5, XRANGE = [10,22], YRANGE = [0,0.5]
device, /close 

; Y-uncertainties
; --------------
set_plot,'PS',/interpolate
device, XOFFSET=0, YOFFSET=0, $
      FILENAME = plotdir + 'sy_' +  chip + '_' + band + '_' + '.eps', XSIZE=15., YSIZE=15., $
      /portrait, /color,BITS_PER_PIXEL=8, encapsulated=1
!P.MULTI=[1,1,0]
!P.CHARSIZE=1.2
!P.THICK=2.0
!P.CHARTHICK=2.0
cgplot, m, syl,  PSYM=3, Color = 'black', XTITLE='[' + band + ']', YTITLE='Statistical uncertainty in y [pixels]',$
      THICK = 0.5, XRANGE = [10,22], YRANGE = [0,0.5]
device, /close 

set_plot,'PS',/interpolate
device, XOFFSET=0, YOFFSET=0, $
      FILENAME = plotdir + 'dy_' +  chip + '_' + band + '_' + '.eps', XSIZE=15., YSIZE=15., $
      /portrait, /color,BITS_PER_PIXEL=8, encapsulated=1
!P.MULTI=[1,1,0]
!P.CHARSIZE=1.2
!P.THICK=2.0
!P.CHARTHICK=2.0
cgplot, m, dyl,  PSYM=3, Color = 'black', XTITLE='[' + band + ']', YTITLE='PSF uncertainty in y [pixels]',$
      THICK = 0.5, XRANGE = [10,22], YRANGE = [0,0.5]
device, /close 

; Photometric uncertainty
; -----------------------
set_plot,'PS',/interpolate
device, XOFFSET=0, YOFFSET=0, $
      FILENAME = plotdir + 'sm_' +  chip + '_' + band + '_' + '.eps', XSIZE=15., YSIZE=15., $
      /portrait, /color,BITS_PER_PIXEL=8, encapsulated=1
!P.MULTI=[1,1,0]
!P.CHARSIZE=1.2
!P.THICK=2.0
!P.CHARTHICK=2.0
cgplot, m, sm,  PSYM=3, Color = 'black', XTITLE='[' + band + ']', YTITLE='Statistical photometric uncertainty [mag]',$
      THICK = 0.5, XRANGE = [10,22], YRANGE = [0,0.5]
device, /close 

set_plot,'PS',/interpolate
device, XOFFSET=0, YOFFSET=0, $
      FILENAME = plotdir + 'dm_' +  chip + '_' + band + '_' + '.eps', XSIZE=15., YSIZE=15., $
      /portrait, /color,BITS_PER_PIXEL=8, encapsulated=1
!P.MULTI=[1,1,0]
!P.CHARSIZE=1.2
!P.THICK=2.0
!P.CHARTHICK=2.0
cgplot, m, dm,  PSYM=3, Color = 'black', XTITLE='[' + band + ']', YTITLE='PSF photometric uncertainty [mag]',$
      THICK = 0.5, XRANGE = [10,22], YRANGE = [0,0.2]
xyouts, 11, 0.15, 'PSF uncertainty: ' + strn(psferr_f)
device, /close 
SAVE, FILENAME = plotdir + 'mags_cal.SAV', m, dm

;plotm = m[multi_obs]
;dplotm = dm[multi_obs]
;SAVE, plotm, FILENAME='plotm_cal.SAV'
;SAVE, dplotm, FILENAME='dplotm_cal.SAV'

;set_plot,'PS',/interpolate
;device, XOFFSET=0, YOFFSET=0, $
;      FILENAME = 'sm_vs_dm_' +  chip + '_' + band + '_' + strn(mag_cal) + '.eps', XSIZE=15., YSIZE=15., $
;      /portrait, /color,BITS_PER_PIXEL=8, encapsulated=1
;!P.MULTI=[1,1,0]
;!P.CHARSIZE=1.2
;!P.THICK=2.0
;!P.CHARTHICK=2.0
;cgplot, sm[multi_obs], dm[multi_obs],  PSYM=3, Color = 'black', XTITLE='jackknife uncertainty [mag]', YTITLE='multi uncertainty [mag]',$
;     THICK = 0.5, XRANGE = [0,0.2], YRANGE = [0,0.2]
;device, /close 

; 5) Assign uniform PSF uncertainty and Write results to file
; ----------------------------------------------------------

; assign uniform PSF uncertainties
sx = sqrt(sx^2 + psferr_x^2)
sy = sqrt(sy^2 + psferr_y^2)
sf = sqrt(sf^2 + (psferr_f*fl)^2)

;sx = sqrt(sx^2 + dxl^2)
;sy = sqrt(sy^2 + dyl^2)
;sf = sqrt(sf^2 + dfl^2)


openw, out, photdir + 'astrophot.txt', /get_lun
;openw, out, photdir + 'astrophot_cal.txt', /get_lun
for i_star = 0, n_detected-1 do begin
  printf, out, xl[i_star], yl[i_star], fl[i_star], sxl[i_star], syl[i_star], sfl[i_star], cl[i_star], nl[i_star], FORMAT='(8F13.3)'
endfor
free_lun, out

END