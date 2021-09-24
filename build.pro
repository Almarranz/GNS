pro build, path, field, chip, i_x0, i_x1, i_y0, i_y1, x_off_tot, y_off_tot, xaxis, yaxis, rebfac, name, border, suffix

 
im =  fltarr(xaxis*rebfac,yaxis*rebfac)
nsubims =  fltarr(xaxis*rebfac,yaxis*rebfac)
exp =  fltarr(xaxis*rebfac,yaxis*rebfac)
sigma =  fltarr(xaxis*rebfac,yaxis*rebfac)

h = readfits(path + name + '0_0.fits.gz')
sz = size(h)
sub_size_x0 = sz[1]
sub_size_y0 = sz[2]
mask_borders = fltarr(sub_size_x0, sub_size_y0)
mask_borders[*,*] = 0
mask_borders[border:sub_size_x0-border-1, border:sub_size_y0-border-1] = 1     

for i_x = i_x0, i_x1 do begin
  
  for i_y = i_y0, i_y1 do begin


     tmp_subfield = fltarr(xaxis*rebfac,yaxis*rebfac)
     tmp_subfield_exp = fltarr(xaxis*rebfac,yaxis*rebfac)
     tmp_subfield_nsubims = fltarr(xaxis*rebfac,yaxis*rebfac)
     tmp_subfield_noise = fltarr(xaxis*rebfac,yaxis*rebfac)
     
     h = readfits(path + name + strn(i_x) + '_' +  strn(i_y) + suffix + '.fits.gz') * mask_borders
     h_noise = readfits(path + name + strn(i_x) + '_' +  strn(i_y) + '_sigma' + suffix  + '.fits.gz') * mask_borders
     h_exp = readfits(path + name + strn(i_x) + '_' +  strn(i_y) + '_wt.fits.gz')  * mask_borders ; no suffix: does not exist for jackknife sub-images

     tmp_subfield[0:sub_size_x0-1,0:sub_size_y0-1] = h            
     tmp_subfield_noise[0:sub_size_x0-1,0:sub_size_y0-1] = h_noise
     tmp_subfield_exp[0:sub_size_x0-1,0:sub_size_y0-1] = h_exp

     accept = where(tmp_subfield_exp gt 0, complement=reject)
     tmp_subfield_nsubims[accept] = 1
     tmp_subfield_nsubims[reject] = 0
      
     tmp_subfield = image_shift(tmp_subfield, x_off_tot[i_y-i_y0,i_x-i_x0], y_off_tot[i_y-i_y0,i_x-i_x0])  
     tmp_subfield_noise = image_shift(tmp_subfield_noise, x_off_tot[i_y-i_y0,i_x-i_x0], y_off_tot[i_y-i_y0,i_x-i_x0]) 
     tmp_subfield_exp = image_shift(tmp_subfield_exp, x_off_tot[i_y-i_y0,i_x-i_x0], y_off_tot[i_y-i_y0,i_x-i_x0]) 
     tmp_subfield_nsubims = image_shift(tmp_subfield_nsubims, x_off_tot[i_y-i_y0,i_x-i_x0], y_off_tot[i_y-i_y0,i_x-i_x0]) 
      
      ; Sub-pixel shifts can lead to strange weights
      ; because of interpolation near edges.
      ; Discard those pixels!
      accept = where(tmp_subfield_nsubims gt 0.99, complement=reject)
      tmp_subfield_nsubims[reject] = 0
      tmp_subfield[reject] = 0
      tmp_subfield_noise[reject] = 0
      tmp_subfield_exp[reject] = 0

      im =  im + tmp_subfield 
      sigma =  sigma + (tmp_subfield_noise)^2
      exp =  exp + tmp_subfield_exp
      nsubims =  nsubims + tmp_subfield_nsubims
      
      
      endfor
   endfor
      
      good = where(nsubims gt 0, complement=bad)
                  
      im[good] = im[good]/nsubims[good]      
      exp[good] = exp[good]/nsubims[good]
      sigma[good] = sigma[good]/nsubims[good]
      im[bad] = 0  
      sigma[bad] = 0
      exp[bad] = 0      
  
      writefits, path + field + '_' + chip + 'holo_' + strn(rebfac) + suffix + '.fits.gz', im, /COMPRESS         
      writefits, path + field + '_' + chip + 'noise_' + strn(rebfac) + suffix + '.fits.gz', sqrt(sigma), /COMPRESS 
      if (suffix eq '') then begin
        writefits, path + field + '_' + chip + 'exp_' + strn(rebfac) +'.fits.gz', exp, /COMPRESS 
        writefits, path + field + '_' + chip + 'nsubims_' + strn(rebfac) +'.fits.gz', nsubims, /COMPRESS 
      endif

   END

