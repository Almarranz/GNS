PRO SELECTGOOD_AM;, chip

; We are going to elimante the bad frames before dejitter.pro, in order to avoind
; uses bad frames with  correl optimace

field = '6'
band = 'H'

common_path = '/Users/alvaromartinez/Desktop/PhD/HAWK/GNS_2/data/GNS/2021/' + band + '/' + field +'/ims/'
pruebas = '/Users/alvaromartinez/Desktop/PhD/HAWK/GNS_2/pruebas/'
cubes = '/Users/alvaromartinez/Desktop/PhD/HAWK/GNS_2/data/GNS/2021/' + band + '/' + field +'/cubes/'
morralla = '/Users/alvaromartinez/Desktop/morralla/'
indir = common_path
outdir = common_path

; indir = pruebas
; outdir = pruebas

n_offset = 70
openr, inp,common_path + 'cube_info.txt', /get_lun
print, inp

; Change format of image cubes
; -----------------------------

for j = 1, n_offset do begin
   readf, inp, i1, i2, i3, i4, i5, i6, i7,i8, FORMAT='(8I)'
   indices = [i1,i2,i3,i4,i5,i6,i7,i8]
   
    ind = indices[where(indices gt 0)] - 1
;   ind = [ind,7]
   
  for chip = 1, 4 do begin
      cube = readfits(indir + 'chip'+ strn(chip) + '_cube' + strn(j) + '.fits.gz', header)
    help, cube[*,*,ind]
      writefits, outdir + 'chip'+ strn(chip) + '_cube' + strn(j) + '.fits.gz', cube[*,*,ind], header, /COMPRESS
      print, 'Finished offset ' + strn(j) + ', chip ' + strn(chip)
  endfor
endfor

free_lun, inp

END
