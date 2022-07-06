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


openr, inp,(common_path + 'cube_info_'+ field +'.txt'), /get_lun

while not EOF(inp) do begin
readf, inp, i0,i1, i2, i3, i4, i5, i6, i7,i8, FORMAT = '(9I)'
indices = [i1,i2,i3,i4,i5,i6,i7,i8]
cube = round(i0)
ind = indices[where(indices gt 0)] - 1

print,cube,ind
stop
help, ind
for chip = 1, 4 do begin
      cube = readfits(indir + 'chip'+ strn(chip) + '_cube' + strn(round(i0)) + '.fits.gz', header)
    help, cube[*,*,ind]
      writefits, outdir + 'chip'+ strn(chip) + '_cube' + strn(round(i0)) + 'test.fits.gz', cube[*,*,ind], header, /COMPRESS
      print, 'Finished offset ' +strn(round(i0)) + ', chip ' + strn(chip)
  endfor

endwhile

END