PRO DARK
field = '6'


in_path = '/Users/alvaromartinez/Desktop/Phd/HAWK/GNS_2/Dark/' + field + '/'


band = 'H'
list = 'list.txt'

NDIT = 7   ; number of sub-integrations in Dark
            ; Sometimes the cubes contain NDIT + 1 
            ; images because the last one is the mean
            ; of all sub-images.
            ; Setting the NDIR correctly assures
            ; discarding the mean.
            

out_path = '/Users/alvaromartinez/Desktop/PhD/HAWK/GNS_2/data/GNS/2021/H/' + field + '/ims/'


readcol, in_path+list, names, FORMAT='(A)'

n_dark = n_elements(names)
darkcube = []
for i = 0, n_dark-1 do begin
  cube = readfits(in_path + names[i])   
  darkcube = [[[]],cube[*,*,*]]
endfor

sz = size(darkcube)
n1 = sz[1]
n2 = sz[2]
dark = fltarr(n1,n2)
dark_sigma = fltarr(n1,n2)

for x = 0, n1-1 do begin
   for y = 0, n2-1 do begin
      data = darkcube[x,y,*]
      dark[x,y] = avg(data)
      dark_sigma[x,y] = stddev(data)/sqrt(3)
   endfor
endfor

writefits, out_path + 'dark.fits', dark
writefits, out_path + 'dark_sigma.fits', dark_sigma

dark_c1 = dark[0:2047,0:2047]
dark_c2 = dark[2048:4095,0:2047]
dark_c3 = dark[2048:4095,2048:4095]
dark_c4 = dark[0:2047,2048:4095]

writefits, out_path + 'dark_ext.fits', dark_c1,header,/app 
writefits, out_path + 'dark_ext.fits', dark_c2,header2,/app 
writefits, out_path + 'dark_ext.fits', dark_c3,header3,/app 
writefits, out_path + 'dark_ext.fits', dark_c4,header4,/app 

print, "dark.pro ended"

END
