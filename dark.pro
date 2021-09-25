PRO DARK
field = '9'

;~ in_path = '/home/data/raw/2015/Dark/2015-06-08/'
in_path = '/Users/amartinez/Desktop/PhD/HAWK/GNS_2/Dark/17-09-2021/'

band = 'H'
list = 'list.txt'

NDIT = 7   ; number of sub-integrations in Dark
            ; Sometimes the cubes contain NDIT + 1 
            ; images because the last one is the mean
            ; of all sub-images.
            ; Setting the NDIR correctly assures
            ; discarding the mean.
            

out_path = '/Users/amartinez/Desktop/PhD/HAWK/GNS_2/data/GNS/2021/H/' + field + '/ims/'
;~ out_path = '/data/GNS/2015/H/' + field + '/ims/'
;out_path = '/data/GNS/2015/'+ band +'/'+field + '/ims/'

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

print, "dark.pro ended"

END
