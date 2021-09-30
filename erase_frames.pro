PRO erase_frames

field = 9
chip = 	1

;~ for chip = chip, chip do begin

band = 'H'
indir = '/Users/amartinez/Desktop/PhD/HAWK/GNS_2/data/GNS/2021/'+band+'/' + strn(field) + '/ims/'
cubos = '/Users/amartinez/Desktop/PhD/HAWK/GNS_2/data/GNS/2021/'+band+'/' + strn(field) + '/cubes_bad_frms/'

pruebas= '/Users/amartinez/Desktop/PhD/HAWK/GNS_2/pruebas/'

output='/Users/amartinez/Desktop/PhD/HAWK/GNS_2/data/GNS/2021/'+band+'/' + strn(field) + '/cubes/
;~ outdir=pruebas


ROBUST=0
;~ j=2

zero=fltarr(4096,4096)

cube=readfits(cubos +'cube34.fits',header)
cube=cube[*,*,[0,1,2,3,5]]
cube=[[[cube]],[[zero]]];added a zero dimesion for the average scripts works
;~ writefits,pruebas+'erased_cube34.fits',cube, header
writefits,output+'cube34.fits',cube, header


	
cube=readfits(cubos +'cube38.fits',header)
cube=cube[*,*,[2,3,4,6]]
cube=[[[cube]],[[zero]]];added a zero dimesion for the average scripts works
;~ writefits,pruebas+'erased_cube38.fits',cube, header
writefits,output+'cube38.fits',cube, header

cube=readfits(cubos +'cube42.fits',header)
cube=cube[*,*,1:6]
cube=[[[cube]],[[zero]]];added a zero dimesion for the average scripts works
;~ writefits,pruebas+'erased_cube42.fits',cube, header
writefits,output+'cube42.fits',cube, header

cube=readfits(cubos +'cube43.fits',header)
cube=cube[*,*,2:6]
cube=[[[cube]],[[zero]]];added a zero dimesion for the average scripts works
;~ writefits,pruebas+'erased_cube43.fits',cube, header
writefits,output+'cube43.fits',cube, header


f3=[10,22,23,40]
frms=n_elements(f3)
print,frms
for i = 0, frms-1 do begin
    print,f3[i]
	cube=readfits(cubos +'cube' + strn(f3[i]) + '.fits',header)
	sz = size(cube)
	nax1 = sz[1]
	nax2 = sz[2]
	nax3 = sz[3] 
	print,nax1,nax2,nax3
	
	cube=cube[*,*,3:nax3-2]
	cube=[[[cube]],[[zero]]];added a zero dimesion for the average scripts works

	;~ writefits,pruebas+'erased_cube'+strn(f3[i])+'.fits',cube, header
	writefits,output+'cube'+strn(f3[i])+'.fits',cube, header
endfor

f4=[7,13,44]
frms=n_elements(f4)
print,frms
for i = 0, frms-1 do begin
    print,f3[i]
	cube=readfits(cubos +'cube' + strn(f4[i]) + '.fits',header)
	sz = size(cube)
	nax1 = sz[1]
	nax2 = sz[2]
	nax3 = sz[3] 
	print,nax1,nax2,nax3
	
	cube=cube[*,*,4:nax3-2]
	cube=[[[cube]],[[zero]]];added a zero dimesion for the average scripts works

	;~ writefits,pruebas+'erased_cube'+strn(f4[i])+'.fits',cube, header
	writefits,output+'cube'+strn(f4[i])+'.fits',cube, header
endfor

f5=[17]
frms=n_elements(f5)
print,frms
for i = 0, frms-1 do begin
    print,f3[i]
	cube=readfits(cubos +'cube' + strn(f5[i]) + '.fits',header)
	sz = size(cube)
	nax1 = sz[1]
	nax2 = sz[2]
	nax3 = sz[3] 
	print,nax1,nax2,nax3
	
	cube=cube[*,*,5:nax3-2]
	cube=[[[cube]],[[zero]]];added a zero dimesion for the average scripts works

	;~ writefits,pruebas+'erased_cube'+strn(f5[i])+'.fits',cube, header
	writefits,output+'cube'+strn(f5[i])+'.fits',cube, header
endfor



































END
	
	
