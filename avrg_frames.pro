PRO avrg_frames

field = 9
chip = 	1

;~ for chip = chip, chip do begin

band = 'H'
indir = '/Users/amartinez/Desktop/PhD/HAWK/GNS_2/data/GNS/2021/'+band+'/' + strn(field) + '/ims/'
outdir = '/Users/amartinez/Desktop/PhD/HAWK/GNS_2/data/GNS/2021/'+band+'/' + strn(field) + '/cubes/'

pruebas= '/Users/amartinez/Desktop/PhD/HAWK/GNS_2/pruebas/'

outdir=pruebas



j=1

cube=readfits(indir +'chip1_cube' + strn(j) + '.fits')
sz = size(cube)
nax1 = sz[1]
nax2 = sz[2]
nax3 = sz[3]
wxy=1000
;~ new_cube=fltarr(nax1+2,nax2+2,nax3)
new_cube=fltarr(nax1,nax2,nax3)

n_frames=nax3-1 ;no the last frame

ref_ima=cube[nax1/2-wxy/2:nax2/2+wxy/2,nax1/2-wxy/2:nax2/2+wxy/2,1];use the second frame as referece cause the frist one look bad (stripes)
writefits, pruebas+'ref_ima.fits',ref_ima
;~ new_cube[1:nax1,1:nax1]=cube[*,*,1]
new_cube[*,*,1]=cube[*,*,1]


for i=0, n_frames-1 do begin 
    if i ne 1 then begin
		frame=cube[nax1/2-wxy/2:nax2/2+wxy/2,nax1/2-wxy/2:nax2/2+wxy/2,i]
	    ;~ writefits, pruebas+'frame'+strn(i)+'.fits',frame
	    
	correl_optimize, ref_ima, frame, x_off, y_off, MAGNIFICATION=4, /NUMPIX ; /NUMPIX is ESSENTIAL
    print, 'Offsets from correlation: ' + strn(x_off) + ', ' + strn(y_off)
    x_off=round(x_off)
	y_off=round(y_off)
	print,'X_off and Y_off rounded: ',x_off,y_off
	shifted=shift(cube[*,*,i],x_off,y_off)
	new_cube[*,*,i]=shifted	
	
	;~ new_cube[*,*,i]=cube[*,*,i]
	;~ new_cube[*,*,i]=shiftnw(new_cube[*,*,i],x_off,y_off)
	
	;~ new_cube[1+x_off:nax1+x_off,1+y_off:nax1+y_off,i]=cube[*,*,i]
	endif
endfor
writefits, pruebas + 'new_cube.fits', new_cube
END
