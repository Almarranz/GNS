PRO avrg_frames

field = 6


;~ for chip = chip, chip do begin

band = 'H'
indir = '/Users/alvaromartinez/Desktop/PhD/HAWK/GNS_2/data/GNS/2021/'+band+'/' + strn(field) + '/ims/'
pruebas= '/Users/alvaromartinez/Desktop/PhD/HAWK/GNS_2/pruebas/'

; indir = pruebas





nam=''
openr, inp, (pruebas + 'bad_cubes_'+strn(field)+'.txt'), /get_lun  ; open input file for reading
while (not (EOF(inp))) do begin
   readf, inp,nam
   print, nam


for chip=1, 4 do begin
		cube=readfits(indir +'chip'+strn(chip)+'_cube' + strn(nam) + '.fits.gz',header)
		;~ cube=readfits(indir +'chip1_cube' + strn(j) + '.fits',header)
		sz = size(cube)
		nax1 = sz[1]
		nax2 = sz[2]
		nax3 = sz[3]
		wxy=1000
		;~ new_cube=fltarr(nax1+2,nax2+2,nax3)
		new_cube=fltarr(nax1,nax2,nax3)

		n_frames=nax3-1 ;no the last frame

		ref_ima=cube[nax1/2-wxy/2:nax2/2+wxy/2,nax1/2-wxy/2:nax2/2+wxy/2,1];use the second frame as referece cause the frist one look bad (stripes)
		new_cube[*,*,1]=cube[*,*,1]


		for i=0, n_frames-1 do begin 
			if i ne 1 then begin ;this make sense if doing on cubes where the frist frame if of lower quality.
				frame=cube[nax1/2-wxy/2:nax2/2+wxy/2,nax1/2-wxy/2:nax2/2+wxy/2,i]
				
				
			correl_optimize, ref_ima, frame, x_off, y_off, MAGNIFICATION=1, /NUMPIX ; /NUMPIX is ESSENTIAL
			print, 'Offsets from correlation: ' + strn(x_off) + ', ' + strn(y_off)
			x_off=round(x_off)
			y_off=round(y_off)
			print,'X_off and Y_off rounded: ',x_off,y_off
			shifted=shift(cube[*,*,i],x_off,y_off)
			new_cube[*,*,i]=shifted	
			
			endif
		endfor
		
		all_frames=new_cube[*,*,0:nax3-2]
		avereged_frame = avg(all_frames,2)
		new_cube[*,*,-1] = avereged_frame
		
		writefits, indir +'chip' + strn(chip) + '_cube' + strn(i) + '.fits.gz', new_cube, header
; 	print,'##################'
; 	print,'Done with chip ',chip
; 	print,'##################'
	
 endfor  
     print,'********************'
 	print,'Done with cube ',nam
 	print,'********************'
endwhile

END
