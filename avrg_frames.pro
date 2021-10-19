PRO avrg_frames

field = 9


;~ for chip = chip, chip do begin

band = 'H'
indir = '/Users/amartinez/Desktop/PhD/HAWK/GNS_2/data/GNS/2021/'+band+'/' + strn(field) + '/ims/'
outdir = '/Users/amartinez/Desktop/PhD/HAWK/GNS_2/data/GNS/2021/'+band+'/' + strn(field) + '/cubes/'

pruebas= '/Users/amartinez/Desktop/PhD/HAWK/GNS_2/pruebas/'

;~ outdir=pruebas


ROBUST=0
;~ j=2
cubes_mod=[7,10,13,17,22,23,34,38,40,42,43,44,47]
;~ cubes_mod=[42,47]


mods=n_elements(cubes_mod)
print, 'Working on ', mods,' cubes'

for chip=1, 4 do begin
	for j=0, mods-1 do begin

		cube=readfits(indir +'chip'+strn(chip)+'_cube' + strn(cubes_mod[j]) + '.fits',header)
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
		;~ writefits, pruebas+'ref_ima.fits',ref_ima
		;~ new_cube[1:nax1,1:nax1]=cube[*,*,1]
		new_cube[*,*,1]=cube[*,*,1]


		for i=0, n_frames-1 do begin 
			if i ne 1 then begin ;this make sense if doing on cubes where the frist frame if of lower quality.
				frame=cube[nax1/2-wxy/2:nax2/2+wxy/2,nax1/2-wxy/2:nax2/2+wxy/2,i]
				;~ writefits, pruebas+'frame'+strn(i)+'.fits',frame
				
			correl_optimize, ref_ima, frame, x_off, y_off, MAGNIFICATION=1, /NUMPIX ; /NUMPIX is ESSENTIAL
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
		; initialize variables
		; ---------------------
		ssa_all = fltarr(nax1,nax2)
		ssa_sigma_all = fltarr(nax1,nax2)
		wt_all = fltarr(nax1,nax2)


		ssa = fltarr(nax1,nax2)
		wt = fltarr(nax1,nax2)
		ssa_sigma = fltarr(nax1,nax2)
		; ---------------------

		all_frames=new_cube[*,*,0:nax3-2]
		;~ writefits, pruebas+'all_frames.fit',all_frames
		lim=1
		if nax3 gt 3 then lim=2 else lim=1
		
		
		   for x = 0, nax1-1 do begin
			 for y = 0, nax2-1 do begin
			  vals = reform(all_frames[x,y,*])
			 ; require at least 3 measurements per pixel
			 ; for a reasonable estimation of sigma
			  good = where(reform(all_frames[x,y,*]) gt 0, complement = bad, ngood)
			  if ngood gt lim then begin 
			  ;~ if ngood gt 2 then begin ; have to changed this couse have cubes with only two frames.
			   vals = vals[good]
			   wt[x,y] = wt[x,y] + n_elements(good)
			  endif else begin
			   vals[*] = 0
			   wt[x,y] = 0
			  endelse
			  if ROBUST gt 0 then begin
				RESISTANT_Mean, vals, Sigma_CUT, vals_mean, vals_sigma, Num_RejECTED
			  endif else begin
				vals_mean = avg(vals)
				vals_sigma = stddev(vals)/sqrt(n_elements(vals)-1)
			  endelse
			  ssa[x,y] = vals_mean
			  ssa_sigma[x,y] = vals_sigma
			 endfor
			endfor
			new_cube[*,*,nax3-1]=ssa
			writefits, indir + 'chip'+strn(chip)+'_cube' + strn(cubes_mod[j]) + '.fits', new_cube,header
			;~ writefits, pruebas +'cube_average/'+ 'ssa_' + strn(cubes_mod[j]) + '.fits', ssa
			;~ writefits, pruebas + 'cube_average/'+'wt_' + strn(cubes_mod[j]) + '.fits', wt
			;~ writefits, pruebas + 'cube_average/'+'ssa_sigma' + strn(cubes_mod[j]) + '.fits', ssa_sigma
			;~ writefits, pruebas + 'cube_average/'+'new_cube_all_' + strn(cubes_mod[j]) + '.fits', new_cube, header
			print, 'Averaged cube: cube ',cubes_mod[j]
			
			;~ wt_all = wt_all + wt
			;~ ssa_all = ssa_all + wt*ssa
			;~ ssa_sigma_all = ssa_sigma_all + wt*ssa_sigma^2
			
			;~ good = where(wt_all gt 0, complement=bad)
			;~ ssa_all[good] = ssa_all[good]/wt_all[good]
			;~ if (bad[0] gt -1) then ssa_all[bad] = 0
			;~ ssa_sigma_all[good] = ssa_sigma_all[good]/wt_all[good]
			;~ if (bad[0] gt -1) then ssa_sigma_all[bad] = 0
			
			;~ new_cube[*,*,nax3-1]=ssa_all
			
			
			;~ writefits, pruebas + '_wt.fits', wt_all
			;~ writefits, pruebas + 'new_cube_all_wt.fits', ssa_all
			;~ writefits, pruebas + '_sigma.fits', sqrt(ssa_sigma_all)
	endfor
	print,'##################'
	print,'Done with chip ',chip
	print,'##################'
	
 endfor  


END
