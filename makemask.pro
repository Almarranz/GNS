PRO MAKEMASK, common_path, bpm_name
  
; PURPOSE Create an array that contains 0 in all unsupported regions of
;         the image. Unsupported regions are, for example, the gaps
;         between the detectors or large clusters of bad pixels.
;         The unsupported regions cannot be repaired by interpolation,
;         as is the case for dead pixels.
;         Therefore, they must be masked.

field = '9'
band = 'H'

common_path = '/Users/amartinez/Desktop/PhD/HAWK/GNS_2/data/GNS/2021/' + band + '/' + field +'/ims/'

;~ common_path = '/data/GNS/2015/' + band + '/' + field +'/ims/'
bpm_name = 'bpm.fits'


debug = 0
bpm = readfits(common_path + bpm_name)
sz = size(bpm)
nax1 = sz[1]
nax2 = sz[2]

; Create a preliminary mask
; (this part is identical to the corresponding part
; in joinflats.pro
mask = fltarr(nax1,nax2)
mask[*,*] = 1
mask[0:10,*] = 0
mask[4073:4095,*] = 0
mask[2044:2051,*] = 0

;
; the edges of tmpmask are 1 pixel broader than
; the edges of the "real" mask
; so that it can be avoided to include the 
; detector edges and gaps into the search for dead pixel
; clusters.
;~ mask[*,*] = 1
;~ mask[0:10,*] = 0
;~ mask[4073:4095,*] = 0 ; why the edge is bigger on the right side of the mask?
;~ mask[2044:2051,*] = 0
;~ mask[*,4088:4095] = 0
;~ mask[*,2044:2051] = 0


tmpmask = fltarr(nax1,nax2)
tmpmask[*,*] = 1
tmpmask[0:11,*] = 0
tmpmask[4073:4095,*] = 0
tmpmask[2043:2052,*] = 0
mask[*,4087:4095] = 0
mask[*,2043:2052] = 0
check = where(tmpmask gt 0,counts)
;print, counts

for i = 0L, counts - 1 do begin
  xy = array_indices(bpm,check[i])
  xx = xy[0]
  yy = xy[1]
  if (bpm[xx,yy] gt 0) then begin
   region = search2d(bpm,xx,yy,0.99,1)
   if (n_elements(region) gt 5) then mask[region] = 0
  endif
endfor

; mask manually an isolated pixel in the middle of a bunch 
; of bad pixels (I found it when running cleancubes.pro)
;~ mask[392,677] = 0
writefits, common_path + 'mask.fits', mask

print, 'makemask.pro ended'

END
