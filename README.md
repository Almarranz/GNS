# GNS 2
### on branch: unwid


GNS pipeline adapted for unwindowed frames. 

0. file_sorter.py. This separte the darks and th flats to be use with esorex recipes. Also generates the .sof files needed for esorex
> Note: at the moment there is no **esorex** not **gasgano** installed at teatime. SO I doing this locally for now
1. Esorex recipes
> esorex hawki_dark_combine dark.sof.

> esorex hawki_twilight_flat_combine flat.sof.

>>Note: the dark.fits file is no longe produce by the dark.pro file, but for esorex dark recipe Instead

2. joint_flats.pro. generates the flat, bpm_H and the mask.
3. sky.pro
4. fullbpm.pro
5. makemask.pro
___
### 5.1 mask_poligon.pro
There are some bright pixels on chip4 that doesnt show on the bpm, neither the mask. This script masks them
>(Rainer doesnt mask these `pixels, so didnt use this script, aperently)
___
6. reduce.pro
___
### 6.1 erase_frames.pro
Erases bad frames and therefore transforms the cubes thrid dimension. Makes the last frame blank.
>(it seems that Raner has replaced this with SELECTGOOD.pro)
___
7. cleancubes.pro
8. makecubes.pro (little change). Apadted for cubes with different number of frames.
___
### 8.1 avrg_frames.pro 
Last frame of cubes which had bad frames (erased) is not the average anymore. This scripts makes the average of first frames and put it in the last slice of the cube. It ONLY does this on cubes with bad frames (the ones used by erase_frames.pro) and leaves the rest untouched.
>(Rainer doesnt use this one I think.)
___
9. dejitter.pro (some changes). The x_off y_off WASN'T applied to the masks. I changed it.
>(it seems Rainer has inroduced a little change here)
___
10. longexposure.pro
11. findstars.pro (Rainer made little change)
12. This step consists in three scripts. Only to be run once per field.
* prep_refim.pro
* extarctpsf.pro
* deepastro.pro
___
13. alignquadrants new verision. 
Instead of clinking on common stars to figure out the x ,y displacement and manually figure the rotation angle between lists, runs astroalign packge on GNS and VVV lists to get the transformation matrix, transforms the list, uses it to find common stars and  then aligns with poliwarp with VVV
*  ### 13. 1 astroaling_VVV.py
>(Rainer has divided this one into 4, one for each chip)
* ### 13.2  aa_alignquadrants.pro
>(Rainer has simplyfied this one. I think he uses the old one)
___
14. mosaic.pro
15. aligframes.pro.(some chnges). I changed the size of the cropping areas to extrac aligned frames of the big  cambas.
(Rainer have changed this one considerably)

16. subcubes.pro. The cube are 2700x2700 and the subcube are 900x900 (insted of 600x600). 



