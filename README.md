# GNS 2
### on branch: unwid


GNS pipeline adapted for unwindowed frames. 

1. dark.pro
2. Create flat
3. sky.pro
4. fullbpm.pro
5. makemask.pro
___
### 5.1 mask_poligon.pro
There are some bright pixels on chip4 that doesnt show on the bpm, neither the mask. This script masks them
___
6. reduce.pro
___
### 6.1 erase_frames.pro
Erases bad frames and therefore transforms the cubes thrid dimension. Makes the last frame blank.
___
7. cleancubes.pro
8. makecubes.pro (little change). Apadted for cubes with different number of frames.
___
### 8.1 avrg_frames.pro 
Last frame of cubes which had bad frames (erased) is not the average anymore. This scripts makes the average of first frames and put it in the last slice of the cube. It ONLY does this on cubes with bad frames (the ones used by erase_frames.pro) and leaves the rest untouched.
___
9. dejitter.pro (some changes). The x_off y_off WASN'T applied to the masks. I changed it.
___
10. longexposure.pro
11. findstars.pro
12. This step consists in three scripts. Only to be run once per field.
* prep_refim.pro
* extarctpsf.pro
* deepastro.pro
___
13. alignquadrants new verision. 
Instead of clinking on common stars to figure out the x ,y displacement and manually figure the rotation angle between lists, runs astroalign packge on GNS list to get the transformation matrix, transforms the list, uses it to find common stars and  then aligns with poliwarp with VVV
*  ### 13. 1 astroaling_VVV.py
* ### 13.2  aa_alignquadrants.pro
___
14. mosaic.pro
15. aligframes.pro.(some chnges). I changed the size of the cropping areas to extrac aligned frames of the big  cambas.
16. subcubes.pro. The cube are 2700x2700 and the subcube are 900x900 (insted of 600x600). 


