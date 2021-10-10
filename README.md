# GNS 2
### on branch: unwid


scripts GNS
1. dark.pro
2. Create flat
3. sky.pro
4. fullbpm.pro
5. makemask.pro
6. reduce.pro
7. cleancubes.pro.
___
Here I added an extra script: erase_frames.pro
Erases bad frames and there fore transfor the cubes thrid dimension.
___
8. makecubes.pro (little change).
___
Here I added an extra script: avrg_frames.pro 
Last frme of cubes which had bad frames (erased) is not the average anymore. This scripts makes the average of last frames. 
Should I make average for only those cubes from  which I erased some frames and left the original ones untouched??
___
9. dejitter.pro (some changes)
in dejitter.pro. the x_off y_off WASN'T applied to the masks. I changed it.
___
10. longexposure.pro
11. findstars.pro
12. This is step consists in three scripts. Only to be run once.
* prep_refim.pro
* extarctpsf.pro
* deepastro.pro
___
13. alignquadrants.pro new verision. Instead of clinking on stars to figure out the x ,y displacement between list, run astroalign packge on GNS lits to get the transformation matrix, transfor the list and then align with poliwarp with VVV
* astroaling_VVV_for_testing.py
* aa_alignquadrants.pro
___


