# GNS
############################

on branch: unwid

############################
scripts GNS
1. dark.pro
2. Create flat
3. sky.pro
4. fullbpm.pro
5. makemask.pro
6. reduce.pro
7. cleancubes.pro.
----------------------------------
Here I added an extra script: erase_frames.pro
Erases bad frames and there fore transfor the cubes thrid dimension.
----------------------------------
8. makecubes.pro (little change).
----------------------------------
Here I added an extra script: avrg_frames.pro 
Last frme of cubes which had bad frames (erased) is not the average anymore. This scripts makes the average of last frames. 
Should I make average for only those cubes from  which I erased some frames and left the original ones untouched??
----------------------------------
9. dejitter.pro (some changes)
---
in dejitter.pro. the x_off y_off WASN'T applied to the masks. I changed it.
---
10. longexposure.pro


