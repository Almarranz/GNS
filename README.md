# GNS 2
### on branch: unwid
#### This scripts are compared to the ones in teatime at June 14th 2022

GNS pipeline adapted for unwindowed frames. 

0. file_sorter.py. This separte the darks and th flats to be use with esorex recipes. Also generates the .sof files needed for esorex
> Note: at the moment there is no **esorex** not **gasgano** installed at teatime. SO I doing this locally for now

> Note2: you need to make a folder with date of the flat in the Flat directory
1. dark.pro. we do the dark with IDL. Genretaes two different darks,one on the canvas one with extentions
1.1 Esorex recipes

> esorex hawki_twilight_flat_combine flat.sof. Run ot in the flat folder, in the date folder

>>Note: the flat hawki_twilight_flat_combine needs the dark with each chip in different extension
2. joint_flats.pro. generates the flat, bpm_H and the mask.
3. sky.pro
4. fullbpm.pro
5. makemask.pro
6. reduce.pro
7. cleancubes.pro
8. makecubes.pro 

### 8.1 images_checking.py 
Creates 4 fits folders with a cube containning all the frames for each chip, for
easier visual inspectiom. It also generates text files matches each slice in the cube
with the corresponding pointing. You have to pass the information of the cubes with bad frames to **cube_info_file.py**

### 8.2 cube_info_file.py
Generates de **cube_info_<field>.txt** file that selectgood will use

### 9. selectgood_AM.pro. 
Eliminates bad frames
> NOTE: I have changed the order and placed selectgood.pro before dejitter.pro. This make sense only in the case where the first cube contains bad frames that can affect to the long exposure frame, which is the one correl_optimaze will use for the fine alignment

### 9.1 avrg_frames.pro. 
Create an average exposure (last frame) for those cubes that have been trimmed.

10. dejitter.pro.
11. longexposure.pro
12. findstars.pro (Rainer made little change)
13. This step consists in three scripts. Only to be run once per field.
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



