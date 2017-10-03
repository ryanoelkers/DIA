README for the DIA routines

The 5 IDL programs and 1 C program make up the image subtraction from Oelkers et al. 2015, AJ, 149, 50. The codes will run with the directories and images provided. Please feel free to tweak the routines how you see fit. If you find a bug please don't hesitate to contact me at ryan.j.oelkers@vanderbilt.edu.

The programs are written to work in this order:

1) clean.pro/.py - This program will apply bias subtraction, flat fielding and background subtraction to all frames. The program will also align the images to the first image in directory. The user will need to update whether the images need bias/flat field subtraction and the directories where the files are located. Also, if you want to sample the background finer than 128x128 pixels it can be updated. Usually I use 32x32 pixels. Additionally, the user will need to update the axs variable on line 22 to be the axis size of the image.

2) mk_master.pro/.py - This program will make the master frame by combining all available iamges. The user will need to update the directories and the axis size of the image. The routine uses a median combine.

3) refphot.pro/.py - This program will do basic aperture photometry on the reference frame. The user should upate the size of the apeture and the directories. This will output a magnitude, flux and star position file.

4) bigdiff.pro/.py - This program is a large wrapper function for the oisdifference.c. The program selects the best isolated stars for the kernel and then runs the subtraction. Additionally, it gets the photometry for the stars on the differenced frame using aperture photometry and the positions from the master frame. These are output as flux files in the differenced image directory. The user will need to update directories, aperture and kernel information. 

5) mkraw.pro/.py - This program will take all of the flux files and combine them into single light curve files. The user should only need to update the directories.
