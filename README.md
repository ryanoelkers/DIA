# DIA
Delta Function Difference Imaging Code from Oelkers et al 2015, AJ, 149, 50

Hello! Thank you for showing interest in using the difference imaging code I have
written based off of Alard & Lupton 1998 and Miller 2008. It is a relatively 
simple code to use and feel free to modify it for the basis of your own codes. I only ask you please 
cite Oelkers et al., 2015, AJ, 149, 50 in addition to the papers above.

Before you can run the code, you will need 6 items. A science image, a reference image,
a text file with the science image name (img.txt), a text file with the reference image
name (ref.txt), a text file with all of the stars you would like to use for the kernel
(refstars.txt) and a parameter file (parms.txt). These files are for the stand alone 
version of the C code, if you use the PYTHON and IDL routines described below, these files
will be generated for you.

The parameter file is 1 line and has 4 items. First is the size of the fwhm box you would
like to use for the stamps (2xfwhm+1). The second is w for the kernel size of 2xw+1. The 
polynomial dimension of the kernel, 0 is constant, 1 is 1st order so on. And finally the 
number of stars you will use for the kernel. So as an example if you want 25x25pix boxes 
around each star for a 1st order 9x9 kernel with 20 reference stars the 1 line parms.txt 
file would look like: 12 4 1 20

You will need the cfitsio library and c math library installed as well as a c compiler, 
I typically use gcc. You will compile to code as normal i.e:

gcc oisdifference.c -L/usr/local/lib/ -I/usr/local/include/ -lcfitsio -lm

and to run you will simply run the .out file: ./a.out

The 5 IDL/Python programs make up a larder wrapper program to use call the C program for the image subtraction. Please feel free to tweak the routines how you see fit. If you find a bug please don't hesitate to contact me at ryan.j.oelkers@vanderbilt.edu.

The programs are written to work in this order:

1) clean.pro/.py - This program will apply bias subtraction, flat fielding and background subtraction to all frames. The program will also align the images to the first image in the directory. The user will need to update whether the images need bias subtraction or flat-fielding and the directories where the files are located. Also, if you want to sample the background finer than 128x128 pixels it can be updated. Usually I use 32x32 pixels. Additionally, the user will need to update the axis variable on line 22 to be the axis size of the image.

2) mk_master.pro/.py - This program will make the master frame by combining all available iamges. The user will need to update the directories and the axis size of the image. The routine uses a median combine.

3) refphot.pro/.py - This program will do basic aperture photometry on the reference frame. The user should upate the size of the aperture and the directories for the images. This will output a magnitude, flux and star position file.

4) bigdiff.pro/.py - This program is a large wrapper function for the oisdifference.c. The program selects the best isolated stars for the kernel and then runs the subtraction. Additionally, it gets the photometry for the stars on the differenced frame using aperture photometry and the positions from the master frame. These are output as flux files in the differenced image directory. The user will need to update directories, aperture and kernel information. 

5) mkraw.pro/.py - This program will take all of the flux files and combine them into single light curve files. The user should only need to update the directories.
