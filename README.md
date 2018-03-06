# DIA
Delta Function Difference Imaging Code from Oelkers & Stassun 2018, submitted to ApJ Letters

Hello! Thank you for showing interest in using the difference imaging code I have written based off of Alard & Lupton 1998 and Miller et al. 2008. I hope that it is relatively simple to use, and I encourage you to modify any of the scripts to form the basis of your own codes. I only ask you please cite Alard & Lupton 1998, Alard 2000 and Miller+2008, Oelkers+2015, and Oelkers & Stassun 2018,  given their large contributions to his work.

The scripts in this repository are written in both PYTHON and IDL and are nearly identical in their outputs. The data products from Oelkers & Stassun 2018 used the IDL scripts (mainly because I can code and de-bug IDL faster), but I have tested all scripts on multiple machines and they should run without issue. If you do find a bug, or would like to suggest a more efficent way to code any script please don't hesitate to email me at ryan.j.oelkers@vanderbilt.edu. The code is currently (3/6/2018) in version 0 status and will be updated prior to and after the first TESS data release.

Alternatively, you can run the differencing code as a stand along C code, but you will need to manually generate 6 files. A science image, a reference image, a text file with the science image name (img.txt), a text file with the reference image name (ref.txt), a text file with all of the stars you would like to use for the kernel (refstars.txt), and a parameter file (parms.txt). If you use the PYTHON and IDL routines described below, these files will be generated for you.

The parameter file is 1 line and has 4 items. First is the size of the fwhm box you would like to use for the stamps (2xfwhm+1). The second is w for the kernel size of 2xw+1. The polynomial dimension of the kernel, 0 is constant, 1 is 1st order so on. And finally the number of stars you will use for the kernel. So as an example if you want 25x25pix boxes around each star for a 1st order 9x9 kernel with 20 reference stars the 1 line parms.txt file would look like: 12 4 1 20

For the C code to run properly, you will need the cfitsio library and c math library installed as well as a c compiler, I typically use gcc. You will compile to code as normal i.e:

gcc oisdifference.c -L/usr/local/lib/ -I/usr/local/include/ -lcfitsio -lm

and to run you will simply run the .out file: ./a.out

The IDL/Python scripts make up a larger pipeline which can be used to clean and align images, generate the master frame, call the C program for the image subtraction, generate light curves, and do a basic detrending to the light curves based on magnitude. 

The programs are written to work in this order:

clean.pro/.py - This program will optionally apply a bias subtraction, flat fielding, background subtraction and align images to the first image in the list, provide a WCS solution exists. The user will need to update whether the images need bias subtraction or flat-fielding and the directories where the files are located. Also, if you want to sample the background finer/courser than 32x32 pixels it can be updated. Usually I use 32x32 pixels. Additionally, the user will need to update the axis variable to be the axis size of the image. It should be noted, the PYTHON routine actually splits the image into 512x512 boxes for the sky subtraction to save memory. This can be changed but beware it may slow down your computer.

mk_master.pro/.py & cmb_tmp.pro/.py- This program will make the master frame by combining all available images. The scripts are run in two steps: mk_master will make temporary master frames which combine 50 images at a time. Then cmb_tmp, will combine the temporary frames into a single master frame. The user will need to update the directories. 

refphot.pro/.py - This program will do basic aperture photometry on the reference frame and should tell you the appropriate optimal aperture for the images. The user should update  the directories for the images. This will output a magnitude, flux and star position file. The code assumes you have a ra/dec position list for the stars, but you can also tweak the script to provide stars via the 'find' routine.

bigdiff.pro/.py - This program is a large wrapper function for oisdifference.c. The program selects the best isolated stars for the kernel and then runs the subtraction. Additionally, it gets the photometry for the stars on the differenced frame using aperture photometry and the positions from the master frame. These are output as flux files in the differenced image directory. The user will need to update directories, aperture and kernel information. Typically, I use a 5x5 pixel kernel with no polynomial function. Usually, the slight increase in precision is not worth the large increase in runtime, but its up to you.

mkraw.pro/.py - This program will take all of the flux files and combine them into single light curve files. The user should only need to update the directories.

detrend.pro/.py - This program will detrend the light curves using the entire data set. It looks for the nearest 1000 stars of similar magntidue, generates a trend based on similar offsets with time and subtracts them from the star. It does not work for every star and may over fit some stars but I have found that usually >95% of the stars have improved rms values. 

To properly run the IDL and Python codes, you will need some additional libraries and routines, they include:
IDL:
IDL Astronomy Library (https://idlastro.gsfc.nasa.gov/)

PYTHON:
numpy (http://www.numpy.org/)
scipy (https://www.scipy.org/)
astropy (http://www.astropy.org/)
FITS_tools (https://pypi.python.org/pypi/FITS_tools)
photutils (https://photutils.readthedocs.io/en/stable/)

C:
cfitsio (https://heasarc.gsfc.nasa.gov/fitsio/fitsio.html)

Please let me know if I have missed anything, you still have questions, or you think you have found a bug.

Thanks!
--Ryan
