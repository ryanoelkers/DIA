#this program will combine images to make a master frame

#if you use this code, please cite Oelkers et al. 2015, AJ, 149, 50

#import the relevant libraries for basic tools
import pyfits
import numpy as np
import scipy
from scipy import stats
from os import path
import math
import time

#import relevant libraries for a list
import glob, os
from os import listdir
from os.path import isfile, join, exists

###UPDATE HERE#####
#useful directories
cdedir = '../Python/' #code directory
caldir = '../calib/' # directory to output the master frame
clndir = '../clean/'# directory for the cleaned images
###END UPDATE###

#get the image list and the number of files which need reduction
os.chdir(clndir) #changes to the raw image direcotory
files = [f for f in glob.glob("*.fits") if isfile(join(clndir, f))] #gets the relevant files with the proper extension
files.sort()
nfiles = len(files)
os.chdir(cdedir) #changes back to the code directory

#iterate through the files
expt = [None]*nfiles
cnt = 0 #counter for the number of images used
for ii in range(0,len(files)):
        if (ii == 0):  # get size on first iteration only
		nx = pyfits.getval(clndir+files[0], 'NAXIS2')
                ny = pyfits.getval(clndir+files[0], 'NAXIS1')
                all_data = np.ndarray(shape=(len(files),nx,ny))

	#read in the image
	img_data = pyfits.getdata(clndir+files[ii])
	expt[ii] = pyfits.getval(clndir+files[ii],'EXPTIME')
	all_data[ii] = img_data 
	cnt = cnt+1
	#median combine the data
        combined_data = np.median(all_data,axis=0)

# Write data to new file    
new_image = pyfits.PrimaryHDU(combined_data)
new_image.header.set('NUMCOMB', cnt)
new_image.header.set('EXPTIME', np.median(expt))
new_image.writeto(caldir+'master.fits',clobber=True)
print "The master frame was created using a median of "+str(cnt)+" images."

del all_data, img_data # clear up some memory

