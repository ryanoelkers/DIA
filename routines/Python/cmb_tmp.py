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
#what field are you looking at?
camera = '2'
ccd = '2'

#useful directories
cdedir = '../code/master/' #code directory
mstdir = '../code/master/frames/' #directory where the cleaned images reside
findir = '../code/master/fin/' #directory for the final master frame
###END UPDATE###

#get the image list and the number of files which need reduction
os.chdir(mstdir) #changes to the raw image direcotory
files = [f for f in glob.glob(camera+"_"+ccd+"_*.fits") if isfile(join(mstdir, f))] #gets the relevant files with the proper extension
files.sort()
nfiles = len(files)
os.chdir(cdedir) #changes back to the code directory

#set up the holder for the final fiel count
nx = pyfits.getval(mstdir+files[0], 'NAXIS2')
ny = pyfits.getval(mstdir+files[0], 'NAXIS1')
all_data = np.ndarray(shape=(nfiles,nx,ny))
expt = np.zeros(nfiles)
num = np.zeros(nfiles)

for ii in range(0,nfiles):

	#read in the image
	img_data = pyfits.getdata(mstdir+files[ii])
	expt[ii] = pyfits.getval(mstdir+files[ii],'EXPTIME')
	num[ii] = pyfits.getval(mstdir+files[ii],'NUMCOMB')

	#add the image to the vector
	all_data[ii] = img_data 

	if (ii % 10 == 0) and (ii > 0):
		print 'Finished with 10 images at '+str(time.strftime("%a %d %b %Y %H:%M:%S"))+'.'

#median combine the data
combined_data = np.median(all_data,axis=0)

# Write data to new file    
new_image = pyfits.PrimaryHDU(combined_data)
new_image.header.set('NUMCOMB', np.sum(num))
new_image.header.set('EXPTIME', np.median(expt))

#print the file with the appropriate counter
new_image.writeto(findir+camera+'_'+ccd+'_master_py.fits',clobber=True)

print "The master frame was created using a median of "+str(np.sum(num))+" images."

print "All done. See ya later alligator!"

del all_data, img_data # clear up some memory

