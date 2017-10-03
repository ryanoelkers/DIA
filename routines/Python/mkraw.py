#this program will make the raw magnitude data files

#if you use this code, please cite Oelkers et al. 2015, AJ, 149, 50
#import the relevant libraries for basic tools
import numpy

#import relevant libraries for a list
import glob, os
from os import listdir
from os.path import isfile, join, exists

####UPDATE INFORMATION HERE####
###be sure to update the number of light curves on line 40, now for simplicity, it is set to 100
#useful directories
cdedir = '../IDL/' # code directory
difdir = '../diff/' # directory to put the differenced images
caldir = '../calib/' # directory with master frame information
lcdir  = '../lightcurve/' # directory to put the light curves

#get the flux lists
os.chdir(difdir) #changes to the raw image direcotory
files = [f for f in glob.glob("*.flux") if isfile(join(difdir, f))] #gets the relevant files with the proper extension
files.sort()
nfiles = len(files)
os.chdir(cdedir) #changes back to the code directory

#read in the master star list
ids, x, y = numpy.loadtxt(caldir+'starlist.txt', unpack = 1, delimiter = ',')

jd = numpy.zeros(nfiles)
mg = numpy.zeros((nfiles,len(ids)))
mge = numpy.zeros((nfiles,len(ids)))

#go through all the files
for ii in range (0, nfiles):
	id1, xx, yy, jds, mags, mages = numpy.loadtxt(difdir+files[ii], unpack = 1, delimiter = ',')
	jd[ii] = jds[0]
	mg[ii,:] = mags
	mge[ii,:] = mages

#openw new files
for ii in range(0, 100):
	if (ii < 10):
		nm = '0000'+str(ii)+'.lc'
	if (ii < 100) and (ii >= 10): 
		nm = '000'+str(ii)+'.lc'
	if (ii < 1000) and (ii >= 100):
		nm = '00'+str(ii)+'.lc'
	if (ii < 10000) and (ii >= 1000): 
		nm = '00'+str(ii)+'.lc'
	if (ii < 100000) and (ii >= 10000): 
		nm = '0'+str(ii)+'.lc'
	if (ii >= 100000):
		nm = str(ii)+'.lc'

	outdir = open(lcdir+nm, 'w')
	for jj in range(0, len(jd)):
		line = str(jd[jj])+','+str(mg[jj,ii])+','+str(mge[jj,ii])+'\n'
		outdir.write(line)
	outdir.close()
