#this program will make the raw magnitude data files

#if you use this code, please cite Oelkers & Stassun 2018
#import the relevant libraries for basic tools
import numpy

#import relevant libraries for a list
import glob, os
from os import listdir
from os.path import isfile, join, exists
from time import strftime

####UPDATE INFORMATION HERE####
###be sure to update the number of light curves on line 40, now for simplicity, it is set to 100
#useful directories
cdedir = '../code/phot/' # code directory
difdir = '../test/dif/' # directory to put the differenced images
caldir = '../code/master/fin/' # directory with master frame information
lcdir  = '../lc/' # directory to put the light curves

camera = '2'
ccd = '2'
####END UPDATE INFORMATION####

#get the flux lists
os.chdir(difdir) #changes to the raw image direcotory
files = [f for f in glob.glob("*.flux") if isfile(join(difdir, f))] #gets the relevant files with the proper extension
files.sort()
nfiles = len(files)
os.chdir(cdedir) #changes back to the code directory

#read in the master frame information to get the star count
ids, xm, ym, ttmg, mstmg, mstmge = numpy.loadtxt(caldir+camera+'_'+ccd+'_master.ap', unpack =1, usecols=(0,1,2,3,4,5))
nstars = len(ids)

jd = numpy.zeros(nfiles)
mg = numpy.zeros((nfiles, nstars))-99.00
er = numpy.zeros((nfiles, nstars))-99.00

print 'Splitting photometry files at '+strftime("%a, %d %b %Y %H:%M:%S")+'.'
#loop through each file to combine the files
for ii in range(0, 10):#nfiles):
	
	#read in the frame file
	idd, mx, my, jdh, mgh, mgeh = numpy.loadtxt(difdir+files[ii], unpack =1, usecols=(0,1,2,3,4,5),dtype = 'string')

	#put the data in the holder vectors
	jd[ii] = jdh[0]
	mg[ii, :] = mgh
	#guard against bad magnitude errors
	bd = numpy.where(mgeh == '********')
	if (len(bd) > 0):
		mgeh[bd]='-9.999999'
	er[ii, :] = mgeh


#loop through the stars make the light curves right now it is set to 100 to simiplify the code
output2=open(lcdir+'starinfo.txtxx', 'w')
for ii in range(0, 10):#nstars): 
	#make a file name from the ticid
	nm=str(long(ids[ii]))+'.lcxx'

	#write out the light curve
	output1=open(lcdir+nm, 'w')
	for jj in range(0, nfiles):
		 output1.write('%f %f %f\n' % (numpy.round(jd[jj], decimals = 6), numpy.round(mg[jj,ii], decimals = 6), numpy.round(er[jj,ii], decimals = 6)))
	output1.close()
	
	#write out the star information for easy access
	output2.write('%s %f %f %f %f\n' %(nm, numpy.round(xm[ii], decimals = 2), numpy.round(ym[ii], decimals = 2), numpy.round(mstmg[ii], decimals = 6), numpy.round(mstmge[ii], decimals = 6)))

output2.close()
print 'All done at '+strftime("%a, %d %b %Y %H:%M:%S")+'. See ya later alligator.'
