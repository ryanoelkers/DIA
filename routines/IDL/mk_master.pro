PRO mk_master

;this program will combine images to make a master frame

;if you use this code, please cite Oelkers et al. 2015, AJ, 149, 50

;;;;UPDATE INFORMATION HERE;;;;;

;useful directories
cdedir = '../IDL/' ;code directory
caldir = '../calib/' ; directory to output the master frame
clndir = '../clean/'; directory for the cleaned images

axs = 4096 ;; udpate here to change the image size
;;;END UPDATE INFORMATION;;;;


print, 'Beginning to make the master frame at ', systime(), '.'
        
;change into the appropriate directory
cd, clndir
spawn, 'ls *.fits', files
nfiles = n_elements(files)
cd, cdedir

;if no images then stop
if nfiles eq 0 then begin
	print, 'No files to combine.'
	stop
endif

;make a file holder
mref = fltarr(axs, axs, nfiles)
exps = fltarr(nfiles)
print, 'Combining images at '+systime()+'.'  

for jj = 0L, nfiles -1 do begin
    	img = readfits(clndir+files[jj], header,/silent) ; read in an image
	exps[jj] = sxpar(header,'EXPTIME')
    	mref[*,*,jj] = img ;place the image in an image array
endfor
        
print, 'Median-Combining the frames at '+systime()+'.'
        
;median the frame and then write the new image
medarr, mref, ref
expt = median(exps)

;write the image and add information to the header
writefits, caldir+'master.fits', ref
hold = readfits(caldir+'master.fits', head, /silent)
sxaddpar, head, 'EXPTIME', expt ; exposure time of the master
sxaddpar, head, 'NUM_COM', nfiles ; number of files that were combiend
writefits, caldir+'master.fits', hold, head

print, 'Finished with the master frame at ', systime(),'.'

print, 'All done! See ya later alligator.'

END
