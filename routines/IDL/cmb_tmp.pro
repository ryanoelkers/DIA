PRO cmb_tmp
;this program will combine the temporary images to make a master frame

;if you use this code, please cite Oelkers & Stassun 2018

;;;;;UPDATE HERE;;;;;;
;what field are you looking at?
camera = '2'
ccd = '2'

;useful directories
cdedir = '../code/master/' ;code directory
mstdir = '../code/master/frames/' ;directory where the temporary master frames reside
findir = '../code/master/fin/' ;directory for the final master frame
;;;;;;END UPDATE;;;;;;;;

;get the image list and the number of files which need reduction
cd, mstdir ;changes to the raw image direcotory
spawn, "ls "+camera+"_"+ccd+"_*.fits", files ;gets the relevant files with the proper extension
nfiles = n_elements(files)
cd, cdedir ;changes back to the code directory

;set up the holder for the final fiel count
img = readfits(mstdir+files[0], /silent)
nxy = size(img)
all_data = fltarr(nxy[1],nxy[2], nfiles)
expt = fltarr(nfiles)
num = fltarr(nfiles)

for ii =0l, nfiles-1 do begin

	;read in the image
	img_data = readfits(mstdir+files[ii], mhead, /silent)
	expt[ii] = sxpar(mhead, 'EXPTIME')
	num[ii] = sxpar(mhead, 'NUMCOMB')

	;add the image to the vector
	all_data[*,*,ii] = img_data 

	if (ii mod 10 eq 0) and (ii gt 0) then print, 'Finished with 10 images at '+systime()+'.'
endfor

;median combine the data
medarr, all_data, combined_data

; Write data to new file    
writefits, 'tmp_mst.fits', combined_data
mast = readfits('tmp_mst.fits', mhead, /silent)
sxaddpar, mhead, 'NUMCOMB', total(num)
sxaddpar, mhead, 'EXPTIME', median(expt)

;print the file
writefits, findir+camera+'_'+ccd+'_master_idl.fits', mast, mhead

print, "The master frame was created using a median of "+strcompress(total(num), /remove_all)+" images."

print, 'All done at '+systime()+'. See ya later alligator!'

END
