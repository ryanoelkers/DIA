PRO mk_master
;this program will combine images to make the temprorary master frames

;if you use this code, please cite Oelkers & Stassun 2018

;;;;;;;;;;UPDATE HERE;;;;;;;;;
;what field are you looking at?
camera = '2'
ccd = '2'
blknum = 50 ;how many images go into each holder?

;useful directories
cdedir = '.../code/master/' ;code directory
clndir = '../clean/' ;directory where the cleaned images reside

;;;;;;;END UPDATE;;;;;;;

;get the image list and the number of files which need reduction
readcol, cdedir+'images.dat', files, format = '(a)', /silent, count = nfiles

;iterate through the files
cnt = 0 ;counter for the number of images used
kk = 0 ;the jumper for file placement

for ii =0l, nfiles-1 do begin

	if (ii eq 0) then begin  ; get size on first iteration only
		img = readfits(clndir+files[0], head, /silent)
		nxy = size(img)
	        all_data = fltarr(nxy[1],nxy[2],blknum)
		expt = fltarr(blknum)
	endif

	;read in the image
	img_data = readfits(clndir+files[ii], head, /silent)
	expt[cnt] = sxpar(head, 'EXPOSURE')

	;add the image to the vector
	all_data[*,*,cnt] = img_data 
	cnt = cnt+1

	if (ii mod 10 eq 0) and (ii gt 0) then print, 'Finished with 10 images at '+systime()+'.'

	if (ii eq nfiles-1) or ((ii+1) mod blknum eq 0) then begin
		
		;median combine the data
		medarr, all_data, combined_data

		; Write data to new file    
		writefits, 'tmp_mst.fits', combined_data
		mast = readfits('tmp_mst.fits', mhead, /silent)
		sxaddpar, mhead, 'NUMCOMB', cnt
		sxaddpar, mhead, 'EXPTIME', median(expt)

		;print the file with the appropriate counter
		if (kk lt 10) then writefits, cdedir+'frames/'+camera+'_'+ccd+'_master_0'+strcompress(kk, /remove_all)+'_idl.fits', mast, mhead
		if (kk ge 10) and (kk lt 100) then writefits, cdedir+'frames/'+camera+'_'+ccd+'_master_'+strcompress(kk, /remove_all)+'_idl.fits', mast, mhead

		print, "The master frame hold was created using a median of "+strcompress(cnt, /remove_all)+" images."
		kk = kk+1
		cnt = 0

	endif
endfor
print, 'All done at '+strcompress(systime(), /remove_all)+'. See ya later alligator!'
END
