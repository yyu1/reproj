;Main program for reprojecting an image.
;The fuzzy version takes a weighted mean of surrounding pixels from input image

;Input file setup - normally, we shouldn't read everything into memory
;But the files here should be under 10 gig in size and the 16GB memory can handle it
;This makes reading individual pixel values much faster

;in_xdim = 43200ULL;  Africa
;in_ydim = 38400ULL;  Africa
;in_file = '/Volumes/Global_250m/output/afr/v3/maxent_afr_hlorey_v3_filled.flt' ;Africa
;in_tl_x = -3335851.5588D     ;top left corner of top left pixel of input image  (in MODIS grid m)
;in_tl_y = 4447802.0784D     ;top left corner of top left pixel of input image  (in MODIS grid m)
;in_x_size = 231.65635825D ;pixel width in meters
;in_y_size = 231.65635825D ;pixel height in meters
;
;;Target file setup
;out_xdim = 28500ULL;  Africa
;out_ydim = 30000ULL;  Africa
;out_file = '/Volumes/Global_250m/output/afr/v3/maxent_afr_hlorey_v3_filled_9.6sec.flt' ;Africa
;out_tl_x = -20.0D    ;top left corner of top left pixel of input image  (in degrees longitude)
;out_tl_y = 40.0D     ;top left corner of top left pixel of input image  (in degrees latitude)
;out_x_size = 0.0026666666666D ;pixel width degrees longitude
;out_y_size = 0.0026666666666D ;pixel height degrees latitude
;background_val = 0


;South America
in_xdim = 38400ULL  
in_ydim = 43200ULL  
in_file = '/Volumes/Global_250m/output/sam/v3/maxent_sam_hlorey_using_newprior.flt' 
in_tl_x = -11119505.1960     ;top left corner of top left pixel of input image  (in MODIS grid m)
in_tl_y = 3335851.5590     ;top left corner of top left pixel of input image  (in MODIS grid m)
in_x_size = 231.65635825D ;pixel width in meters
in_y_size = 231.65635825D ;pixel height in meters
out_xdim = 34125ULL  
out_ydim = 34125ULL 
out_file = '/Volumes/Global_250m/output/sam/v3/maxent_sam_hlorey_using_newprior_9.6sec_fuzzy.flt'
out_tl_x = -120.0D    ;top left corner of top left pixel of input image  (in degrees longitude)
out_tl_y = 30.0D     ;top left corner of top left pixel of input image  (in degrees latitude)
out_x_size = 0.0026666666666D ;pixel width degrees longitude
out_y_size = 0.0026666666666D ;pixel height degrees latitude
background_val = 0

;South America slice
;in_xdim = 38400ULL  
;in_ydim = 43200ULL  
;in_file = '/Volumes/Global_250m/output/sam/v3/maxent_sam_hlorey_using_newprior.flt' 
;in_tl_x = -11119505.1960     ;top left corner of top left pixel of input image  (in MODIS grid m)
;in_tl_y = 3335851.5590     ;top left corner of top left pixel of input image  (in MODIS grid m)
;in_x_size = 231.65635825D ;pixel width in meters
;in_y_size = 231.65635825D ;pixel height in meters
;out_xdim = 13000ULL  
;out_ydim = 4000ULL 
;out_file = '/Volumes/Global_250m/output/sam/v3/maxent_sam_hlorey_using_newprior_9.6sec_slice_fuzzy2.flt'
;out_tl_x = -82.32084179   ;top left corner of top left pixel of input image  (in degrees longitude)
;out_tl_y = 1.91818712   ;top left corner of top left pixel of input image  (in degrees latitude)
;out_x_size = 0.0026666666666D ;pixel width degrees longitude
;out_y_size = 0.0026666666666D ;pixel height degrees latitude
;background_val = 0


;allocate memory
print, 'Allocating memory...'
in_image = fltarr(in_xdim,in_ydim)
print, 'reading input file...'
openr, 1, in_file
readu, 1, in_image
close, 1

;Go through target file one line at a time, and grab nearst neighbor pixel from source (use center pixel location)
print, 'converting image...'
openw, 1, out_file

out_line = fltarr(out_xdim)
out_xcoords = fltarr(out_xdim)
out_ycoords = fltarr(out_xdim)
in_ind_x = ulon64arr(out_xdim)
in_ind_y = ulon64arr(out_xdim)
xcoeff = fltarr(out_xdim)
for i=0ULL, out_xdim-1 do begin
	xcoeff[i] = i+0.5
endfor

weights = [[0.05,0.05,0.05],[0.05,0.6,0.05],[0.05,0.05,0.05]]

for j=0ULL, out_ydim-1 do begin
	out_line[*] = 0
	if (j mod 10000 eq 0) then print, 'line ', j

	;set up values for middle of pixel coordinates of current output line
	out_xcoords[*] = out_tl_x + out_x_size*xcoeff
	out_ycoords[*] = out_tl_y - (j+0.5)*out_y_size

	ddlonlat2modxy, out_xcoords, out_ycoords, in_xcoords, in_ycoords

	;convert in_xcoords and  in_ycoords to image indices
	in_ind_x[*] = floor((in_xcoords-in_tl_x)/in_x_size, /l64)
	in_ind_y[*] = floor((in_tl_y-in_ycoords)/in_y_size, /l64)

	;Look for out of bounds pixels
	x_obb_index = where((in_ind_x lt 1) or (in_ind_x ge in_xdim-1), xcount)
	y_obb_index = where((in_ind_y lt 1) or (in_ind_y ge in_ydim-1), ycount)

	;temporarily setting out of bounds indices to 0,0
	if (xcount gt 0) then begin
		in_ind_x[x_obb_index]=1
	endif
	if (ycount gt 0) then begin
		in_ind_y[y_obb_index]=1
	endif

	;get values for outline
	for i=0ULL, out_xdim-1 do begin
		in_win = in_image[in_ind_x[i]-1:in_ind_x[i]+1,in_ind_y[i]-1:in_ind_y[i]+1]
		if (in_win[1,1] gt 0) then begin
			index = where(in_win gt 0, count)
			;if we have a single pixel, that borders at least 2 pixels of 0 (water) and is over 2 stdev from mean of the pixels that do have value, we drop it back to the mean of the bordering pixels with value.
			;this is to remove those artifacts of bright pixels next to water
			if (count le 7) then begin
				tmp_mean = mean(in_win[index])
				tmp_stdev = stddev(in_win[index])
				if (abs(in_win[1,1] - tmp_mean ) gt (1.5*tmp_stdev)) then in_win[1,1] = tmp_mean
			endif
			out_line[i] = total(in_win[index] * weights[index])/total(weights[index])
		endif
	endfor

	;set out of bounds pixels of background value
	if (xcount gt 0) then out_line[x_obb_index] = background_val
	if (ycount gt 0) then out_line[y_obb_index] = background_val

	writeu, 1, out_line
endfor

print, 'Done!'

close, 1

end
