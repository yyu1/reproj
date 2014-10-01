;Main program for reprojecting an image.

;Input file setup - normally, we shouldn't read everything into memory
;But the files here should be under 10 gig in size and the 16GB memory can handle it
;This makes reading individual pixel values much faster

in_xdim = 43200ULL;  Africa
in_ydim = 38400ULL;  Africa
in_file = '/Volumes/Global_250m/output/afr/v3/maxent_afr_hlorey_v3_filled.flt' ;Africa
in_tl_x = -3335851.5588D     ;top left corner of top left pixel of input image  (in MODIS grid m)
in_tl_y = 4447802.0784D     ;top left corner of top left pixel of input image  (in MODIS grid m)
in_x_size = 231.65635825D ;pixel width in meters
in_y_size = 231.65635825D ;pixel height in meters


;Target file setup
out_xdim = 28500ULL;  Africa
out_ydim = 30000ULL;  Africa
out_file = '/Volumes/Global_250m/output/afr/v3/maxent_afr_hlorey_v3_filled_9.6sec.flt' ;Africa
out_tl_x = -20.0D    ;top left corner of top left pixel of input image  (in degrees longitude)
out_tl_y = 40.0D     ;top left corner of top left pixel of input image  (in degrees latitude)
out_x_size = 0.0026666666666D ;pixel width degrees longitude
out_y_size = 0.0026666666666D ;pixel height degrees latitude
background_val = 0


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

for j=0ULL, out_ydim-1 do begin
	if (j mod 10000 eq 0) then print, 'line ', j

	;set up values for middle of pixel coordinates of current output line
	out_xcoords[*] = out_tl_x + out_x_size*xcoeff
	out_ycoords[*] = out_tl_y - (j+0.5)*out_y_size

	ddlonlat2modxy, out_xcoords, out_ycoords, in_xcoords, in_ycoords

	;convert in_xcoords and  in_ycoords to image indices
	in_ind_x[*] = floor((in_xcoords-in_tl_x)/in_x_size, /l64)
	in_ind_y[*] = floor((in_tl_y-in_ycoords)/in_y_size, /l64)

	;Look for out of bounds pixels
	x_obb_index = where((in_ind_x lt 0) or (in_ind_x ge in_xdim), xcount)
	y_obb_index = where((in_ind_y lt 0) or (in_ind_y ge in_ydim), ycount)

	;temporarily setting out of bounds indices to 0,0
	if (xcount gt 0) then begin
		in_ind_x[x_obb_index]=0
	endif
	if (ycount gt 0) then begin
		in_ind_y[y_obb_index]=0
	endif

	;get values for outline
	for i=0ULL, out_xdim-1 do out_line[i] = in_image[in_ind_x[i],in_ind_y[i]]

	;set out of bounds pixels of background value
	if (xcount gt 0) then out_line[x_obb_index] = background_val
	if (ycount gt 0) then out_line[y_obb_index] = background_val

	writeu, 1, out_line
endfor

print, 'Done!'

close, 1

end
