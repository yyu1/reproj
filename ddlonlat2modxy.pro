;this procedure maps a lon lat in degrees to a pixel on the MODIS sinusoidal 1km 
;global grid   . returned x,y index are 0 based indices

PRO ddlonlat2modxy, ddlon, ddlat, modx, mody
  ;43200 x 21600
  xdim = 43200L
  ydim = 21600L
  pix_size = 926.625433d0   ;modis 1km pixel size in m
  PI = 3.141592653589793238d0
	R_earth = 6371007.181   ;MODIS sinusoidal sphere radius in m

  ;convert from degrees to radians
  rad_lon = ddlon / 180.0d0 * PI
  rad_lat = ddlat / 180.0d0 * PI

	modx = R_earth * rad_lon * cos(rad_lat)
	mody = R_earth * rad_lat

end

