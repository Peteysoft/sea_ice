@time_lib

function get_ice_conc, t1, clon, clat, hemi, path=path, $
		xind=xind, yind=yind, gridpath=gridpath

  ;common get_ice_con_buf, ice, longrid, latgrid, t0

  if n_elements(path) ne 1 then path="."
  if strmid(path, strlen(path)-1, 1) ne '/' then path=path+'/'
  
  if n_elements(gridpath) ne 1 then $
	  	gridpath="/vapor/pmills/Police/"
  if strmid(gridpath, strlen(gridpath)-1, 1) ne '/' then gridpath=gridpath+'/'
  
  t=convert_date(t1)
  tstring=string(t.day+100L*(t.month+100L*t.year), format="(i8.8)")
  
  if hemi eq 1 then begin
    grid_file="LongitudeLatitudeGrid-n6250-Arctic.hdf"
    asi_file="asi-n6250-"+tstring+"-v5.hdf"
  endif else begin
    grid_file="LongitudeLatitudeGrid-s6250-Antarctic.hdf"
    asi_file="asi-s6250-"+tstring+"-v5.hdf"
  endelse

  ;read the gridding for the sea-ice concentrations:
  id=hdf_sd_start(gridpath+grid_file)
  sd=hdf_sd_select(id, 0)
  hdf_sd_getdata, sd, longrid
  sd=hdf_sd_select(id, 1)
  hdf_sd_getdata, sd, latgrid
  hdf_sd_end, id

  ;read the sea-ice concentrations:
  id=hdf_sd_start(path+asi_file)
  sd=hdf_sd_select(id, 0)
  hdf_sd_getdata, sd, ice
  hdf_sd_end, id

  s=size(ice, /dim)
  nx=s[0]
  ny=s[1]

  mapll, x0, y0, latgrid[0, 0], longrid[0, 0]+45, hemi
  mapll, x1, y1, latgrid[nx-1, ny-1], longrid[nx-1, ny-1]+45, hemi

  mapll, x, y, clat, clon+45, hemi

  xind=(x-x0)*(nx-1)/(x1-x0)
  yind=(y-y0)*(ny-1)/(y1-y0)

  con=interpolate(ice, xind, yind)

  ;stop

  return, con

end

