@time_lib

function get_skin_temp, date1, lon1, lat1
  common get_skin_temp_buf, nc, year, lon, lat, t, tind_old, slice1, slice2

  date=convert_date(date1)
  if n_elements(year) eq 0 then begin
    year=-1
    tind_old=-1
  endif
  if date.year ne year then begin
    if n_elements(nc) eq 1 then ncdf_close, nc
    year=date.year

    yrstr=string(year, format="(i4.4)")
    filename="/freax/storage/home/pmills/smos/data/skt.sfc.gauss."+yrstr+".nc"
    nc=ncdf_open(filename)

    ncdf_varget, nc, "lon", lon
    ncdf_varget, nc, "lat", lat
    ncdf_varget, nc, "time", t
  endif

  nlon=n_elements(lon)
  nlat=n_elements(lat)

  refdate=convert_date('1/1/14')
  ;refdate={time_str, year, 1, 1, 0, 0, 0}
  tval=day_value(date, refdate)*24.

  nt=n_elements(t)
  tindd=interpol(findgen(nt), t, tval)
  tind=long(tindd)

  if tind ne tind_old then begin
    ncdf_attget, nc, "skt", "add_offset", add_offset
    ncdf_attget, nc, "skt", "scale_factor", scale_factor
    ncdf_varget, nc, "skt", slice1, offset=[0, 0, tind], count=[nlon, nlat, 1]
    ncdf_varget, nc, "skt", slice2, offset=[0, 0, tind+1], count=[nlon, nlat, 1]
    slice1=add_offset+scale_factor*slice1
    slice2=add_offset+scale_factor*slice2
    tind_old=tind
  endif

  xind=interpol(findgen(nlon), lon, lon1)
  yind=interpol(findgen(nlat), lat, lat1)

  skt1=interpolate(slice1, xind, yind)
  skt2=interpolate(slice2, xind, yind)

  frac=tindd-tind
  skt=skt1*(1-frac)+skt2*frac

  ;stop

  return, skt

end
    

