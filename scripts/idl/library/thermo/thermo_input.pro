@time_lib
@wv_formulas

;get the input to the thermodynamic model:

function thermo_input, t0, lon, lat, nt
  ;use 2 m height:
  h=2.
  alpha=0.7
  solar=1366.

  sixhour=convert_date('6:0')
  date=time_array(t0, sixhour, nt)
  np=n_elements(lon)

  data=replicate({date:{time_str}, rh:0., t2m:0., v2m:0., cc:0., p2m:0., sw:0.}, $
	  	np, nt)

  for i=0L, nt-1 do begin
    print, time_string(date[i])
    data[*, i].date=replicate(date[i], np)
    get_ecmwf_field, date[i], "ZZ", zz, longrid, latgrid
    nlon=n_elements(longrid)
    nlat=n_elements(latgrid)
    nz=(size(zz, /dim))[2]
    get_ecmwf_field, date[i], "TT", tt
    get_ecmwf_field, date[i], "UU", uu
    get_ecmwf_field, date[i], "VV", vv
    get_ecmwf_field, date[i], "TT", tt
    get_ecmwf_field, date[i], "PP", pp
    get_ecmwf_field, date[i], "SQ", sq
    get_ecmwf_field, date[i], "CC", cc

    lonind=interpol(findgen(nlon), longrid, lon)
    latind=interpol(findgen(nlat), latgrid, lat)
    
    z1=interpolate(zz[*, *, nz-1], lonind, latind)
    z2=interpolate(zz[*, *, nz-2], lonind, latind)
    zind=h/(z2-z1)

    t1=interpolate(tt[*, *, nz-1], lonind, latind)
    t2=interpolate(tt[*, *, nz-2], lonind, latind)
    data[*, i].t2m=t1+(t2-t1)*zind

    v1=interpolate(vv[*, *, nz-1], lonind, latind)
    v2=interpolate(vv[*, *, nz-2], lonind, latind)
    v2m=v1+(v2-v1)*zind

    u1=interpolate(uu[*, *, nz-1], lonind, latind)
    u2=interpolate(uu[*, *, nz-2], lonind, latind)
    u2m=u1+(u2-u1)*zind

    data[*, i].v2m=sqrt(v2m*v2m+u2m*u2m)

    p1=interpolate(pp[*, *, nz-1], lonind, latind)
    p2=interpolate(pp[*, *, nz-2], lonind, latind)
    p2m=p1+(p2-p1)*zind
    data[*, i].p2m=p2m

    sq1=interpolate(sq[*, *, nz-1], lonind, latind)
    sq2=interpolate(sq[*, *, nz-2], lonind, latind)
    sq2m=sq1+(sq2-sq1)*zind

    ccprof=interpolate(transpose(cc), latind, lonind)
    data[*, i].cc=1
    for j=0, nz-1 do begin
      data[*, i].cc=data[*, i].cc*(1-ccprof[j, *])
    endfor
    data[*, i].cc=1-data[*, i].cc

    data[*, i].rh=100*sq1*1.6*p2m/eq_vp(data[*, i].t2m)

    data[*, i].sw=solar*sw_flux(lat, date[i])*(1-0.62*data[*, i].cc)*(1-alpha)

  endfor

return, data

end

