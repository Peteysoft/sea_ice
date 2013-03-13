@time_lib
@wv_formulas

function dtempdt, t, temp
  common dtempdt_common, data, thice, tind, dz, nz, npt

  rho0=interpol(reform(thice[1, *]), reform(thice[0, *]), temp[*, *]-273.2)
  ;exponential density profile
  ;from Rist et al. (2002) JGR 107 (B1) 10.1029/2000JB000058
  rho=918.-559.*exp(-dz*findgen(nz)/32.5)
  rho=rebin(rho, nz, npt)

  ;density relative to maximum:
  rf=rho/rho0

  k=interpol(reform(thice[2, *]), reform(thice[0, *]), temp[*, *]-273.2)
  cp=interpol(reform(thice[3, *]), reform(thice[0, *]), temp[*, *]-273.2)

  temp=reform(temp, nz, npt, /overwrite)

  dtdz=-(2*temp[1:nz-2, *]-temp[0:nz-3, *]-temp[2:nz-1, *])/dz/dz

  nt=n_elements(tind)
  tint=interpol(findgen(nt), tind, t)
  tintl=long(tint)
  if tintl ge nt-1 then tintl=nt-2
  tfrac=tint-tintl
  t2m=(1-tfrac)*data[*, tintl].t2m+tfrac*data[*, tintl+1].t2m
  v2m=(1-tfrac)*data[*, tintl].v2m+tfrac*data[*, tintl+1].v2m
  rh=(1-tfrac)*data[*, tintl].rh+tfrac*data[*, tintl+1].rh
  cc=(1-tfrac)*data[*, tintl].cc+tfrac*data[*, tintl+1].cc
  qsw=(1-tfrac)*data[*, tintl].sw+tfrac*data[*, tintl+1].sw

  q=calc_q(t2m, v2m, rh, cc, qsw, temp[0, *])

  dtdt=[transpose(q/rf/cp/dz), dtdz*k/rf/cp, transpose(fltarr(npt))]/1000000.

  ind=where(finite(dtdt) ne 1, cnt)
  if cnt ne 0 then stop
  ;stop

  return, dtdt

end


function thermo_model, data, t0, nz, dz, dt, flt_bnd=flt, $
		save_t=save_t, date=date
  common dtempdt_common, data1, thice, tind, dz1, nz1, npt

  dz1=dz
  nz1=nz
  data1=data

  ;heat capacity of ice:
  openr, lun, "ice_thermal_properties.txt", /get_lun
  thice=fltarr(4, 15)
  readf, lun, thice
  free_lun, lun

  sz=size(data, /dim)
  npt=sz[0]

  ;tind=day_value(data[0, *].date, data[0, 0].date)*86400
  tind=time_class_diff(data[0, *].date, data[0, 0].date)*86400
  nt=long(tind[sz[1]-1]/dt)+1
  t=[findgen(nt)*dt, tind[sz[1]-1]]

  nd=n_elements(date)
  if nd eq sz[0] then begin
    ;intind=day_value(date, data[0, 0].date)*86400
    intind=time_class_diff(date, data[0, 0].date)*86400
    save_t=fltarr(nz, nd)
  endif

  ;tswp=ptr_new()
  temp1=t0

  for i=1, nt-1 do begin
    print, i
    dtdt=dtempdt(t[i-1], temp1)

    temp2=rk4(temp1, dtdt, t[i-1], dt, "dtempdt")

    ;stop 

    ind=where(finite(temp1) eq 0, cnt)
    if cnt gt 0 then stop

    ;dt=day_value(data[0, i+1].date, data[i-1].date)*86400

    ;temp[*, *, i+1]=temp[*, *, i-1]+dtdt*dt
    ;floating boundary conditions:
    if keyword_set(flt) then temp2[nz-1, *]=temp2[nz-2, *]
    ;temp[0, *, 1]=tsurf

    if nd eq sz[0] then begin
      for j=0, nd-1 do begin
        if intind[j] gt t[i-1] and intind[j] gt t[i] then begin
       	  lind=long(intind[j])
	  frac=intind[j]-lind
          save_t[*, j]=temp1[*, j]*(1-frac)+frac*temp1[*, j]
	endif
      endfor
    endif

    temp1=temp2

  endfor

  return, temp1
end


