;not stable
function thermo_model2, data, t0, nz, dz

  ;heat capacity of ice:
  openr, lun, "ice_thermal_properties.txt", /get_lun
  thice=fltarr(4, 15)
  readf, lun, thice
  free_lun, lun

  sz=size(data, /dim)

  ;tind=day_value(data[0, *].date, data[0, 0].date)*86400
  ;nt=long(tind[sz[1]-1]/dt)
  ;t=[findgen(nt)*dt, tind[sz[1]-1]]

  temp=fltarr(nz, sz[0], sz[1])
  temp[*, *, 0]=t0

  ;Euler step to begin with, then leapfrog scheme:
  rho=interpol(reform(thice[1, *]), reform(thice[0, *]), t0-273.2)*1000
  k=interpol(reform(thice[2, *]), reform(thice[0, *]), t0-273.2)
  cp=interpol(reform(thice[2, *]), reform(thice[0, *]), t0-273.2)

  dtdz=-(2*temp[1:nz-2, *, 0]-temp[0:nz-3, *, 0]-temp[2:nz-1, *, 0])/dz/dz

  ;q=calc_q(data[*, 0].t2m, data[*, 0].v2m, data[*, 0].rh, data[*, 0].cc, $
;	  	data[*, 0].sw, temp[0, *, 0])

  ;just assume that the surface equilibrates:
  tsurf=fltarr(sz[0])
  for j=0, sz[0]-1 do begin
    tsurf[j]=derive_surf_temp(dz, data[j, 0].t2m, data[j, 0].v2m, data[j, 0].rh, data[j, 0].cc, $
	  	data[j, 0].sw, tw=temp[1, j, 0])
  endfor

  dtdt=[transpose(fltarr(sz[0])), -dtdz*k/rho/cp, transpose(fltarr(sz[0]))]

  dt=day_value(data[0, 1].date, data[0, 0].date)*86400

  temp[*, *, 1]=temp[*, *, 0]+dtdt*dt
  temp[0, *, 1]=tsurf

  for i=1, sz[1]-1 do begin
    rho=interpol(reform(thice[1, *]), reform(thice[0, *]), temp[*, *, i]-273.2)*1000
    k=interpol(reform(thice[2, *]), reform(thice[0, *]), temp[*, *, i]-273.2)
    cp=interpol(reform(thice[2, *]), reform(thice[0, *]), temp[*, *, i]-273.2)

    dtdz=-(2*temp[1:nz-2, *, i]-temp[0:nz-3, *, i]-temp[2:nz-1, *, i])/dz/dz

    ;q=calc_q(t2m, v2m, rh, cc, qsw, temp[0, *, i])
    for j=0, sz[0]-1 do begin
      tsurf[j]=derive_surf_temp(dz, data[j, i].t2m, data[j, i].v2m, data[j, i].rh, data[j, i].cc, $
	  	data[j, i].sw, tw=temp[1, j, i])
    endfor

    dtdt=[transpose(fltarr(sz[0])), -dtdz*k/rho/cp, transpose(fltarr(sz[0]))]

    ind=where(finite(dtdt) eq 0, cnt)
    if cnt gt 0 then stop

    dt=day_value(data[0, i+1].date, data[i-1].date)*86400

    temp[*, *, i+1]=temp[*, *, i-1]+dtdt*dt
    temp[nz-1, *, i+1]=temp[nz-2, *, i+1]
    temp[0, *, i+1]=tsurf
  endfor

  return, temp
end


