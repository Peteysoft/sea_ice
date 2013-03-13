@ice_model_util

;implements "brine correction:" because of decreasing temperatures in
;the ice sheet, some of the brine will be ejected
;because temperature is a monotonic function of time (for the const.
;weather model) this can be calculated from initial and final temperatures
;only

pro brine_correction, d, s, g, ts, s1, skim, slost, weight=weight, $
		temp=temp, stnew=st1, stwnew=st1a, time=t
  ;crude approximation for sea temp.
  ;(should really use sea salinity...):
  tw=271.

  nt=n_elements(d)

  dz=d-[0, d[0:nt-2]]

  s2=rebin(transpose(s), nt, nt)
  dz2=rebin(transpose(dz), nt, nt)
  ;age of each layer:
  t=[0, total(dz/g, /cum)]

  ;memory is cheap:
  temp=fltarr(nt, nt)
  for tind=0, nt-1 do begin
    ;temperature of each layer as a function of time:
    temp[tind, 0:tind]=ts[tind]+(tw-ts[tind])*d[0:tind]/d[tind]
  endfor

  ;brine salinity:
  temp_c=temp-abs_zero()
  sb=brine_salinity(temp_c)

  ;ice density:
  rhoi=0.917-1.403e-4*temp_c
  ;brine density:
  rhob=1.+0.0008*sb

  ;brine volume:
  vb=s2*rhoi/(sb*rhob+s2*(rhoi-rhob))

  ;total density:
  rho=vb*rhob+(1-vb)*rhoi

  ;mass loss for each layer:
  rho0=rebin(transpose(rho[lindgen(nt), lindgen(nt)]), nt, nt)
  ;relative volume loss:
  dv=rho0/rho-1

  ;need to compare this with the brine volume:
  ind=where(dv gt vb, cnt, complement=ind2)

  ;mass loss:
  dm2=fltarr(nt, nt)
  if cnt gt 0 then begin
    dm1=fltarr(nt, nt)
    vbnew=fltarr(nt, nt)
    vbnew[ind2]=vb[ind2]-dv[ind2]
    dm1[ind]=rhob[ind]*dz2[ind]*vb[ind]
    dm2[ind]=dz2[ind]*rhoi[ind]*(dv[ind]-vb[ind])
    dm1[ind2]=rhob[ind2]*dz2[ind2]*dv[ind2]
  endif else begin
    dm1=rhob*dz2*dv
    vbnew=vb-dv
  endelse
  dm=dm1+dm2

  ;total mass loss for each time step:
  dmt=fltarr(nt)
  slost=fltarr(nt)	;average salinity of lost brine
  if keyword_set(weight) then begin
    w=1-findgen(nt)/(nt-1)
    for tind=0, nt-1 do begin
      dmt[tind]=total(w*dm[tind, 0:tind])
      slost[tind]=total(w*dm1[tind, 0:tind]*sb[tind, 0:tind])/dmt[tind]
    endfor
  endif else begin
    for tind=0, nt-1 do begin
      dmt[tind]=total(dm[tind, 0:tind])
      slost[tind]=total(dm1[tind, 0:tind]*sb[tind, 0:tind])/dmt[tind]
    endfor
  endelse

  ;should work out the whole formula at some point...
  s1=vbnew*sb*rhob/(vbnew*(rhob-rhoi)+rhoi)

  ;thickness of top layer:
  ts_c=ts-abs_zero()
  sbskim=brine_salinity(ts_c)
  rhob_skim=1.+0.0008*sbskim
  rhoi_skim=0.917-1.403e-4*ts_c
  vbskim=slost*rhoi/(sbskim*rhob_skim+slost*(rhoi_skim-rhob_skim))
  rhoskim=vbskim*rhob_skim+(1-vbskim)*rhoi_skim
  if keyword_set(weight) then begin
    skim=dmt/rhoskim
  endif else begin
    skim=dmt/rhoskim/2
  endelse

  ;join salinity of top layer and the rest:
  ;s1=[[slost], [s1]]
  indb=where(s1 lt 0, cntb)
  if cntb gt 0 then print, "err"

  st1a=fltarr(nt)
  st1=fltarr(nt)

  for tind=0, nt-1 do begin
    dz_tmp=dz[0:tind]
    st1a[tind]=total(s1[tind, *]*dz_tmp)/total(dz_tmp)
    st1[tind]=(skim[tind]*slost[tind]+total(s1[tind, *]*dz_tmp))/(total(dz_tmp)+skim[tind])
  endfor

  ;stop

end

