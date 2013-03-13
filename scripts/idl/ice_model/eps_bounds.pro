@ice_model_util

pro eps_bounds, nu, t, s, eps_f1, eps_f2, $
		nf=nf, ave=ave, std=std, ntrial=ntrial, $
		eps=eps

    if n_elements(nf) ne 1 then nf=100

    abs0=abs_zero()

    t1=t-abs0

    ;use more accurate formula for brine volume:
    v=brine_volume_from_salinity(t1, s, sb)
    sb=sb[0]
    v=v[0]

    ;eps1=complex(3.1884+9.1e-4*t, 0.)
    ;again, more accurate formula for eps of pure ice (here it probably doesn't make a difference...):
    eps1=calc_eps_pi(nu, t1)   

    eps2=eps_br(nu, t1, sb)
	    
    z=v*findgen(nf)/(nf-1)
    c1=(1-v)/(1/(1-eps1/eps2)-z)
    eps_f1=eps2*(1-c1)
    z=(1-v)*findgen(nf)/(nf-1)
    c2=v/(1/(1-eps1/eps2)-z)
    eps_f2=eps1/(1-c2)
			  
    if arg_present(ave) or arg_present(std) or arg_present(eps) then begin
      if n_elements(ntrial) ne 1 then ntrial=100
      df1=eps_f1[1:nf-1]-eps_f1[0:nf-2]
      df2=eps_f2[1:nf-1]-eps_f2[0:nf-2]
      ds1=sqrt(float(df1)^2+imaginary(df1)^2)
      ds2=sqrt(float(df2)^2+imaginary(df2)^2)
      s1=[0, total(ds1, /cum)]
      s2=[0, total(ds2, /cum)]

      eps=complexarr(ntrial)
      sran=(s1[nf-1]+s2[nf-1])*randomu(seed, ntrial)
      indg=where(sran gt s1[nf-1], cntg, complement=indl, ncompl=cntl)
      if cntg gt 0 then begin
        eps[indg]=interpol(eps_f2, s2, sran[indg]-s1[nf-1])
      endif
      if cntl gt 0 then begin
        eps[indl]=interpol(eps_f1, s1, sran[indl])
      endif

      ave=total(eps)/ntrial
      std=stddev(eps)
    endif

end

