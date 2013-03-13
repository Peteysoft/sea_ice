function create_salinity_prof, x, bulks, type
  case type of
    'C': begin
		c0=11.63
		c1=-18.46
		c2=-1.53
		c3=18.78
         end
    'S': begin
		c0=11.24
		c1=-32.33
		c2=56.82
		c3=-30.42
         end
    '?': begin
		c0=1.89
		c1=29.19
		c2=-77.49
		c3=56.38
	 end
    'I': begin
		c0=8.32
		c1=-4.16
		c2=-11.28
		c3=8.80
         end
  endcase

  ;normalize the curves so that integrated bulk salinity is 1.:
  norm=c0+c1/2.+c2/3.+c3/4.
  ;print, norm

  nz=n_elements(x)

  ;assume ascending order:
  d=x[nz-1]

  x1=x/d
  x2=x1*x1
  x3=x2*x1

  s=bulks*(c0+c1*x1+c2*x2+c3*x3)/norm

  return, s

end

