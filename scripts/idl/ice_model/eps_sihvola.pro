
;effective permittivity mixture formulas from Sihvola

function eps_Sihvola, eps1, eps2, v, dfac

  eps=eps1+v*(eps2-eps1)*eps1/(eps1+dfac*(eps2-eps1))/(1-dfac*v*(eps2-eps1)/(eps1+(eps2-eps1))*dfac)

  return, eps

end

