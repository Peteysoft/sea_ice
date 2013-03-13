;@ridged_ice_model2
@ridged_ice_model3

resolve_all

path="/freax/storage/home/pmills/smos/campaign_data/Police2007/"

s=0.65
tice=-3.

vb=s*(-49.185/tice+0.532)/1000
;eps=complex(3.2+9.*vb, 0.03+4.*vb)
eps=complex(4., 0.1)

stop

;goto, skip
doplot=0
ntest=5000
niter=10000

  tfile=path+"ice_thickness/20070313_P2A_P4X.dat"
  lfile=path+"laser/200703131539.fb"
  rfile=path+"corrected/20070313_TrackAndCircles/P4X_to_P2A/07217460.r31"
  eps=complex(4.0, 0.1)
  ridged_ice_model, rfile, tfile, lfile, "model_results/model_run_P4X_to_P2A_const_eps=4.0+0.1i.idlsav", $
		ntest=ntest, eps_ice0=eps, tice=273.2+tice, doplot=doplot, niter=niter
;  eps=complex(4.6, 0.43)
;  ridged_ice_model, rfile, tfile, lfile, "model_results/model_run_P4X_to_P2A_const_eps=4.6+0.43i.idlsav", $
;		ntest=ntest, eps_ice=eps, tice=273.2+tice, doplot=doplot
  ridged_ice_model, rfile, tfile, lfile, "model_results/model_run_P4X_to_P2A.parm_p5.idlsav", $
		ntest=ntest, tice=273.2+tice, eps_parm0="eps_parm_p5", doplot=doplot, $
		niter=niter
;  ridged_ice_model, rfile, tfile, lfile, "model_run_P4X_to_P2A.dvary.idlsav", $
;		ntest=ntest, tice=273.2+tice, /dvary;, /doplot
;skip:

end


