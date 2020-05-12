#include <math.h>
#include <getopt.h>

#include "mcc_ice_model.h"

#define DATELEN 25
#define LINELEN 200

int main (int argc, char **argv) {
  long seed;
  FILE * lin;

  char *topofile;
  char *initfile;

  char line[LINELEN];

  mcc_ice_model *ice;

  //parameters:
  real_t ti;		//ice temp.
  real_t tw;		//water temp.
  real_t sw;		//sea salinity
  real_t nu;		//frequency
  real_t ts;		//sky temp.

  real_t *x0, *y0;	//initial conditions
  real_t *zangle;
  real_t *offset;

  long ninit;		//number of initial conditions
  long niter;		//number of iterations

  real_t *s, *h, *z;	//ice surface
  long np;

  int ncon;
  char c;

  real_t fov;

  real_t *tbv;
  real_t *tbh;

  real_t *verr;
  real_t *herr;

  complex_t *eps;
  real_t epsp, epspp;

  char * logfile=NULL;

  int err;

  //exit(0);
  //
  
  fov=13.16;
  niter=2000;

  //set defaults to match our particular problem:
  ti=272.;		//ice temperature
  tw=272.;		//water temperature
  ts=5.;		//sky temperature
  sw=5.;		//sea salinity
  nu=1.4;		//frequency

  while ((c = getopt(argc, argv, "v:w:y:s:t:f:n:a:")) != -1) {
    switch (c) {
      case ('n'):
        ncon=sscanf(optarg, "%ld", &niter);
        if (ncon != 1) {
          fprintf(stderr, "Warning: garbled command option argument: -n %s", optarg);
	  exit(2);
        }
        break;
      case ('f'):
        ncon=sscanf(optarg, "%lg", &fov);
        if (ncon != 1) {
          fprintf(stderr, "Warning: garbled command option argument: -f %s", optarg);
	  exit(2);
        }
        break;
      case ('w'):
        ncon=sscanf(optarg, "%lg", &tw);
        if (ncon != 1) {
          fprintf(stderr, "Warning: garbled command option argument: -w %s", optarg);
	  exit(2);
        }
        break;
      case ('t'):
        ncon=sscanf(optarg, "%lg", &ti);
        if (ncon != 1) {
          fprintf(stderr, "Warning: garbled command option argument: -t %s", optarg);
	  exit(2);
        }
        break;
      case ('y'):
        ncon=sscanf(optarg, "%lg", &ts);
        if (ncon != 1) {
          fprintf(stderr, "Warning: garbled command option argument: -y %s", optarg);
	  exit(2);
        }
        break;
      case ('s'):
        ncon=sscanf(optarg, "%lg", &sw);
        if (ncon != 1) {
          fprintf(stderr, "Warning: garbled command option argument: -s %s", optarg);
	  exit(2);
        }
        break;
      case ('v'):
        ncon=sscanf(optarg, "%lg", &nu);
        if (ncon != 1) {
          fprintf(stderr, "Warning: garbled command option argument: -s %s", optarg);
	  exit(2);
        }
        break;
      case ('a'):
	logfile=optarg;
        break;
      case ('?'):
             fprintf(stderr, "Unknown option: -%c -- ignored\n", optopt);
             break; 
      default:
             fprintf(stderr, "Error parsing command line\n");
             exit(2);
             break;
    }
  } 
    
  argc-=optind;
  argv+=optind;

  if (argc != 2) {
    printf("\nUsage: ice_main [-n niter] [-f fov] topofile initfile\n\n");
    printf("where:\n");
    printf("  topofile   ascii file containing ice topography\n");
    printf("  initfile   ascii file containing initial conditions\n\n");
    printf("options:\n");
    printf("  -n         number of iterations (%ld)\n", niter);
    printf("  -f         instrument fov (%5.1f) [deg.]\n", fov);
    printf("  -v         frequency (%5.1f) [GHz]\n", nu);
    printf("  -t         ice temperature (%6.1f) [GHz]\n", ti);
    printf("  -w         water temperature (%6.1f) [K]\n", tw);
    printf("  -s         water salinity (%5.1f) [psu]\n", sw);
    printf("  -y         sky temperature (%4.1f) [K]\n", ts);
    printf("  -a         name of log file\n");
    printf("\n\n");
    exit(1);
  }

  topofile=argv[0];
  initfile=argv[1];

  lin=fopen(topofile, "r");
  fgets(line, LINELEN, lin);
  sscanf(line, "%ld", &np);

  s=new real_t[np];
  h=new real_t[np];
  z=new real_t[np];

  eps=new complex_t[np];

  for (long i=0; i<np; i++) {
    fgets(line, LINELEN, lin);
    //printf("%s\n", line);
    sscanf(line, "%lf %lf %lf %lf %lf", s+i, h+i, z+i, &epsp, &epspp);
    eps[i]=complex_t(epsp, epspp);
    //printf("%g %g %g\n", s[i], h[i], z[i]);
  }
  fclose(lin);

  lin=fopen(initfile, "r");
  fgets(line, LINELEN, lin);
  sscanf(line, "%ld", &ninit);

  x0=new real_t[ninit];
  y0=new real_t[ninit];
  zangle=new real_t[ninit];
  offset=new real_t[ninit];

  for (long i=0; i<ninit; i++) {
    fgets(line, LINELEN, lin);
    sscanf(line, "%lg %lg %lg %lg", x0+i, y0+i, zangle+i, offset+i);
  }
  fclose(lin);

  ice=new mcc_ice_model();

  if (logfile!=NULL) ice->setlogfile(logfile);
  ice->set_fov(fov);

  ice->set_parm(nu, ti, tw, sw, ts);

  err=ice->set_surface(s, h, eps, np, s, z, eps, np);
  if (err != 0) {
    exit(err);
  }

  //ice->print_fit("pissfit.txt", 0.01);

  tbv=new real_t[ninit];
  tbh=new real_t[ninit];

  verr=new real_t[ninit];
  herr=new real_t[ninit];

  for (long i=0; i<ninit; i++) {
    tbv[i]=ice->tb(x0[i], y0[i], zangle[i], offset[i], 0, niter, verr[i]);
    tbh[i]=ice->tb(x0[i], y0[i], zangle[i], offset[i], 1, niter, herr[i]);

    printf("%f +/- %f\n", tbv[i], verr[i]);
    printf("%f +/- %f\n", tbh[i], herr[i]);

  }

  delete ice;

  delete [] tbv;
  delete [] tbh;

  delete [] verr;
  delete [] herr;

  delete [] s;
  delete [] h;
  delete [] z;

  delete [] x0;
  delete [] y0;
  delete [] zangle;
  delete [] offset;

  delete [] eps;

}
    

