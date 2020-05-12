#include <stdio.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_min.h>

typedef double real;

struct Vant_param {
  int n;
  real *freq;
  real *ar;		//Vant etal mixture model
  real *ai;
  real *br;
  real *bi;
  //fitted parameters:
  real c;
  real d;
  real e;
  real f;
};

double Debye_cost(double tau, void *param) {
  Vant_param *p2=(Vant_param *) param;
  real fp[p2->n];
  gsl_multifit_linear_workspace *work;
  gsl_matrix *cov;
  gsl_matrix *A;
  gsl_vector *b;
  gsl_vector *x;
  double chisq;
  int err;

  for (int i=0; i<p2->n; i++) fp[i]=p2->freq[i]*tau;

  A=gsl_matrix_alloc(p2->n*4, 4);
  gsl_matrix_set_zero(A);
  b=gsl_vector_alloc(p2->n*4);

  for (int i=0; i<p2->n; i++) {
    real fp2=fp[i]*fp[i];
    gsl_vector_set(b, i*4, p2->ar[i]*(1+fp2));
    gsl_vector_set(b, i*4+1, p2->br[i]*(1+fp2));
    gsl_vector_set(b, i*4+2, p2->ai[i]*(1+fp2));
    gsl_vector_set(b, i*4+3, p2->bi[i]*(1+fp2));
    gsl_matrix_set(A, i*4, 0, 1.);
    gsl_matrix_set(A, i*4, 2, 1.);
    gsl_matrix_set(A, i*4+1, 1, 1.);
    gsl_matrix_set(A, i*4+1, 3, 1.);
    gsl_matrix_set(A, i*4+2, 2, fp[i]);
    gsl_matrix_set(A, i*4+3, 3, fp[i]);
  }

  /*
  for (int i=0; i<p2->n*4; i++) {
    printf("%lg %lg %lg %lg\n", gsl_matrix_get(A, i, 0), gsl_matrix_get(A, i, 1),
		    gsl_matrix_get(A, i, 2), gsl_matrix_get(A, i, 3));
  }
  printf("\n");

  for (int i=0; i<p2->n*4; i++) {
    printf("%lg\n", gsl_vector_get(b, i));
  }
  */

  x=gsl_vector_alloc(4);
  cov=gsl_matrix_alloc(4, 4);
  work=gsl_multifit_linear_alloc(p2->n*4, 4);

  err=gsl_multifit_linear(A, b, x, cov, &chisq, work);

  p2->c=gsl_vector_get(x, 0);
  p2->d=gsl_vector_get(x, 1);
  p2->e=gsl_vector_get(x, 2);
  p2->f=gsl_vector_get(x, 3);

  gsl_multifit_linear_free(work);
  gsl_matrix_free(A);
  gsl_vector_free(b);
  gsl_vector_free(x);
  gsl_matrix_free(cov);

  printf("%lg\n", chisq);

  return chisq;
}

int main(int argc, char **argv) {
  FILE *fs;
  Vant_param param;
  real tau1, tau2, tau3;
  gsl_min_fminimizer *minimizer;
  gsl_function F;
  double a, b, m;
  int iter=0;
  int maxiter=1000;
  int status;
  double tol=0.001;

  fs=fopen(argv[1], "r");
  fscanf(fs, "%d", &param.n);
  param.freq=new real[param.n];
  param.ar=new real[param.n];
  param.ai=new real[param.n];
  param.br=new real[param.n];
  param.bi=new real[param.n];
  for (int i=0; i<param.n; i++) {
    fscanf(fs, "%lf %lf %lf %lf %lf", param.freq+i, param.ar+i, param.ai+i, param.br+i, param.bi+i);
  }
  fclose(fs);
  param.c=0;
  param.d=0;
  param.e=0;
  param.f=0;

  tau1=atof(argv[2]);
  tau3=atof(argv[3]);

  tau2=(tau1+tau3)/2;

  F.function=&Debye_cost;
  F.params=&param;

  minimizer=gsl_min_fminimizer_alloc(gsl_min_fminimizer_brent);
  gsl_min_fminimizer_set(minimizer, &F, tau2, tau1, tau3);

  do {
    iter++;
    status = gsl_min_fminimizer_iterate (minimizer);

    m = gsl_min_fminimizer_x_minimum (minimizer);
    a = gsl_min_fminimizer_x_lower (minimizer);
    b = gsl_min_fminimizer_x_upper (minimizer);

    status = gsl_min_test_interval (a, b, 0.0, tol);

    if (status == GSL_SUCCESS) printf ("Converged:\n");

    printf ("%5d [%.7f, %.7f] %.7f %.7f\n",
              iter, a, b, m, b-a);
  } while (status == GSL_CONTINUE && iter < maxiter);

  printf("\n");
  printf("%lg %lg %lg %lg\n", param.c, param.d, param.e, param.f);
  printf("\n");

  for (int i=0; i<param.n; i++) {
    real denom=1+param.freq[i]*param.freq[i]*m*m;
    printf("%lg %lg %lg %lg %lg\n", param.freq[i], param.c+param.e/denom, param.d+param.f/denom,
		    param.e*param.freq[i]*m/denom, param.f*param.freq[i]*m/denom);
  }

  gsl_min_fminimizer_free(minimizer);
  delete [] param.freq;
  delete [] param.ar;
  delete [] param.ai;
  delete [] param.br;
  delete [] param.bi;

}

