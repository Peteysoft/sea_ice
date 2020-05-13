#include <stdio.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_min.h>

typedef double real;

struct Vant_param2 {
  int n;
  real *freq;
  real *ar;		//Vant etal mixture model
  real *ai;
  //fitted parameters:
  real c;
  real d;
};

//fit each component, Vb and non-Vb, separately:
double Debye_cost2(double tau, void *param) {
  Vant_param2 *p2=(Vant_param2 *) param;
  real fp[p2->n];
  gsl_multifit_linear_workspace *work;
  gsl_matrix *cov;
  gsl_matrix *A;
  gsl_vector *b;
  gsl_vector *x;
  double chisq;
  int err;

  for (int i=0; i<p2->n; i++) fp[i]=p2->freq[i]*tau;

  A=gsl_matrix_alloc(p2->n*2, 2);
  gsl_matrix_set_zero(A);
  b=gsl_vector_alloc(p2->n*2);

  for (int i=0; i<p2->n; i++) {
    real fp2=fp[i]*fp[i];
    gsl_vector_set(b, i*2, p2->ar[i]*(1+fp2));
    gsl_vector_set(b, i*2+1, p2->ai[i]*(1+fp2));
    gsl_matrix_set(A, i*2, 0, 1.);
    gsl_matrix_set(A, i*2, 1, 1.);
    gsl_matrix_set(A, i*2+1, 1, fp[i]);
  }

  /*
  for (int i=0; i<p2->n*2; i++) {
    printf("%lg %lg\n", gsl_matrix_get(A, i, 0), gsl_matrix_get(A, i, 1));
  }
  printf("\n");

  for (int i=0; i<p2->n*2; i++) {
    printf("%lg\n", gsl_vector_get(b, i));
  }
  */

  x=gsl_vector_alloc(2);
  cov=gsl_matrix_alloc(2, 2);
  work=gsl_multifit_linear_alloc(p2->n*2, 2);

  err=gsl_multifit_linear(A, b, x, cov, &chisq, work);

  p2->c=gsl_vector_get(x, 0);
  p2->d=gsl_vector_get(x, 1);

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
  int n;
  real **data;
  Vant_param2 param;
  real tau1, tau2, tau3;
  gsl_min_fminimizer *minimizer;
  gsl_function F;
  double a, b, m;
  int iter=0;
  int maxiter=1000;
  int status;
  double tol=0.001;

  fs=fopen(argv[1], "r");
  fscanf(fs, "%d", &n);
  data=new real *[n];
  data[0]=new real[n*5];
  for (int i=0; i<n; i++) {
    data[i]=data[0]+i*5;
    fscanf(fs, "%lf %lf %lf %lf %lf", data[i], data[i]+1, data[i]+2, data[i]+3, data[i]+4);
  }
  fclose(fs);
  param.c=0;
  param.d=0;

  tau1=atof(argv[2]);
  tau3=atof(argv[3]);

  tau2=(tau1+tau3)/2;

  param.n=n;
  param.freq=new real[n];
  param.ar=new real[n];
  param.ai=new real[n];

  F.function=&Debye_cost2;
  F.params=&param;

  for (int i=0; i<param.n; i++) {
    printf("%lf %lf %lf %lf %lf\n", data[i][0], data[i][1], data[i][2],
		    data[i][3], data[i][4]);
    param.freq[i]=data[i][0];
    param.ar[i]=data[i][1];
    param.ai[i]=data[i][2];
  }

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
  printf("%lg %lg\n", param.c, param.d);
  printf("\n");

  for (int i=0; i<param.n; i++) {
    real denom=1+param.freq[i]*param.freq[i]*m*m;
    printf("%lg %lg %lg\n", param.freq[i], param.c+param.d/denom,
		    param.d*param.freq[i]*m/denom);
  }

  for (int i=0; i<param.n; i++) {
    param.ar[i]=data[i][3];
    param.ai[i]=data[i][4];
  }

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
  printf("%lg %lg\n", param.c, param.d);
  printf("\n");

  for (int i=0; i<param.n; i++) {
    real denom=1+param.freq[i]*param.freq[i]*m*m;
    printf("%lg %lg %lg\n", param.freq[i], param.c+param.d/denom,
		    param.d*param.freq[i]*m/denom);
  }


  gsl_min_fminimizer_free(minimizer);
  delete [] param.freq;
  delete [] param.ar;
  delete [] param.ai;

}

