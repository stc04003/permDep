#include <R.h>
#include <Rmath.h>
#include <math.h>

void mysample(int *x, int *n, int *N, 
	      //output
	      double *out) {
  int t = 0; 
  int m = 0; 
  double u;
  while (m < *n) {
    u = runif(0, 1);
    if ( (*N - t) * u >= *n - m ) { t += 1; }
    else {
      out[m] = x[t];
      t += 1; 
      m += 1;
    }
  }
  out;
}

void kendallTrun(double *x, double *y, double *Delta, int *n, 
		 // output
		 double *out) {
  int i, j, k;
  double *bb = Calloc(*n * (*n - 1), double);
  double M = 0.0;
  double v = 0.0;
  double v1 = 0.0;
  double v2 = 0.0;
  double h12 = 0.0;
  double mu;
  double tmp;
  for (i = 0; i < (*n - 1); i++) {
    for (j = i + 1; j < *n; j++) {
      if (fmax(x[i], x[j]) <= fmin(y[i], y[j]) &&
	  Delta[i] *  Delta[j] + Delta[i] * (y[i] <= y[j]) + Delta[j] * (y[j] <= y[i]) > 0) {
	tmp = (y[i] - y[j]) * (x[i] - x[j]);
	h12 += (tmp > 0) - (tmp < 0);
	bb[i * (*n - 1) + j - 1] = ((tmp > 0) - (tmp < 0)) / sqrt((*n - 1) * (*n - 2));
	bb[j * (*n - 1) + i] = ((tmp > 0) - (tmp < 0)) / sqrt((*n - 1) * (*n - 2));
	M += 1;
      }
    }
  }
  for (i = 0; i < *n; i++) {
    for (j = 0; j < (*n - 1); j++) {
      v1 += bb[i * (*n - 1) + j];
      v2 += bb[i * (*n - 1) + j] * bb[i * (*n - 1) + j];
    }
    v += (v1 * v1 - v2) / *n;
    v1 = 0.0; 
    v2 = 0.0;
  }
  mu = 2 * M / (*n * (*n - 1)); 
  out[0] = h12 / M;
  // out[1] = 4 * (v / (*n * (*n - 1) * (*n - 2))) / (mu * mu * *n);
  out[1] = 4 * v / (mu * mu * *n);
  out[2] = M;
  Free(bb);
  out;
}


/* Jing's version with minor modifications*/
/* There are problems with se as it directly implements \xi rather than using formula from Sec. 6 */

/* void kendallTrun(double *x, double *y, double *Delta, int *n,  */
/* 		 // output */
/* 		 double *out) { */
/*   double *bb = Calloc(*n, double);  */
/*   int i, j, k; */
/*   for (i = 0; i < *n; i++) { */
/*     bb[i] = 0; */
/*   } */
/*   double h12 = 0.0; */
/*   double M = 0.0; */
/*   double n2 = 0.0; */
/*   double tau; */
/*   double tmp; */
/*   double bbSum; */
/*   for (i = 0; i < (*n - 1); i++) { */
/*     for (j = i + 1; j < *n; j++) { */
/*       //      if (fmax(x[i], x[j]) <= fmin(y[i], y[j]) & Delta[i] * Delta[j] + Delta[i] * (y[i] < y[j]) + Delta[j] * (y[j] < y[i]) > 0) { */
/*       if (fmax(x[i], x[j]) <= fmin(y[i], y[j]) & Delta[i] *  Delta[j] + Delta[i] * (y[i] < y[j]) + Delta[j] * (y[j] < y[i]) > 0) { */
/* 	tmp = (y[i] - y[j]) * (x[i] - x[j]); */
/* 	h12 += (tmp > 0) - (tmp < 0); */
/* 	bb[i] += (tmp > 0) - (tmp < 0); */
/* 	bb[j] += (tmp > 0) - (tmp < 0); */
/* 	M += 1; */
/*       } */
/*       n2 += 1; */
/*     } */
/*   } */
/*   tau = h12 / M; */
/*   for (i = 0; i < *n; i++) {  */
/*     bbSum += bb[i] * bb[i]; */
/*   } */
/*   out[0] = tau; */
/*   out[1] = 1/(2 * n2 * (*n-2)) * (bbSum - 2 * M) - M * M * tau * tau / (n2 * n2); */
/*   out[1] = 4 * out[1] / (*n * M * M / (n2 * n2)); */
/*   out[2] = M; */
/*   Free(bb); */
/*   out; */
/* } */


void kendallTrunWgt(double *x, double *y, double *Delta, int *n, double *scX, double *scT, 
		    // output
		    double *out) {
  double h12 = 0.0;
  int i, j, k;
  double M = 0.0;
  double n2 = 0.0;
  double tau;
  double tmp;
  double wgt;
  for (i = 0; i < (*n - 1); i++) {
    for (j = i + 1; j < *n; j++) {
      if (fmax(x[i], x[j]) <= fmin(y[i], y[j]) &&
	  Delta[i] *  Delta[j] + Delta[i] * (y[i] <= y[j]) + Delta[j] * (y[j] <= y[i]) > 0 &&
	  scX[i] * scX[j] * scT[i] * scT[j] > 0) {
	wgt = scX[i] * scX[j] / (scT[i] * scT[j]);
	tmp = (y[i] - y[j]) * (x[i] - x[j]);
	h12 += ((tmp > 0) - (tmp < 0)) / wgt;
	M += 1 / wgt;
      }
    }
  }
  tau = h12 / M;
  out[0] = tau;
  out[1] = M;
  out;
}

