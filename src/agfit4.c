#include "R.h"
#include "Rinternals.h"
#include <R_ext/Utils.h>  

#ifdef USING_R
/* typedef int Sint; */
#define S_EVALUATOR    /* Turn this into a "blank line" in R */
#else

typedef long Sint;
#ifdef  defineVar
#undef  defineVar
#endif
#define defineVar(a,b,c) ASSIGN_IN_FRAME(a,b, INTEGER_VALUE(c))

#ifdef  eval
#undef  eval
#endif
#define eval(a, b)  EVAL_IN_FRAME(a, INTEGER_VALUE(b))

#ifdef asInteger
#undef asInteger
#endif
#define asInteger(a) INTEGER_VALUE(a)
#ifdef asReal
#undef asReal
#endif
#define asReal(a) NUMERIC_VALUE(a)
#endif

#ifdef USING_R
#define ALLOC(a,b)  R_alloc(a,b)
#else
#define ALLOC(a,b)  S_alloc(a,b)
#endif

#ifdef USING_R
void cox_callback(int which, double *coef, double *first, double *second,
		  double *penalty, int *flag, int p, SEXP fexpr, SEXP rho);
#endif

#include <math.h>

double **dmatrix(double *array, int ncol, int nrow){
  int i;
  double **pointer;
  pointer = (double **) ALLOC(nrow, sizeof(double *));
  for (i=0; i<nrow; i++) {
    pointer[i] = array;
    array += ncol;
  }
  return(pointer);
}

int cholesky2(double **matrix, int n, double toler){
  double temp;
  int  i,j,k;
  double eps, pivot;
  int rank;
  int nonneg;
  nonneg=1;
  eps =0;
  for (i=0; i<n; i++) {
    if (matrix[i][i] > eps)  eps = matrix[i][i];
    for (j=(i+1); j<n; j++)  matrix[j][i] = matrix[i][j];
  }
  eps *= toler;
  //
  rank =0;
  for (i=0; i<n; i++) {
    pivot = matrix[i][i];
    if (pivot < eps) {
      matrix[i][i] =0;
      if (pivot < -8*eps) nonneg= -1;
    }
    else  {
      rank++;
      for (j=(i+1); j<n; j++) {
	temp = matrix[j][i]/pivot;
	matrix[j][i] = temp;
	matrix[j][j] -= temp*temp*pivot;
	for (k=(j+1); k<n; k++) matrix[k][j] -= temp*matrix[k][i];
      }
    }
  }
  return(rank * nonneg);
}

void chsolve2(double **matrix, int n, double *y){
  register int i,j;
  register double temp;
  for (i=0; i<n; i++) {
    temp = y[i] ;
    for (j=0; j<i; j++)
      temp -= y[j] * matrix[i][j] ;
    y[i] = temp ;
  }
  for (i=(n-1); i>=0; i--) {
    if (matrix[i][i]==0)  y[i] =0;
    else {
      temp = y[i]/matrix[i][i];
      for (j= i+1; j<n; j++)
	temp -= y[j]*matrix[j][i];
      y[i] = temp;
    }
  }
}

void chinv2(double **matrix , int n){
  register double temp;
  register int i,j,k;
  for (i=0; i<n; i++){
    if (matrix[i][i] >0) {
      matrix[i][i] = 1/matrix[i][i];   /*this line inverts D */
      for (j= (i+1); j<n; j++) {
	matrix[j][i] = -matrix[j][i];
	for (k=0; k<i; k++)     /*sweep operator */
	  matrix[j][k] += matrix[j][i]*matrix[i][k];
      }
    }
  }
  for (i=0; i<n; i++) {
    if (matrix[i][i]==0) {  /* singular row */
      for (j=0; j<i; j++) matrix[j][i]=0;
      for (j=i; j<n; j++) matrix[i][j]=0;
    }
    else {
      for (j=(i+1); j<n; j++) {
	temp = matrix[j][i]*matrix[j][j];
	if (j!=i) matrix[i][j] = temp;
	for (k=i; k<j; k++)
	  matrix[i][k] += temp*matrix[j][k];
      }
    }
  }
}


SEXP agfit4(SEXP surv2,      SEXP covar2,    SEXP strata2,
            SEXP weights2,  SEXP offset2,    SEXP ibeta2,
            SEXP sort12,     SEXP sort22,    SEXP method2,
            SEXP maxiter2,   SEXP  eps2,     SEXP tolerance2,
            SEXP doscale2) { 

    int i,j,k,person;
    int indx2, istrat, p;
    int ksave, nrisk, ndeath;
    int nused, nvar;
 
    double **covar, **cmat, **imat;  /*ragged array versions*/
    double *a, *oldbeta, *maxbeta;
    double *scale;
    double *a2, **cmat2;
    double *eta;
    double  denom, zbeta, risk;
    double  time;
    double  temp, temp2;
    double  newlk =0;
    int     halving;    /*are we doing step halving at the moment? */
    double  tol_chol, eps;
    double  meanwt;
    int itemp, deaths;
    double efron_wt, d2, meaneta;

    /* inputs */
    double *start, *stop, *event;
    double *weights, *offset;
    int *sort1, *sort2, maxiter;
    int *strata;
    double method;  /* saving this as double forces some double arithmetic */
    int doscale;

    /* returned objects */
    SEXP imat2, means2, beta2, u2, loglik2;
    double *beta, *u, *loglik, *means;
    SEXP sctest2, flag2, iter2;
    double *sctest;
    int *flag, *iter;
    SEXP rlist;
    static const char *outnames[]={"coef", "u", "imat", "loglik", "means",
                                   "sctest", "flag", "iter", ""};
    int nprotect;  /* number of protect calls I have issued */

    /* get sizes and constants */
    nused = nrows(covar2);
    nvar  = ncols(covar2);
    method= asInteger(method2);
    eps   = asReal(eps2);
    tol_chol = asReal(tolerance2);
    maxiter = asInteger(maxiter2);
    doscale = asInteger(doscale2);
  
    /* input arguments */
    start = REAL(surv2);
    stop  = start + nused;
    event = stop + nused;
    weights = REAL(weights2);
    offset = REAL(offset2);
    sort1  = INTEGER(sort12);
    sort2  = INTEGER(sort22);
    strata = INTEGER(strata2);

    /*
    ** scratch space
    **  nvar: a, a2, newbeta, maxbeta, scale
    **  nvar*nvar: cmat, cmat2
    **  n: eta
    */
    eta = (double *) R_alloc(nused + 5*nvar + 2*nvar*nvar, sizeof(double));
    a = eta + nused;
    a2 = a +nvar;
    maxbeta = a2 + nvar;
    scale  = maxbeta + nvar;
    oldbeta = scale + nvar;   

    /*
    **  Set up the ragged arrays
    **  covar2 might not need to be duplicated, even though
    **  we are going to modify it, due to the way this routine was
    **  was called.  In this case NAMED(covar2) will =0
    **  Changed NAMED(covar2) to MAYBE_REFERENCED(covar2) on 04/17/2018 
    **  as suggested by Tomas Kalibera
    */
    PROTECT(imat2 = allocVector(REALSXP, nvar*nvar));
    nprotect =1;
    if (MAYBE_REFERENCED(covar2)>0) {
        PROTECT(covar2 = duplicate(covar2)); 
        nprotect++;
        }
    covar= dmatrix(REAL(covar2), nused, nvar);
    imat = dmatrix(REAL(imat2),  nvar, nvar);
    cmat = dmatrix(oldbeta+ nvar,   nvar, nvar);
    cmat2= dmatrix(oldbeta+ nvar + nvar*nvar, nvar, nvar);

    /*
    ** create the output structures
    */
    PROTECT(rlist = mkNamed(VECSXP, outnames));
    nprotect++;
    beta2 = SET_VECTOR_ELT(rlist, 0, duplicate(ibeta2));
    beta  = REAL(beta2);

    u2 =    SET_VECTOR_ELT(rlist, 1, allocVector(REALSXP, nvar));
    u = REAL(u2);

    SET_VECTOR_ELT(rlist, 2, imat2);
    loglik2 = SET_VECTOR_ELT(rlist, 3, allocVector(REALSXP, 2)); 
    loglik  = REAL(loglik2);

    means2 = SET_VECTOR_ELT(rlist, 4, allocVector(REALSXP, nvar));
    means  = REAL(means2);

    sctest2 = SET_VECTOR_ELT(rlist, 5, allocVector(REALSXP, 1));
    sctest =  REAL(sctest2);
    flag2  =  SET_VECTOR_ELT(rlist, 6, allocVector(INTSXP, 1));
    flag   =  INTEGER(flag2);
    iter2  =  SET_VECTOR_ELT(rlist, 7, allocVector(INTSXP, 1));
    iter = INTEGER(iter2);
    
    /*
    ** Subtract the mean from each covar, as this makes the variance
    **  computation much more stable
    */
    temp2 =0;
    for (i=0; i<nused; i++) temp2 += weights[i];  /* sum of weights */

    for (i=0; i<nvar; i++) {
        maxbeta[i] = 0;  /* temporary, save the max abs covariate value */
        temp=0;
        for (person=0; person<nused; person++)
            temp += weights[i] * covar[i][person];
        temp /= temp2;
        means[i] = temp;
        for (person=0; person<nused; person++)
            covar[i][person] -=temp;
        if (doscale ==1) { /* also scale the regression */
            temp =0;
            for (person=0; person<nused; person++)
                temp += weights[person] * fabs(covar[i][person]);
            if (temp >0) temp = temp2/temp;
            else temp = 1.0;  /* rare case of a constant covariate */
            scale[i] = temp;
            for (person=0; person<nused; person++) {
                covar[i][person] *= temp;
                if (fabs(covar[i][person]) > maxbeta[i])
                    maxbeta[i] = fabs(covar[i][person]);
            }
        } else {
            /* scaling is only turned off during debugging 
               still, cover the case */
            for (person=0; person<nused; person++) {
                if (fabs(covar[i][person]) > maxbeta[i])
                    maxbeta[i] = fabs(covar[i][person]);
            }
        }
    }
 
    if (doscale ==1) {
        for (i=0; i<nvar; i++) beta[i] /= scale[i]; /* rescale initial betas */
        }
    else {for (i=0; i<nvar; i++) scale[i] = 1.0;}
    
    /* 
    **        Set the max beta.  For safety we want beta*x < log(max double) so 
    **  the the risk score exp(x beta) never becomes either infinite or 0.
    **  This limit is around 700 for most hardware.  Since exp(23) > population
    **  of the earth, any beta*x over 20 is a silly relative risk for a Cox
    **  model, however.  
    **   We want to cut off huge values, but not take action very often since
    **  doing so can mess up the iteration in general.
    **  One of the case-cohort papers suggests using anoffset of -100 to 
    **  indicate "no risk", meaning that x*beta values of 50-100 can occur 
    **  in "ok" data sets.  Compromise.
    */
    for (i=0; i<nvar; i++) maxbeta[i] = 200/maxbeta[i];

    ndeath =0;
    for (i=0; i<nused; i++) ndeath += event[i];
    
    /* First iteration, which has different ending criteria */
    for (i=0; i<nvar; i++) {
        u[i] =0;
        a[i] =0;
        for (j=0; j<nvar; j++) {
            imat[i][j] =0 ;
            cmat[i][j] =0;
        }
    }

    for (person=0; person<nused; person++) {
        zbeta = 0;      /* form the term beta*z   (vector mult) */
        for (i=0; i<nvar; i++)
            zbeta += beta[i]*covar[i][person];
        eta[person] = zbeta + offset[person];
    }

    /*
    **  'person' walks through the the data from 1 to n,
    **     sort1[0] points to the largest stop time, sort1[1] the next, ...
    **  'time' is a scratch variable holding the time of current interest
    **  'indx2' walks through the start times.  It will be smaller than 
    **    'person': if person=27 that means that 27 subjects have stop >=time,
    **    and are thus potential members of the risk set.  If 'indx2' =9,
    **    that means that 9 subjects have start >=time and thus are NOT part
    **    of the risk set.  (stop > start for each subject guarrantees that
    **    the 9 are a subset of the 27). 
    **  Basic algorithm: move 'person' forward, adding the new subject into
    **    the risk set.  If this is a new, unique death time, take selected
    **    old obs out of the sums, add in obs tied at this time, then
    **    add terms to the loglik, etc.
    */
    istrat=0;
    indx2 =0;
    denom =0;
    meaneta =0;
    nrisk =0;
    newlk =0;
    for (person=0; person<nused;) {
        p = sort1[person];
        if (event[p]==0){
            nrisk++;
            meaneta += eta[p];
            risk = exp(eta[p]) * weights[p];
            denom += risk;
            for (i=0; i<nvar; i++) {
                a[i] += risk*covar[i][p];
                for (j=0; j<=i; j++)
                    cmat[i][j] += risk*covar[i][p]*covar[j][p];
                }
            person++;
            /* nothing more needs to be done for this obs */
        }
        else {
            time = stop[p];
            /*
            ** subtract out the subjects whose start time is to the right
            */
            for (; indx2<strata[istrat]; indx2++) {
                p = sort2[indx2];
                if (start[p] < time) break;
                nrisk--;
                meaneta -= eta[p];
                risk = exp(eta[p]) * weights[p];
                denom -= risk;
                for (i=0; i<nvar; i++) {
                    a[i] -= risk*covar[i][p];
                    for (j=0; j<=i; j++)
                        cmat[i][j] -= risk*covar[i][p]*covar[j][p];
                    }
                }

            /*
            ** compute the averages over subjects with
            **   exactly this death time (a2 & c2)
            ** (and add them into a and cmat while we are at it).
            */
            efron_wt =0;
            meanwt =0;
            for (i=0; i<nvar; i++) {
                a2[i]=0;
                for (j=0; j<nvar; j++) {
                    cmat2[i][j]=0;
                    }
                }
            deaths=0;
            for (k=person; k<strata[istrat]; k++) {
                p = sort1[k];
                if (stop[p] < time) break;
                risk = exp(eta[p]) * weights[p];
                denom += risk;
                nrisk++;
                meaneta += eta[p];

                for (i=0; i<nvar; i++) {
                    a[i] += risk*covar[i][p];
                    for (j=0; j<=i; j++)
                        cmat[i][j] += risk*covar[i][p]*covar[j][p];
                    }
                if (event[p]==1) {
                    deaths += event[p];
                    efron_wt += risk*event[p];
                    meanwt += weights[p];
                    for (i=0; i<nvar; i++) {
                        a2[i]+= risk*covar[i][p];
                        for (j=0; j<=i; j++)
                            cmat2[i][j] += risk*covar[i][p]*covar[j][p];
                        }
                }
                }
            ksave = k;
                
            /* 
            ** If the average eta value has gotton out of hand, fix it.
            ** We must avoid overflow in the exp function (~750 on Intel)
            ** and want to act well before that, but not take action very often.  
            ** One of the case-cohort papers suggests an offset of -100 meaning
            ** that etas of 50-100 can occur in "ok" data, so make it larger than this.
            */
            if (fabs(meaneta) > (nrisk *110)) {  
                meaneta = meaneta/nrisk;
                for (i=0; i<nused; i++) eta[i] -= meaneta;
                temp = exp(-meaneta);
                denom *= temp;
                for (i=0; i<nvar; i++) {
                    a[i] *= temp;
                    a2[i] *= temp;
                    for (j=0; j<nvar; j++) {
                        cmat[i][j]*= temp;
                        cmat2[i][j] *= temp;
                    }
                }
                meaneta =0;
            }
                
            /*
            ** Add results into u and imat for all events at this time point
            */
            meanwt /= deaths;
            itemp = -1;
            for (; person<ksave; person++) {
                p = sort1[person];
                if (event[p]==1) {
                    itemp++;
                    temp = itemp*method/(double) deaths;
                    d2 = denom - temp*efron_wt;
                    newlk +=  weights[p]*eta[p] -meanwt *log(d2);

                    for (i=0; i<nvar; i++) {
                        temp2 = (a[i] - temp*a2[i])/d2;
                        u[i] += weights[p]*covar[i][p] - meanwt*temp2;
                        for (j=0; j<=i; j++)
                            imat[j][i] += meanwt* (
                                        (cmat[i][j] - temp*cmat2[i][j])/d2-
                                           temp2*(a[j]-temp*a2[j])/d2);
                        }
                    }
                }
        }

        if (person == strata[istrat]) {
            istrat++;
            denom =0;
            meaneta=0;
            nrisk =0;
            indx2 = person;
            for (i=0; i<nvar; i++) {
                a[i] =0;
                for (j=0; j<nvar; j++) {
                    cmat[i][j]=0;
                    }
                }
        }
    }   /* end  of accumulation loop */
    loglik[0] = newlk;   /* save the loglik for iteration zero  */
    loglik[1] = newlk;

    /* Calculate the score test */
    for (i=0; i<nvar; i++) /*use 'a' as a temp to save u0, for the score test*/
        a[i] = u[i];
    *flag = cholesky2(imat, nvar, tol_chol);
    chsolve2(imat,nvar,a);        /* a replaced by  a *inverse(i) */
    *sctest=0;
    for (i=0; i<nvar; i++)
        *sctest +=  u[i]*a[i];

    if (maxiter ==0) {
        *iter =0;
        loglik[1] = newlk;
        chinv2(imat, nvar);
        for (i=0; i<nvar; i++) {
            beta[i] *= scale[i];  /* return to original scale */
            u[i] /= scale[i];
            imat[i][i] *= scale[i] * scale[i];
            for (j=0; j<i; j++) {
                imat[j][i] *= scale[i] * scale[j];
                imat[i][j] = imat[j][i];
            }
        }
        UNPROTECT(nprotect);
        return(rlist);
    }
    else {  
        /* Update beta for the next iteration
        **  Never complain about convergence on this first step or impose step
        **  halving.  That way someone can force one iter at a time.
        */
        for (i=0; i<nvar; i++) {
            oldbeta[i] = beta[i];
            beta[i] = beta[i] + a[i];
        }
    }
    /* main loop */
    halving =0 ;             /* =1 when in the midst of "step halving" */
    for (*iter=1; *iter<= maxiter; (*iter)++) {
        for (i=0; i<nvar; i++) {
            u[i] =0;
            a[i] =0;
            for (j=0; j<nvar; j++) {
                imat[i][j] =0 ;
                cmat[i][j] =0;
            }
        }

        for (person=0; person<nused; person++) {
            zbeta = 0;      /* form the term beta*z   (vector mult) */
            for (i=0; i<nvar; i++)
                zbeta += beta[i]*covar[i][person];
            eta[person] = zbeta + offset[person];
        }

        /*
        **  'person' walks through the the data from 1 to n,
        **     sort1[0] points to the largest stop time, sort1[1] the next, ...
        **  'time' is a scratch variable holding the time of current interest
        **  'indx2' walks through the start times.  It will be smaller than 
        **    'person': if person=27 that means that 27 subjects have stop >=time,
        **    and are thus potential members of the risk set.  If 'indx2' =9,
        **    that means that 9 subjects have start >=time and thus are NOT part
        **    of the risk set.  (stop > start for each subject guarrantees that
        **    the 9 are a subset of the 27). 
        **  Basic algorithm: move 'person' forward, adding the new subject into
        **    the risk set.  If this is a new, unique death time, take selected
        **    old obs out of the sums, add in obs tied at this time, then
        **    add terms to the loglik, etc.
        */
        istrat=0;
        indx2 =0;
        denom =0;
        meaneta =0;
        nrisk =0;
        newlk =0;
        for (person=0; person<nused;) {
            p = sort1[person];
            if (event[p]==0){
                nrisk++;
                meaneta += eta[p];
                risk = exp(eta[p]) * weights[p];
                denom += risk;
                for (i=0; i<nvar; i++) {
                    a[i] += risk*covar[i][p];
                    for (j=0; j<=i; j++)
                        cmat[i][j] += risk*covar[i][p]*covar[j][p];
                    }
                person++;
                /* nothing more needs to be done for this obs */
            }
            else {
                time = stop[p];
                /*
                ** subtract out the subjects whose start time is to the right
                */
                for (; indx2<strata[istrat]; indx2++) {
                    p = sort2[indx2];
                    if (start[p] < time) break;
                    nrisk--;
                    meaneta -= eta[p];
                    risk = exp(eta[p]) * weights[p];
                    denom -= risk;
                    for (i=0; i<nvar; i++) {
                        a[i] -= risk*covar[i][p];
                        for (j=0; j<=i; j++)
                            cmat[i][j] -= risk*covar[i][p]*covar[j][p];
                        }
                    }

                /*
                ** compute the averages over subjects with
                **   exactly this death time (a2 & c2)
                ** (and add them into a and cmat while we are at it).
                */
                efron_wt =0;
                meanwt =0;
                for (i=0; i<nvar; i++) {
                    a2[i]=0;
                    for (j=0; j<nvar; j++) {
                        cmat2[i][j]=0;
                        }
                    }
                deaths=0;
                for (k=person; k<strata[istrat]; k++) {
                    p = sort1[k];
                    if (stop[p] < time) break;
                    risk = exp(eta[p]) * weights[p];
                    denom += risk;
                    nrisk++;
                    meaneta += eta[p];

                    for (i=0; i<nvar; i++) {
                        a[i] += risk*covar[i][p];
                        for (j=0; j<=i; j++)
                            cmat[i][j] += risk*covar[i][p]*covar[j][p];
                        }
                    if (event[p]==1) {
                        deaths += event[p];
                        efron_wt += risk*event[p];
                        meanwt += weights[p];
                        for (i=0; i<nvar; i++) {
                            a2[i]+= risk*covar[i][p];
                            for (j=0; j<=i; j++)
                                cmat2[i][j] += risk*covar[i][p]*covar[j][p];
                            }
                    }
                    }
                ksave = k;
                    
                /* 
                ** If the average eta value has gotton out of hand, fix it.
                ** We must avoid overflow in the exp function (~750 on Intel)
                ** and want to act well before that, but not take action very often.  
                ** One of the case-cohort papers suggests an offset of -100 meaning
                ** that etas of 50-100 can occur in "ok" data, so make it larger than this.
                */
                if (fabs(meaneta) > (nrisk *110)) {  
                    meaneta = meaneta/nrisk;
                    for (i=0; i<nused; i++) eta[i] -= meaneta;
                    temp = exp(-meaneta);
                    denom *= temp;
                    for (i=0; i<nvar; i++) {
                        a[i] *= temp;
                        a2[i] *= temp;
                        for (j=0; j<nvar; j++) {
                            cmat[i][j]*= temp;
                            cmat2[i][j] *= temp;
                        }
                    }
                    meaneta =0;
                }
                    
                /*
                ** Add results into u and imat for all events at this time point
                */
                meanwt /= deaths;
                itemp = -1;
                for (; person<ksave; person++) {
                    p = sort1[person];
                    if (event[p]==1) {
                        itemp++;
                        temp = itemp*method/(double) deaths;
                        d2 = denom - temp*efron_wt;
                        newlk +=  weights[p]*eta[p] -meanwt *log(d2);

                        for (i=0; i<nvar; i++) {
                            temp2 = (a[i] - temp*a2[i])/d2;
                            u[i] += weights[p]*covar[i][p] - meanwt*temp2;
                            for (j=0; j<=i; j++)
                                imat[j][i] += meanwt* (
                                            (cmat[i][j] - temp*cmat2[i][j])/d2-
                                               temp2*(a[j]-temp*a2[j])/d2);
                            }
                        }
                    }
            }

            if (person == strata[istrat]) {
                istrat++;
                denom =0;
                meaneta=0;
                nrisk =0;
                indx2 = person;
                for (i=0; i<nvar; i++) {
                    a[i] =0;
                    for (j=0; j<nvar; j++) {
                        cmat[i][j]=0;
                        }
                    }
            }
        }   /* end  of accumulation loop */

        *flag = cholesky2(imat, nvar, tol_chol);
        if (fabs(1-(loglik[1]/newlk))<= eps  && halving==0){ /* all done */
            loglik[1] = newlk;
            chinv2(imat, nvar);
            for (i=0; i<nvar; i++) {
                beta[i] *= scale[i];  /* return to original scale */
                u[i] /= scale[i];
                imat[i][i] *= scale[i] * scale[i];
                for (j=0; j<i; j++) {
                    imat[j][i] *= scale[i] * scale[j];
                    imat[i][j] = imat[j][i];
                }
            }
            UNPROTECT(nprotect);
            return(rlist);
        }

        if (*iter < maxiter) { /*update beta */
            if (newlk < loglik[1])   {    /*it is not converging ! */
                halving =1;
                for (i=0; i<nvar; i++)
                    beta[i] = (oldbeta[i] + beta[i]) /2; /*half of old increment */
            }
            else {
                halving=0;
                loglik[1] = newlk;
                chsolve2(imat,nvar,u);

                for (i=0; i<nvar; i++) {
                    oldbeta[i] = beta[i];
                    beta[i] = beta[i] +  u[i];
                    if (beta[i]> maxbeta[i]) beta[i] = maxbeta[i];
                    else if (beta[i] < -maxbeta[i]) beta[i] = -maxbeta[i];
                }
            }
        }  
        R_CheckUserInterrupt();  /* be polite -- did the user hit cntrl-C? */
    } /*return for another iteration */
    loglik[1] = newlk;
    chinv2(imat, nvar);
    for (i=0; i<nvar; i++) {
        beta[i] *= scale[i];  /* return to original scale */
        u[i] /= scale[i];
        imat[i][i] *= scale[i] * scale[i];
        for (j=0; j<i; j++) {
            imat[j][i] *= scale[i] * scale[j];
            imat[i][j] = imat[j][i];
        }
    }
    UNPROTECT(nprotect);
    return(rlist);
}
