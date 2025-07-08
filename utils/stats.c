/**+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

This file is part of FORCE - Framework for Operational Radiometric 
Correction for Environmental monitoring.

Copyright (C) 2013-2022 David Frantz

FORCE is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

FORCE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with FORCE.  If not, see <http://www.gnu.org/licenses/>.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**/

/** The t-score calculations were obtained from SurfStat Australia's 
+++ t-distribution Calculator, online available at:
+++ https://surfstat.anu.edu.au/surfstat-home/tables/t.php
+++ t-distribution Calculator (C) 2017 Keith Dear
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**/

/**+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
This file contains functions for statistics
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**/


#include "stats.h"


/** One-pass variance and covariance estimation
+++ This function implements a one-pass estimation of variance and covari-
+++ ance based on recurrence formulas. It can be used to estimate mean of 
+++ x and y, variance or standard deviation of x and y, covariance between
+++ x and y. This function can also be used to compute one-pass linear re-
+++ gressions. Use this function in a loop.
+++-----------------------------------------------------------------------
+++ P. P�bay. SANDIA REPORT SAND2008-6212 (2008). Formulas for Robust, 
+++ One-Pass Parallel Computation of Co- variances and Arbitrary-Order 
+++ Statistical Moments.
+++-----------------------------------------------------------------------
--- x:      current x-value
--- y:      current y-value
--- mx:     last estimate of mean of x (is updated)
--- my:     last estimate of mean of y (is updated)
--- vx:     last estimate of variance of x (is updated)
--- vy:     last estimate of variance of y (is updated)
--- cv:     last estimate of covariance between x and y (is updated)
--- n:      number of observations
+++ Return: void
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**/
void covar_recurrence(double   x, double   y, 
                      double *mx, double *my, 
                      double *vx, double *vy, 
                      double *cv, double n){
double oldmx = *mx, oldmy = *my;
double oldvx = *vx, oldvy = *vy;
double oldcv = *cv;

  *mx = oldmx + (x-oldmx)/n;
  *my = oldmy + (y-oldmy)/n;
  *vx = oldvx + (x-oldmx)*(x-*mx);
  *vy = oldvy + (y-oldmy)*(y-*my);
  *cv = oldcv + (n-1)/n*(x-oldmx)*(y-oldmy);

  return;
}


/** One-pass covariance estimation
+++ This function implements a one-pass estimation of covariance based on
+++ recurrence formulas. It can be used to estimate mean of x and y, and 
+++ covariance between x and y. Use this function in a loop.
+++-----------------------------------------------------------------------
+++ P. P�bay. SANDIA REPORT SAND2008-6212 (2008). Formulas for Robust, 
+++ One-Pass Parallel Computation of Co- variances and Arbitrary-Order 
+++ Statistical Moments.
+++-----------------------------------------------------------------------
--- x:      current x-value
--- y:      current y-value
--- mx:     last estimate of mean of x (is updated)
--- my:     last estimate of mean of y (is updated)
--- cv:     last estimate of covariance between x and y (is updated)
--- n:      number of observations
+++ Return: void
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**/
void cov_recurrence(double   x, double   y, 
                    double *mx, double *my, 
                    double *cv, double n){
double oldmx = *mx, oldmy = *my;
double oldcv = *cv;

  *mx = oldmx + (x-oldmx)/n;
  *my = oldmy + (y-oldmy)/n;
  *cv = oldcv + (n-1)/n*(x-oldmx)*(y-oldmy);

  return;
}


/** One-pass skewness and kurtosis estimation
+++ This function implements a one-pass estimation of skewness and kurto-
+++ sis based on recurrence formulas. It can be used to estimate mean of 
+++ x, variance or standard deviation of x, skewness of x, kurtosis of x. 
+++ Use this function in a loop.
+++-----------------------------------------------------------------------
+++ P. P�bay. SANDIA REPORT SAND2008-6212 (2008). Formulas for Robust, 
+++ One-Pass Parallel Computation of Co- variances and Arbitrary-Order 
+++ Statistical Moments.
+++-----------------------------------------------------------------------
--- x:      current x-value
--- mx:     last estimate of mean of x (is updated)
--- vx:     last estimate of variance of x (is updated)
--- sx:     last estimate of skewness of x (is updated)
--- kx:     last estimate of kurtosis of x (is updated)
--- n:      number of observations
+++ Return: void
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**/
void kurt_recurrence(double   x, double *mx, 
                     double *vx, double *sx,
                     double *kx, double n){
double delta, delta_n, delta_n2, tmp;

  delta = x-(*mx);
  delta_n = delta/n;
  delta_n2 = delta_n*delta_n;
  tmp = delta*delta_n*(n-1);

  *mx = *mx + delta_n;
  *kx = *kx + tmp*delta_n2*(n*n-3*n+3) + 6*delta_n2*(*vx) - 4*delta_n*(*sx);
  *sx = *sx + tmp*delta_n*(n-2) - 3*delta_n*(*vx);
  *vx  = *vx + tmp;

  return;
}


/** One-pass skewness and kurtosis estimation
+++ This function implements a one-pass estimation of skewness based on 
+++ recurrence formulas. It can be used to estimate mean of x, variance or
+++ standard deviation of x, skewness of x. Use this function in a loop.
+++-----------------------------------------------------------------------
+++ P. P�bay. SANDIA REPORT SAND2008-6212 (2008). Formulas for Robust, 
+++ One-Pass Parallel Computation of Co- variances and Arbitrary-Order 
+++ Statistical Moments.
+++-----------------------------------------------------------------------
--- x:      current x-value
--- mx:     last estimate of mean of x (is updated)
--- vx:     last estimate of variance of x (is updated)
--- sx:     last estimate of skewness of x (is updated)
--- n:      number of observations
+++ Return: void
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**/
void skew_recurrence(double   x, double *mx, 
                     double *vx, double *sx,
                     double n){
double delta, delta_n, tmp;

  delta = x-(*mx);
  delta_n = delta/n;
  tmp = delta*delta_n*(n-1);

  *mx = *mx + delta_n;
  *sx = *sx + tmp*delta_n*(n-2) - 3*delta_n*(*vx);
  *vx  = *vx + tmp;

  return;
}


/** One-pass variance estimation
+++ This function implements a one-pass estimation of variance based on 
+++ recurrence formulas. It can be used to estimate mean of x, variance or
+++ standard deviation of x. Use this function in a loop.
+++-----------------------------------------------------------------------
+++ P. P�bay. SANDIA REPORT SAND2008-6212 (2008). Formulas for Robust, 
+++ One-Pass Parallel Computation of Co- variances and Arbitrary-Order 
+++ Statistical Moments.
+++-----------------------------------------------------------------------
--- x:      current x-value
--- mx:     last estimate of mean of x (is updated)
--- vx:     last estimate of variance of x (is updated)
--- n:      number of observations
+++ Return: void
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**/
void var_recurrence(double   x, double *mx, 
                    double *vx, double n){
double oldmx = *mx;
double oldvx = *vx;

  *mx = oldmx + (x-oldmx)/n;
  *vx = oldvx + (x-oldmx)*(x-*mx);

  return;
}


/** Compute kurtosis
+++ This function computes kurtosis based on the estimates of the recur-
+++ rence formulas above.
+++-----------------------------------------------------------------------
+++ P. P�bay. SANDIA REPORT SAND2008-6212 (2008). Formulas for Robust, 
+++ One-Pass Parallel Computation of Co- variances and Arbitrary-Order 
+++ Statistical Moments.
+++-----------------------------------------------------------------------
--- var:    recurrence estimate of variance
--- kurt:   recurrence estimate of kurtosis
--- n:      number of observations
+++ Return: kurtosis
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**/
double kurtosis(double var, double kurt, double n){

  return(kurt/(n*variance(var, n+1)*variance(var, n+1)));
}


/** Compute skewness
+++ This function computes skewness based on the estimates of the recur-
+++ rence formulas above.
+++-----------------------------------------------------------------------
+++ P. P�bay. SANDIA REPORT SAND2008-6212 (2008). Formulas for Robust, 
+++ One-Pass Parallel Computation of Co- variances and Arbitrary-Order 
+++ Statistical Moments.
+++-----------------------------------------------------------------------
--- var:    recurrence estimate of variance
--- skew:   recurrence estimate of skewness
--- n:      number of observations
+++ Return: skewness
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**/
double skewness(double var, double skew, double n){

  return(skew/(n*pow(standdev(var, n+1),3)));
}


/** Compute variance
+++ This function computes variance based on the estimates of the recur-
+++ rence formulas above.
+++-----------------------------------------------------------------------
+++ P. P�bay. SANDIA REPORT SAND2008-6212 (2008). Formulas for Robust, 
+++ One-Pass Parallel Computation of Co- variances and Arbitrary-Order 
+++ Statistical Moments.
+++-----------------------------------------------------------------------
--- var:    recurrence estimate of variance
--- n:      number of observations
+++ Return: variance
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**/
double variance(double var, double n){

  return(var/(n-1));
}


/** Compute standard deviation
+++ This function computes standard deviation based on the estimates of 
+++ the recurrence formulas above.
+++-----------------------------------------------------------------------
+++ P. P�bay. SANDIA REPORT SAND2008-6212 (2008). Formulas for Robust, 
+++ One-Pass Parallel Computation of Co- variances and Arbitrary-Order 
+++ Statistical Moments.
+++-----------------------------------------------------------------------
--- var:    recurrence estimate of variance
--- n:      number of observations
+++ Return: standard deviation
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**/
double standdev(double var, double n){

  return(sqrt(variance(var, n)));
}


/** Compute covariance
+++ This function computes covariance based on the estimates of the recur-
+++ rence formulas above.
+++-----------------------------------------------------------------------
+++ P. P�bay. SANDIA REPORT SAND2008-6212 (2008). Formulas for Robust, 
+++ One-Pass Parallel Computation of Co- variances and Arbitrary-Order 
+++ Statistical Moments.
+++-----------------------------------------------------------------------
--- cov:    recurrence estimate of covariance
--- n:      number of observations
+++ Return: covariance
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**/
double covariance(double cov, double n){

  return(cov/(n-1));
}


/** Compute slope of linear regression
+++ This function computes the slope of a linear regression.
--- cov:    covariance between x and y
--- varx:   variance of x
--- slope:  slope (returned)
+++ Return: void
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**/
void linreg_slope(double cov, double varx, double *slope){

  *slope = cov/varx;

  return;
}


/** Compute intercept of linear regression
+++ This function computes the intercept of a linear regression.
--- slope:     slope of the regression
--- mx:        mean of x
--- my:        mean of y
--- intercept: intercept (returned)
+++ Return:    void
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**/
void linreg_intercept(double slope, double mx, double my, double *intercept){

  *intercept = my-slope*mx;

  return;
}


/** Compute regression coefficients of linear regression
+++ This function computes the slope and intercept of a linear regression.
--- mx:        mean of x
--- my:        mean of y
--- cov:       covariance between x and y
--- varx:      variance of x
--- slope:     slope of the regression (returned)
--- intercept: intercept (returned)
+++ Return:    void
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**/
void linreg_coefs(double mx, double my, double cov, double varx, 
                  double *slope, double *intercept){

  linreg_slope(cov, varx, slope);
  linreg_intercept(*slope, mx, my, intercept);

  return;
}


/** Compute correlation coefficient of a linear regression
+++ This function computes the r of a linear regression.
--- cov:    covariance between x and y
--- varx:   variance of x
--- vary:   variance of y
--- r:      correlation coefficient (returned)
+++ Return: void
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**/
void linreg_r(double cov, double varx, double vary, double *r){

  *r = cov/sqrt(varx*vary);

  return;
}


/** Compute coefficient of determination of a linear regression
+++ This function computes the R^2 of a linear regression.
--- cov:    covariance between x and y
--- varx:   variance of x
--- vary:   variance of y
--- rsq:    coefficient of determination (returned)
+++ Return: void
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**/
void linreg_rsquared(double cov, double varx, double vary, double *rsq){

  *rsq = cov*cov/(varx*vary);

  return;
}


/** Predict a value based on a linear regression
+++ This function predicts y values.
--- x:         x value
--- slope:     slope of the regression
--- intercept: intercept
--- y:         predicted y value (returned)
+++ Return:    void
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**/
void linreg_predict(double x, double slope, double intercept, double *y){

  *y = intercept + slope*x;

  return;
}

