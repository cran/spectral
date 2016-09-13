#include <Rcpp.h>
using namespace Rcpp;

//#include <omp.h>
// [[Rcpp::plugins(openmp)]]

//' Lomb-Scargle estimation function
//'
//' calculates the Lomb-Scargle estimation
//'
//' This wrapper function calls compiled C++ code, which does the job.
//'
//' @param x spatial vector
//' @param y data vector
//' @param omega frequency vector
// [[Rcpp::export]]

List lmb(NumericVector x, NumericVector y, NumericVector omega)
{
  // declare variables
  double XC,XS,CC,SS,CS,co,si,tauL,ct,st,R,I;
  R_len_t  i,j;

  int n  = x.size();
  int fn = omega.size();
  NumericVector amp(fn);
  NumericVector phase(fn);

//#pragma omp parallel for private(i,j,XC,XS,CC,SS,CS,co,si,tauL,ct,st,R,I)
  for(i = 0; i < fn; i++) {
    XC = 0;
    XS = 0;
    CC = 0;
    SS = 0;
    CS = 0;
    for(j=0; j < n; j++) {
      co = cos(omega[i]*x[j]);
      si = sin(omega[i]*x[j]);
      XC = XC + y[j] * co;
      XS = XS + y[j] * si;
      CC = CC + co*co;
      SS = SS + si*si;
      CS = CS + co*si;
    }
    tauL = atan(2.0 * CS / (CC - SS)) / (2.0*omega[i]);
    ct   = cos(omega[i]*tauL);
    st   = sin(omega[i]*tauL);

    R = ct*XC + st*XS;
    I = ct*XS - st*XC;

    amp[i] = sqrt(2.0 / n * (((R*R)/ (ct*ct * CC + 2.0 * ct * st * CS + st*st * SS)) + ((I*I)/ (ct*ct * SS - 2.0 * ct * st * CS + st*st * CC))));

    phase[i] = -omega[i] * tauL - atan2(I, R);
  }

  return List::create( _["amp"] = amp, _["phase"] = phase);
}


