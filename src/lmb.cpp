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
  double xc,xs,cc,ss,cs,co,si,tauL,ct,st,re,im;
  R_len_t  i,j;

  int n  = x.size();
  int fn = omega.size();
  NumericVector amp(fn);
  NumericVector phase(fn);

//#pragma omp parallel for private(i,j,xc,xs,cc,ss,cs,co,si,tauL,ct,st,re,im)
  for(i = 0; i < fn; i++) {
    xc = 0;
    xs = 0;
    cc = 0;
    ss = 0;
    cs = 0;
    for(j=0; j < n; j++) {
      co = cos(omega[i]*x[j]);
      si = sin(omega[i]*x[j]);
      xc = xc + y[j] * co;
      xs = xs + y[j] * si;
      cc = cc + co*co;
      ss = ss + si*si;
      cs = cs + co*si;
    }
    tauL = atan(2.0 * cs / (cc - ss)) / (2.0*omega[i]);
    ct   = cos(omega[i]*tauL);
    st   = sin(omega[i]*tauL);

    re = ct*xc + st*xs;
    im = ct*xs - st*xc;

    amp[i] = sqrt(2.0 / n * (((re*re)/ (ct*ct * cc + 2.0 * ct * st * cs + st*st * ss)) + ((im*im)/ (ct*ct * ss - 2.0 * ct * st * cs + st*st * cc))));

    phase[i] = -omega[i] * tauL - atan2(im, re);
  }

  return List::create( _["amp"] = amp, _["phase"] = phase);
}


