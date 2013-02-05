#include <RcppGSL.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppGSL)]]
// [[Rcpp::export]]
NumericMatrix RcppGibbs(int N, int thn) {
  int i,j;
  NumericMatrix mat(N, 2);
  RNGScope scope;         // Initialize Random number generator

  // The rest of the code follows the R version
  double x=0, y=0;

  for (i=0; i<N; i++) {
    for (j=0; j<thn; j++) {
      x = ::Rf_rgamma(3.0,1.0/(y*y+4));
      y = ::Rf_rnorm(1.0/(x+1),1.0/sqrt(2*x+2));
    }
    mat(i,0) = x;
    mat(i,1) = y;
  }
  return mat;             // Return to R
}

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
// [[Rcpp::export]]
NumericMatrix GSLGibbs(int N, int thin) {
  int i, j;
  gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
  double x=0, y=0;
  NumericMatrix mat(N, 2);
  for (i=0; i<N; i++) {
    for (j=0; j<thin; j++) {
      x = gsl_ran_gamma(r,3.0,1.0/(y*y+4));
      y = 1.0/(x+1)+gsl_ran_gaussian(r,1.0/sqrt(2*x+2));
    }
    mat(i,0) = x;
    mat(i,1) = y;
  }
  gsl_rng_free(r);

  return mat;           // Return to R
}

#include <boost/random.hpp>
#include <boost/generator_iterator.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/gamma_distribution.hpp>

typedef boost::mt19937 RNGType;     // select a generator, MT good default
RNGType rng(123456);			// instantiate and seed

boost::normal_distribution<> n01(0.0, 1.0);
boost::variate_generator< RNGType, boost::normal_distribution<> > rngNormal(rng, n01);

boost::gamma_distribution<> g3(3.0); // older Boost took one arg "alpha", divide draw by "beta"
boost::variate_generator< RNGType, boost::gamma_distribution<> > rngGamma(rng, g3);
// [[Rcpp::export]]
NumericMatrix BoostGibbs(int N, int thin) {
  int i, j;
  double x=0, y=0;
  NumericMatrix mat(N, 2);

  for (i=0; i<N; i++) {
    for (j=0; j<thin; j++) {
      x = rngGamma()/(1.0/(y*y+4));     // dividing by beta gives us Gamma(alpha, beta)
      y = 1.0/(x+1) + rngNormal()*(1.0/sqrt(2*x+2));  // scale by sigma and move by mu
    }
    mat(i,0) = x;
    mat(i,1) = y;
  }
  return mat;           // Return to R
}
