#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace std;
using namespace Rcpp;
#include<iostream>

// [[Rcpp::export]]
double nll_Rcpp(arma::mat allprobs, arma::mat gamma, arma::rowvec delta, int T_y)
{
	double lscale;
	int t=1;
	double sumfoo;
	arma::rowvec foo;
	foo=delta % allprobs.row(0);
	sumfoo=sum(foo);
	lscale=log(sumfoo);
	foo=foo/sumfoo;
	for(t=1;t<T_y;t++)
	{		
		foo=(foo*gamma)%allprobs.row(t);
		sumfoo=sum(foo);
		lscale=lscale+log(sumfoo);
		foo=foo/sumfoo;
	}
	return -lscale;
}
