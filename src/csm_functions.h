// load Rcpp
#include <Rcpp.h>
using namespace Rcpp;

int get_reaction (double p_lam1, double p_lam2, double p_lam3, double p_mu1, double p_mu2) ;
NumericMatrix speciation_good (NumericMatrix tree, int number_good, double t, int s) ;
NumericMatrix speciation_nasc (NumericMatrix tree, int number_nasc, double t, int s) ;
NumericMatrix completion_nasc (NumericMatrix tree, int number_nasc, double t, int s) ;
NumericMatrix extinction_good (NumericMatrix tree, int number_good, double t, int s) ;
NumericMatrix extinction_nasc (NumericMatrix tree, int number_nasc, double t, int s) ;
NumericMatrix resize_numeric_matrix (NumericMatrix nm, int n_rows, int n_cols) ;
NumericMatrix update_trait (NumericMatrix trait, int trait_vars, NumericMatrix tree, int s, double pp) ;
