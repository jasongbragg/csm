// load Rcpp
#include <Rcpp.h>
#include "csm_functions.h"

using namespace Rcpp;

// [[Rcpp::export]]
List csm(double lam1, double lam2, double lam3, double mu1, double mu2, int crown_or_stem, double pp, double t_step) {


    // double pp = 0.01;
    // int crown_or_stem = 0;

    // max time, all normalized to 1
    double t_max = 1.0;

    // max number of Gillespie increments
    int ginc_max = 100000;

    // number of data columns
    int tree_vars   = 5;
    int trace_vars  = 7;
    int trait_vars  = 5;

    // make arrays with row numbers sufficient 
    // for all gillespie increments
    NumericMatrix tree(ginc_max, tree_vars);
    NumericMatrix trace(ginc_max, trace_vars);
    NumericMatrix trait(ginc_max, trait_vars);

    // Gillespie step random
    double v_rand;
 
    // Gillespie time step
    double t_ginc;
    //double t_step;
    //t_step = 0.01;

    // initialize time
    // note: first speciation event will occur 
    // at this time. Will be reaction = 1.
    double t=0;
    double tt=0;

    // initialize species index
    int s = 1;

    // initialize tree
    // zeroeth species born good at time -1
    tree(0,0) = 0;
    tree(0,1) = -1;
    tree(0,2) = -1;
    tree(0,3) = 2;    
    tree(0,4) = 0;    

    // initialize living species counters
    int number_good = 1;
    int number_nasc = 0;
    int number_tot = number_good + number_nasc;

    // declare Gillespie reaction probs
    double p_lam1;
    double p_lam2;
    double p_lam3;
    double p_mu1;
    double p_mu2;
    double p_tot;

    int reaction;
    CharacterVector terminate_status(1);

    // initialize count of Gillespie increments
    int n_ginc = 0;

    // main loop for speciation
    while (t < t_max && number_tot > 0 && n_ginc < ginc_max) {

        // calculate reaction probs
        // pl1 good spec | pl2 nasc spec | pl3 complete nasc | mu1 ext good | mu2 ext nasc 
        // Rprintf(" Gillespie step   : %f\n", t);

        p_lam1 = lam1 * (double) number_good;
        p_lam2 = lam2 * (double) number_nasc;
        p_lam3 = lam3 * (double) number_nasc;
        p_mu1  = mu1  * (double) number_good;
        p_mu2  = mu2  * (double) number_nasc;
        p_tot  = p_lam1 + p_lam2 + p_lam3 + p_mu1 + p_mu2;
        

        // next Gillespie time step
        n_ginc = n_ginc + 1;
        v_rand = R::runif(0,1);
        t_ginc = -1.0 * std::log(v_rand) / p_tot;
        //Rprintf(" Gillespie step: %f\n", t_ginc);

        if (crown_or_stem == 1 && n_ginc == 1 ) {
           t      = t + t_ginc;
        }


        // implement reaction
        reaction = get_reaction(p_lam1, p_lam2, p_lam3, p_mu1, p_mu2);

        if (reaction == 1) { tree = speciation_good( tree, number_good, t, s );  number_nasc++; s++;}
        if (reaction == 2) { tree = speciation_nasc( tree, number_nasc, t, s );  number_nasc++; s++;}
        if (reaction == 3) { tree = completion_nasc( tree, number_nasc, t, s );  number_good++; number_nasc--;}
        if (reaction == 4) { tree = extinction_good( tree, number_good, t, s );  number_good--;}
        if (reaction == 5) { tree = extinction_nasc( tree, number_nasc, t, s );  number_nasc--;}

        number_tot = number_good + number_nasc;

        // next Gillespie time step
        //n_ginc = n_ginc + 1;
        //v_rand = R::runif(0,1);
        //t_ginc = -1.0 * std::log(v_rand) / p_tot;
        //Rprintf(" Gillespie step: %f\n", t_ginc);

        //if (crown_or_stem == 1 && n_ginc == 1 ) {
        //   t      = t + t_ginc;
        //}


        // newly made species traits
        if (reaction == 1 || reaction == 2) {
           int parent  = tree((s-1), 2);
           //Rprintf(" s id   : %i\n", s);  
           //Rprintf(" number   : %i\n", number_tot);  
           //Rprintf(" parent id   : %i\n", parent);  

           for (int i=0; i<trait_vars; ++i) {
              trait((s-1), i) = trait(parent, i) ;
           }
        }

        // updates to traits
        if ( tt <= (t + t_ginc) && tt <= t_max ) {
           
           while ( tt <= (t + t_ginc) && tt <= t_max ) {
              trait = update_trait(trait, trait_vars, tree, s, pp);
              //Rprintf(" uniform step   : %f\n", tt); 
              tt = tt + t_step;
           } 
        }

        // keep a trace -- useful for debug
        trace(n_ginc,0) = reaction;
        trace(n_ginc,1) = number_good;
        trace(n_ginc,2) = number_nasc;
        trace(n_ginc,3) = number_tot;
        trace(n_ginc,4) = t;
        trace(n_ginc,5) = tt;
        trace(n_ginc,6) = t_ginc;
        // end trace

        // update time, with Gillespie step
        t      = t + t_ginc;

    } // end main loop for speciation

    // make some nicer objects to return
    NumericMatrix rtrace = resize_numeric_matrix(trace, (n_ginc+1), trace_vars); 
    NumericMatrix rtree  = resize_numeric_matrix(tree, s, tree_vars); 
    NumericMatrix rtrait = resize_numeric_matrix(trait, s, trait_vars); 

    if ( t >= t_max )           { terminate_status(0) = "timeout"; }
    if ( number_tot <= 0 )      { terminate_status(0) = "extinction";}
    if ( n_ginc >= ginc_max )   { terminate_status(0) = "overstep"; }


    // bundle stuff into a list and return
    List result =
        List::create(terminate_status,
                     rtree, rtrace, rtrait);

    return result;
  
}


