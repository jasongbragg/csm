#' Performs optimization of morphospace threshold value
#'
#' @param target_model   - a csm model object [Required]
#' @param reference_nCAD - normalized cryptic ancestor depths vector [Required]
#' @param start_mt_val   - starting value for optimization [Required]
#' @return a list with optim convergence status (0 is good), mt, and value of nCAD diff 
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export
#' @examples
#' \dontrun{
#' morph_thresh_optim(test_csm, obs_nCAD)
#' }


fit_cryptic_proportion <- function(target_model, reference_model, start_mt_val=1) {

    rm_summ <- CSMsumm(reference_model,cryptic=TRUE,tinc=0.0005)

    get_cryptic_proportion <- function(model_summ) {
       return( model_summ$final_cryptic_pairs / (model_summ$final_lineages * (model_summ$final_lineages-1)) )
    }

    ref_cryptic_proportion <- get_cryptic_proportion(rm_summ)

    optim_function <- function(x) {

       target_model  <- CSMcryptic(target_model, morph_thresh=x)

       if (is.list(target_model$cryptic)) {
          tm_summ   <- CSMsumm(target_model,cryptic=TRUE,tinc=0.0005)
          target_cryptic_proportion  <- get_cryptic_proportion(tm_summ)
       } else {
          target_cryptic_proportion  <- 0
       } 
       cryp_prop_diff     <- sum(abs((target_cryptic_proportion - ref_cryptic_proportion)))
       return(cryp_prop_diff)

    }

    opt_out <- optim(start_mt_val, optim_function,method="Brent",lower=0,upper=20)
    opt_morpho_thresh <- opt_out$par
    min_prop_diff     <- opt_out$value
    status            <- opt_out$convergence

    return(list(status=status, opt_morpho_thresh=opt_morpho_thresh, min_prop_diff=min_prop_diff))
}
