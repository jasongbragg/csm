#' Performs optimization of morphospace threshold value
#'
#' @param target_model   - a csm model object [Required]
#' @param reference_nCAD - normalized cryptic ancestor depths vector [Required]
#' @param start_mt_val   - starting value for optimization [Required]
#' @return a list with optim convergence status (0 is good), mt, and value of nCAD diff 
#' @author Jason Bragg (jasongbragg@gmail.com)


 morph_thresh_optim <- function(target_model, reference_nCAD, start_mt_val=1) {

    optim_function <- function(x) {

       target_model  <- CSMcryptic(target_model, morph_thresh=x)

       if (is.list(target_model$cryptic)) {
          target_summ   <- CSMsumm(target_model,cryptic=TRUE,tinc=0.0005)
          target_nCAD   <- target_summ$norm_cryptic_depths
       } else {
          target_nCAD <- rep(0,(1/0.0005)+1)
       } 
       nCAD_diff     <- sum(abs((target_nCAD-reference_nCAD)))
       return(nCAD_diff)

    }

    opt_out <- optim(start_mt_val, optim_function,method="Brent",lower=0,upper=20)
    opt_morpho_thresh <- opt_out$par
    min_nCAD_diff     <- opt_out$value
    status            <- opt_out$convergence

    return(list(status=status, opt_morpho_thresh=opt_morpho_thresh, min_nCAD_diff=min_nCAD_diff))
 }
