#' Normalize CSM simulation to have crown age equal to 1
#'
#' @param sim - speciation of good species [Required]
#' @return an object 
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export
#' @examples
#' \dontrun{
#' csm_sim <- CSMcryptic(csm_sim, morph_thresh=0.05)
#' }

CSMnorm <- function(csm_sim) {

   run_status <- csm_sim$run_status

   if (run_status != "timeout") {
      csm_out <- csm_sim
   } else {

      if ( is.null(csm_sim$cryptic) ) {

         csm_new             <- csm_sim
         mbt                 <- max(branching.times(csm_sim$csm_phy))
         csm_new$csm_phy$edge.length <- csm_sim$csm_phy$edge.length / mbt
         csm_new[[ "norm_crown_1" ]] <- TRUE
         csm_out <- csm_new

      } else {
         cat(" Normalization not performed if cryptic lineages in simulation object \n")
         csm_out <- csm_sim
      }

   }
   return(csm_out)
}
