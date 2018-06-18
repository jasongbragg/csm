#' Performs a simulation of a protracted birth-death model
#' with phenotypic evolution suitable for identifying 
#' species that are 'cryptic'
#'
#' (Based on models in: )
#'
#' @param lam1 - speciation of good species [Required]
#' @param lam2 - speciation of nascent species [Required]
#' @param lam3 - rate of speciation completion [Required]
#' @param mu1  - extinction of good species [Required]
#' @param mu2  - extinction of nascent species [Required]
#' @param crown_or_stem  - extinction of nascent species 
#' @param pp             - if nascent, p(shared parent phenotype change) 
#' @return an object 
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export
#' @examples
#' \dontrun{
#' simulation <- CSMsim(5,1.5,1.5,0.4,0.4,1,0.01)
#' }

CSMsim <- function(lam1, lam2, lam3, mu1, mu2, crown_or_stem=0, pp=1.0, trait_step=0.01, out="phylo") {

   require(ape)

   csm_out <- csm(lam1, lam2, lam3, mu1, mu2, crown_or_stem=crown_or_stem, pp=pp, t_step=trait_step)

   run_status <- csm_out[[1]]

   if (run_status == "extinction") {
      #cat("extinction \n")
      return(list(run_status=run_status))
   } 

   if (run_status == "overstep") {
      #cat("overstepping \n")
      return(list(run_status=run_status))
   } 

   if (run_status == "timeout" & out == "phylo") {

      # get the phylogeny in ape phylo class
      tax_tree   <- csm_out[[2]]
      i_extant   <- which(tax_tree[,4] > 0 )

      if (length(i_extant) <= 2) {
         return(list(run_status="TwoOrFewerExtant"))
      } else {

         csm_phy    <- csm2phylo(tax_tree)

         # get traits
         tax_traits <- csm_out[[4]]

         # remove extinct
         i_extinct  <- which(tax_tree[,4] == 0 )
         if (length(i_extinct) == 0) {
            # none extinct
         } else {
            tips_to_drop <- csm_phy$tip.label[i_extinct] 
            csm_phy      <- drop.tip(csm_phy, tips_to_drop)
            tax_traits   <- tax_traits[ -i_extinct, ]
         }

         
         rownames(tax_traits) <- csm_phy$tip.label

         csm_sim <- list(run_status=run_status, csm_phy=csm_phy, csm_traits=tax_traits)

         return(csm_sim)
      }
   }

   if (run_status == "timeout" & out == "flat") {
      return(csm_out)
   }

}
