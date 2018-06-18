#' Apply thresholds for crypsis, identify cryptic lineages in CSM simulation
#'
#' @param csm_sim - speciation of good species [Required]
#' @param csm_sim - morph_thresh
#' @return an object 
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export
#' @examples
#' \dontrun{
#' csm_sim <- CSMcryptic(csm_sim, morph_thresh=0.05)
#' }

CSMcryptic <- function(csm_sim, morph_thresh, monophyly=FALSE) {

   run_status <- csm_sim$run_status

   if (run_status != "timeout") {
      csm_sim <- csm_sim
   } else {

      csm_phylo   <- csm_sim$csm_phy
      tax_traits  <- csm_sim$csm_traits
      trait_dist  <- dist(tax_traits)

      cryptic_pairs <- which( ( as.matrix( trait_dist )) < morph_thresh, arr.ind=TRUE)
      cryptic_pairs <- cryptic_pairs[ which(cryptic_pairs[,1] < cryptic_pairs[,2]), ]


      if (length(cryptic_pairs) == 0) {
         csm_sim[[ "cryptic" ]] <- "None"
      } else {

         if (length(cryptic_pairs) == 2) { 
            cryptic_pairs <- matrix(cryptic_pairs,ncol=2)
         }

         cryptic_taxa  <- cbind( rownames(as.matrix(trait_dist))[ cryptic_pairs[,1] ], rownames(as.matrix(trait_dist))[ cryptic_pairs[,2] ]  )

         count_cryptic <- table(c(cryptic_taxa[,1], cryptic_taxa[,2]))


         anc_depths    <- cophenetic.phylo(csm_phylo)/2

         cryp_depths   <- mat.or.vec(nrow(cryptic_taxa),1)

         for (i in 1:nrow(cryptic_taxa)) {
            cryp_depths[i] <- anc_depths[cryptic_taxa[i,1],cryptic_taxa[i,2]] 
         }

         if (monophyly) {

            cryp_status    <- mat.or.vec(nrow(cryptic_taxa),3)

            for (i in 1:nrow(cryptic_taxa)) {

               ca       <- getMRCA(phy=csm_phylo, c(cryptic_taxa[i,1],cryptic_taxa[i,2])) 
               ca_tree  <- extract.clade(csm_phylo, ca)
               cai1 <- 0
               cai2 <- 0

               if (length(ca_tree$tip.label) == 2) {
                  # do nothing, status = 0, 0
               } else {
                  
                  cryp_status[i,1] <- length(ca_tree$tip.label) - 2

                  # in clade, how many cryptic with each
                  cl1 <- cryptic_taxa[i,1]
                  cl2 <- cryptic_taxa[i,2]

                  icl1   <- which( cryptic_taxa[,1] == cl1 | cryptic_taxa[,2] == cl1)
                  cwith1 <- unique(c(cryptic_taxa[icl1,1], cryptic_taxa[icl1,2])) 
                  cint1  <- intersect(cwith1, ca_tree$tip.label)                  
                  cai1   <- length(cint1) -2

                  icl2    <- which( cryptic_taxa[,1] == cl2 | cryptic_taxa[,2] == cl2)
                  cwith2 <- unique(c(cryptic_taxa[icl2,1], cryptic_taxa[icl2,2]))
                  cint2  <- intersect(cwith2, ca_tree$tip.label)
                  cai2   <- length(cint2) -2 
               }
               cryp_status[i,2:3] <- c(cai1, cai2)
            }
         }

         cryptic_lineages <- list(pairs=cryptic_taxa, depths=cryp_depths, counts=count_cryptic)
         if (monophyly) {
            cryptic_lineages[[ "monophyly" ]] <- cryp_status
         }
         csm_sim[[ "cryptic" ]] <- cryptic_lineages
      } 
   }

   return(csm_sim)
}
