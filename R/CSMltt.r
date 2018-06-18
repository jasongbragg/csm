#' Makes a lineage through time plot using output
#' from a CSM simulation 
#'
#' @param csm_sim - speciation of good species [Required]
#' @param cryptic - speciation of nascent species [Required]
#' @return an object 
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export
#' @examples
#' \dontrun{
#' plots <- CSMltt(csm_sim, cryptic=TRUE)
#' }

CSMltt <- function(csm_sim, cryptic=FALSE, as_pdf=FALSE, pdf_file=NULL) {

   phy <- csm_sim$csm_phy
   pch <- 19 
   cex <- 1.3
   lin_col <- "blue"        

   brtimes <- rev( sort( as.numeric(branching.times(phy)) ) )
   brtimes <- brtimes/max(brtimes)
   lins    <- log2(2:(length(brtimes)+1))                        # log(2:(length(bt)+1));
   bt      <- max(brtimes)-brtimes                               # convert to time t=0 at root

   if (as_pdf) {
      pdf(file=pdf_file)
   } 
   #plot.new()
   par(mfrow=c(2,1), mar=c(6,6,2,2))
 
   crtimes <- sort(as.numeric(1-csm_sim$cryptic$depths))
   plot(crtimes,log2(1:length(crtimes)), xlim=c(0,1), pch=pch, cex=cex, col="blue", xlab="",ylab="")
   mtext(side=1, text="time", line=3.5, cex=1.2)
   mtext(side=2, text="cryptic pairs [log2]", line=3.5, cex=1.2)

   ml <- max(lins)
   yaxis <- ceiling(ml * 1.05 / 2 )*2

   plot.new()
   plot.window(xlim=c(0,1), ylim=c(0.5, yaxis))
   points(bt, lins, pch=pch, cex=cex, col=lin_col)		
 
   yinc <- 1
   axis(1, at=seq(-0.2, 1.2, by=0.2))
   axis(2, at=seq(0, yaxis, by=yinc))
   mtext(side=1, text="time", line=3.5, cex=1.2)
   mtext(side=2, text="lineages [log2]", line=3.5, cex=1.2)

   if (as_pdf) {
      dev.off()
   } 
   
}
