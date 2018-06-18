

get_next_oldest <- function(current_node, taxon, tree) {

   ytree <- matrix(tree[ -(1:current_node), ],ncol=6)
   children_of_taxon <- which( ytree[,3] == taxon)
   if (length( children_of_taxon) == 0 ) {
      next_node  <- taxon
      next_time  <- 1
   } else {
      i_next <- children_of_taxon[1]
      next_node  <- ytree[i_next,6]
      next_time  <- ytree[i_next,2]
   }

   next_oldest = list(node=next_node, time=next_time)
   return(next_oldest)
}


csm2phylo <- function( flat_tree ) {

   Ntips     <- nrow(flat_tree) 
   Nnode     <- Ntips - 1
   Nedges    <- Ntips + Nnode -1
   tip.label <- as.character(flat_tree[,1])

   #cat(Ntips,"\n")
   tree <- cbind(flat_tree[-1,], (Nnode+(1:Nnode)))

   colnames(tree) <- c("child", "time", "parent", "status", "time_ext", "new_node")

   edge         <- mat.or.vec(Nedges,2)
   edge.length  <- mat.or.vec(Nedges,1)

   c <- 1
   for (i in 1:nrow(tree)) {

      focal_node     <- tree[i, "new_node"]
      child_at_node  <- tree[i, "child"]
      parent_at_node <- tree[i, "parent"]

      i_child_next  <- get_next_oldest(i, child_at_node, tree)
      i_parent_next <- get_next_oldest(i, parent_at_node,  tree)

      edge[c,1]      <- i_child_next$node
      edge[c,2]      <- focal_node
      edge.length[c] <- i_child_next$time - tree[i, "time"] 
      c <- c+1

      edge[c,1] <- i_parent_next$node
      edge[c,2] <- focal_node
      edge.length[c] <- i_parent_next$time - tree[i, "time"] 
      c <- c+1

   }

   edge  <- matrix(sapply(edge,as.integer),ncol=2)
   Nnode <- as.integer(Nnode)

   edge=edge+1L
   edge=edge[,2:1]
   atre <- list(edge=edge, tip.label=tip.label, edge.length=edge.length, Nnode=Nnode) 

   class(atre) <- "phylo"
   atre <- as.phylo(atre)

   return(atre)
}

