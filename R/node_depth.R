#' @importFrom phangorn getRoot Ancestors Descendants
node_depth_evonet <- function(x, ...){
   x <- ape::reorder.phylo(x)
   root <- getRoot(x)
   max_nodes <- max(x$edge)
   nTip <- Ntip(x)
   desc <- Descendants(x, seq_len(max_nodes), "children")
   anc <- Ancestors(x)

   pa <- vector("list", max_nodes)
   ind <- which(x$edge[,2] > Ntip(x))
   pa[x$edge[ind,2]] <- x$edge[ind,1]
   for(i in seq_len(nrow(x$reticulation))){
      pa[[x$reticulation[i,2]]] <- sort( c(pa[[x$reticulation[i,2]]],
                                           x$reticulation[i,1] ) )
#       pa[[x$reticulation[i,1] ]] <- numeric(0)
   }
   ind <- which(lengths(pa) > 0)
   depth <- numeric(max_nodes)
   depth[root] <- 1
   done <- logical(max_nodes)
   done[root] <- TRUE
   candidates <- desc[[root]]
   candidates <- candidates[candidates>nTip]
   d <- 1

   while(length(candidates)>0){
      active <- sapply(candidates, function(x) all(done[pa[[x]]]))
      tmp <- which(active)[1]  #sample(active,1)

      candidates <- c(candidates, desc[[candidates[tmp] ]])
      candidates <- candidates[candidates>nTip]
      d <- d+1
      done[candidates[tmp]] <- TRUE
      depth[candidates[tmp]] <- d
      candidates <- candidates[-tmp]
   }
   depth <- d+2 - depth
   depth[seq_len(nTip)] <- 1
   depth
}
