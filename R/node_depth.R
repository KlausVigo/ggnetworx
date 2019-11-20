#' These functions return the depths or heights of nodes and tips.
#'
#' @title Depth of Nodes
#' @param x an object of class "evonet"
#' @param \dots Further arguments passed to or from other methods.
#' @return a vector with the depth of the nodes
#' @seealso \code{\link[ape]{node.depth}}
#' @examples
#' z <- ape::read.evonet(text = "((1,((2,(3,(4)Y#H1)g)e,
#' (((Y#H1, 5)h,6)f)X#H2)c)a,((X#H2,7)d,8)b)r;")
#' nd <- node.depth.evonet(z)
#' z$edge.length <- nd[z$edge[,1]] - nd[z$edge[,2]]
#' ggevonet(z)
#'
#' @importFrom phangorn getRoot Ancestors Descendants
#' @export
node.depth.evonet <- function(x, ...){
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
      active <- vapply(candidates, function(x) all(done[pa[[x]]]), FALSE)
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
