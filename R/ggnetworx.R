
##' @method fortify evonet
##' @importFrom ggplot2 fortify
##' @export
fortify.evonet <- function(model, data,
                           layout   ="rectangular",
                           ladderize=FALSE,
                           right    =FALSE,
                           mrsd     =NULL,
                           as.Date  =FALSE, ...){
    class(model) <- "phylo"
#  ggtree:::fortify.phylo
    df <- fortify(model, ladderize=ladderize)

    hybridEdge <- logical(nrow(df))
    hybridEdge[grep("#", df$label)] <- TRUE
    df <- cbind(df, hybridEdge=hybridEdge)

    reticulation <- model$reticulation
    df.ret <- df[reticulation[,1], , drop=FALSE]
#    df.ret <- df[reticulation[,2], , drop=FALSE]
    df.ret[,c("node", "parent")] <- reticulation
    df.ret[, "hybridEdge"] <- TRUE
    df <- rbind(df, df.ret)
    df
}


#' drawing phylogenetic tree from phylo object
#'
#'
#' @title ggevonet
#' @param tr phylo object
#' @param mapping aes mapping
#' @param layout one of 'rectangular', 'slanted'
#' @param open.angle open angle, only for 'fan' layout
#' @param mrsd most recent sampling date
#' @param as.Date logical whether using Date class in time tree
#' @param yscale y scale
#' @param yscale_mapping yscale mapping for category variable
#' @param ladderize logical
#' @param right logical
#' @param branch.length variable for scaling branch, if 'none' draw cladogram
#' @param ndigits number of digits to round numerical annotation variable
#' @param ... additional parameter
#' @return tree
#' @seealso \code{\link[ape]{evonet}}, \code{\link[ggtree]{ggtree}}
#' @importFrom magrittr %<>%
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 annotate
#' @importFrom ggplot2 scale_x_reverse
#' @importFrom ggplot2 ylim
#' @importFrom ggplot2 coord_flip
#' @importFrom ggplot2 coord_polar
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 aes_
#' @importFrom ggtree geom_tree2
#' @importFrom ggtree theme_tree
#' @importFrom phangorn coords
#' @importFrom ape read.evonet
#' @author Klaus Schliep
#' @examples
#' (enet <- ape::read.evonet(text="((a:2,(b:1)#H1:1):1,(#H1,c:1):2);"))
#' ggevonet(enet) + geom_tiplab()
#' @export
ggevonet <- function (tr, mapping=NULL, layout="rectangular", open.angle=0,
          mrsd=NULL, as.Date=FALSE, yscale="none", yscale_mapping=NULL,
          ladderize=FALSE, right=FALSE, branch.length="branch.length",
          ndigits=NULL, ...)
{
    layout %<>% match.arg(c("rectangular", "slanted"))
# , "fan", "circular", "radial", "unrooted", "equal_angle", "daylight"
    if (yscale != "none") {
        layout <- "slanted"
    }
    if (is.null(mapping)) {
        mapping <- aes_(~x, ~y)
    }
    else {
        mapping <- modifyList(aes_(~x, ~y), mapping)
    }
    p <- ggplot(tr, mapping=mapping, layout=layout, mrsd=mrsd, as.Date=as.Date,
                yscale=yscale, yscale_mapping=yscale_mapping,
                ladderize=ladderize, right=right, branch.length=branch.length,
                ndigits=ndigits, ...)

    p <- p + geom_tree2(layout=layout, ...)
    p <- p + theme_tree()
    class(p) <- c("ggtree", class(p))
    return(p)
}


##' @method fortify networx
##' @importFrom phangorn getRoot
##' @export
fortify.networx <- function(model, data,
                           layout   ="unrooted",
                           ladderize=FALSE,
                           right    =FALSE,
                           mrsd     =NULL,
                           as.Date  =FALSE, ...){
#    root <- getRoot(model)
    nTips <- length(model$tip.label)
    label <- character(nrow(model$edge))
    isTip <- logical(nrow(model$edge))  # edge leading to tip
    # 1:nTips may not correspond to tips
    if(!is.null(model$translate)){
        ind <- match(model$translate$node, model$edge[,2])
        label[ind] <- model$translate$label
    }
    else{
        ind <- match(1:nTips, model$edge[,2])
        label[ind] <- model$tip.label
    }
    isTip[ind] <- TRUE
    df <- data.frame(node=model$edge[,2], parent=model$edge[,1],
                     branch.length=model$edge.length,
                     split=model$splitIndex,
                     label=label, isTip=isTip)
    if(!is.null(model$.plot)) coord <- model$.plot$vertices
    else coord <- coords(model, dim="2D")
    df <- cbind(df, x=coord[df$node,1], y=coord[df$node,2],
                xend=coord[df$parent,1], yend=coord[df$parent,2])
    angle <- atan2(df$y - df$yend, df$x - df$xend) * 360 / (2*pi)
    angle[angle<0] <- angle[angle<0] + 360
    df <- cbind(df, angle=angle)
#    rownames(df) <- 1:nrow(df)
    df
}


#' drawing phylogenetic tree from phylo object
#'
#'
#' @title ggsplitnet
#' @param tr phylo object
#' @param mapping aes mapping
#' @param layout one of 'rectangular', 'slanted', 'fan', 'circular', 'radial',
#' 'equal_angle' or 'daylight'
#' @param open.angle open angle, only for 'fan' layout
#' @param mrsd most recent sampling date
#' @param as.Date logical whether using Date class in time tree
#' @param yscale y scale
#' @param yscale_mapping yscale mapping for category variable
#' @param ladderize logical
#' @param right logical
#' @param branch.length variable for scaling branch, if 'none' draw cladogram
#' @param ndigits number of digits to round numerical annotation variable
#' @param ... additional parameter
#' @return tree
#' @seealso \code{\link[ggtree]{ggtree}}, \code{\link[phangorn]{networx}},
#' \code{\link[phangorn]{consensusNet}}, \code{\link[phangorn]{neighborNet}}
#' @references Schliep, K., Potts, A. J., Morrison, D. A. and Grimm, G. W.
#' (2017), Intertwining phylogenetic trees and networks.
#' \emph{Methods Ecol Evol}. \bold{8}, 1212--1220. doi:10.1111/2041-210X.12760
#' @importFrom utils modifyList
#' @importFrom magrittr %<>%
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 annotate
#' @importFrom ggplot2 scale_x_reverse
#' @importFrom ggplot2 ylim
#' @importFrom ggplot2 coord_flip
#' @importFrom ggplot2 coord_polar
#' @importFrom ggtree theme_tree
#' @author Klaus Schliep
#' @examples
#' data(yeast, package="phangorn")
#' dm <- phangorn::dist.ml(yeast)
#' nnet <- phangorn::neighborNet(dm)
#' ggsplitnet(nnet) + geom_tiplab2()
#'
#'
#' @export
ggsplitnet <- function (tr, mapping=NULL, layout="slanted", open.angle=0,
        mrsd=NULL, as.Date=FALSE, yscale="none", yscale_mapping=NULL,
        ladderize=FALSE, right=FALSE, branch.length="branch.length",
        ndigits=NULL, ...)
{
    layout %<>% match.arg(c("slanted"))
# "rectangular", "fan", "circular", "radial", "unrooted", "equal_angle",
#    "daylight"
    if (is.null(mapping)) {
        mapping <- aes_(~x, ~y)
    }
    else {
        mapping <- modifyList(aes_(~x, ~y), mapping)
    }
    p <- ggplot(tr, mapping=mapping, layout=layout, mrsd=mrsd, as.Date=as.Date,
                yscale=yscale, yscale_mapping=yscale_mapping,
                ladderize=ladderize, right=right, branch.length=branch.length,
                ndigits=ndigits, ...)

    p <- p + geom_tree3(layout=layout, ...)
    p <- p + theme_tree()
    class(p) <- c("ggtree", class(p))
    return(p)
}





##' add tree layer
##'
##'
##' @title geom_tree3
##' @param layout one of 'rectangular', 'slanted', 'circular', 'radial' or
##' 'unrooted'
##' @param ... additional parameter
##' @return tree layer
##' @importFrom ggplot2 geom_segment
##' @importFrom ggplot2 aes
##' @export
##' @author Klaus Schliep
geom_tree3 <- function(layout="rectangular", ...) {
    x <- y <- xend <- yend <- parent <- NULL
    lineend <- "round"
    if (layout == "rectangular" || layout == "fan" || layout == "circular") {
        list(
            geom_segment(aes(x    = x,
                             xend = xend,
                             y    = y,
                             yend = y),
                         lineend  = lineend, ...),

            geom_segment(aes(x    = xend,
                             xend = xend,
                             y    = y,
                             yend = yend),
                         lineend  = lineend, ...)
        )
    } else if (layout == "slanted" || layout == "radial" ||
               layout == "unrooted") {
        geom_segment(aes(x    = x,
                         xend = xend,
                         y    = y,
                         yend = yend),
                     lineend  = lineend, ...)
    }
}



##' pipe
##' @importFrom magrittr %>%
##' @name %>%
##' @export
##' @rdname pipe
##' @param lhs left hand side
##' @param rhs right hand side
##' @usage lhs \%>\% rhs
##' @seealso
##' \link[magrittr]{pipe}
NULL
