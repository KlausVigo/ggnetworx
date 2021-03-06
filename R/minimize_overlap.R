#' @title minimize_overlap
#' reduces reticulation lines crossing over in plots
#' @param x Tree of class 'evonet'
#' @return  A Tree with rotated nodes of class 'evonet'
#' @importFrom ape node.height rotate Ntip
#' @importFrom phangorn Ancestors
#' @author L. Francisco Henao Diaz
#' @examples
#' fishnet <- ape::read.evonet(text="(Xalvarezi,Xmayae,((Xsignum,((Xmonticolus,(Xclemenciae_F2,#H25:::0.032485577051040423):0.979116879342896):1.4858786525152954,(((((((((Xgordoni,Xmeyeri):0.26111415598145865,Xcouchianus):3.5325554148283564,Xvariatus):0.6445452395880061,Xevelynae):0.41219807813960146,(Xxiphidium,#H24:0.0::0.16538333427567822):1.4020042304402076):0.2921320715527958,Xmilleri):0.4688709909949548,Xandersi):0.6492623942884177,Xmaculatus):1.0200542801748702,(((Xmontezumae,(Xcortezi,(Xbirchmanni_GARC,Xmalinche_CHIC2):0.8943004029517725):0.6903622798656042):0.4926275928330048,((Xnigrensis,Xmultilineatus):1.4508210422790195,(Xpygmaeus,Xcontinens):1.9686498919562234):2.428065080719157):2.0435103751744865)#H24:0.6476450483433508::0.8346166657243218):0.7816554661308311):1.9841381634366537):0.31614675520900004,(Xhellerii)#H25:::0.9675144229489596):0.28226883225049243);")
#' fishnet$edge.length <- NULL
#' new_tre <- minimize_overlap(fishnet)
#'
#' par(mfrow=c(1,2))
#' ggevonet(fishnet)
#' ggevonet(new_tre)
#'
#' net2 <- ape::read.evonet(text="(15,(1,((14,(#H1,(((12,13),(11,#H3)),(7,
#'    ((10)#H3,(8,9)))))),((((2,3))#H2,(6,(5,(#H2,4)))))#H1)));")
#' # Cui et al. 2013 Evol.
#' new_net2 <- minimize_overlap(net2)
#' ggevonet(net2)
#' ggevonet(new_net2)
#'
#' @export
minimize_overlap <- function(x){
    if(class(x)[1]!='evonet') stop("x should be an 'evonet' class")
    n_iter <- round(x$Nnode*3/4)
#    r_hist <- c()
    for(j in seq_len(n_iter)){
        h <- node.height(x)
        best_r <- sum(abs(h[x$reticulation[,1]]- h[x$reticulation[,2]]))
        best_c <- -1
#        r_hist[j] <- best_r
        nodes2rot <- intersect(sort(unique(unlist(Ancestors(x,
                        c(x$reticulation))))), which(tabulate(x$edge[,1]) > 1) )
        for(i in seq_along(nodes2rot)){
            nh <- node.height(ape::rotate(x, nodes2rot[i]))
            best_nr <- sum(abs(nh[x$reticulation[,1]] - nh[x$reticulation[,2]]))
            if(best_nr < best_r){
                best_c <- nodes2rot[i]
                best_r <- best_nr
            }
        }
        if(best_c > 0) x <- ape::rotate(x, best_c)
        else break()
    }
    x
}

