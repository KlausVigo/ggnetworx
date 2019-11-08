library(phangorn)
x <- allCircularSplits(6)
z <- as.networx(x)
plot(z, "equal angle", show.tip.label = FALSE, edge.color = "grey")
