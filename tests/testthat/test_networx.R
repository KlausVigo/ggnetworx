context("networx")


data(yeast, package="phangorn")
dm <- phangorn::dist.ml(yeast)
nnet <- phangorn::neighborNet(dm)

(enet <- ape::read.evonet(text="((a:2,(b:1)#H1:1):1,(#H1,c:1):2);"))

splits_network <-  ggsplitnet(nnet) + geom_tiplab2()
reticulate_network <- ggevonet(enet) + geom_tiplab()


# Visual tests ------------------------------------------------------------
#test_that("visual appearance", {
#    testthat::skip_on_cran()

#    vdiffr::expect_doppelganger("Basic splits network", splits_network)
#    vdiffr::expect_doppelganger("Basic reticulate network", reticulate_network)
#})
