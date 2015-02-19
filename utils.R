load_dependencies <- function() {
  dependencies = c("skmeans", "ggplot2","pls", "maps", "maptools", "RColorBrewer", 
                   "classInt", "gpclib", "cluster", "rgeos", "PBSmapping", "scales",
                   "plyr", "spdep", "xtable", "pracma","grid", "dplyr")
  installed_packages = rownames(installed.packages())
  lapply(dependencies, function(x){
    if (!(x %in% installed_packages)){
      install.packages(x)
    }
    library(x, character.only=TRUE)
  })
}


gg_color_n <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[n]
}


toDF = function(json) {
  r = lapply(json, function(x) {
    x[sapply(x, is.null)] <- NA
    unlist(x)
  })
  r = do.call("rbind", r)
  r = as.data.frame(r)
}

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}