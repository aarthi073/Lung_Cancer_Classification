#see what functions and data a package has
babynames::births
#use library() to load it into memory
library(babynames)
#now you have direct access
births
#library loads a package, referring to where the package is contained (like a folder)
#package is a collection/bundle of functions
#lapply() allows you to load multiple packages
packages <- c("dplyr", "ggplot2", "tidyr")
lapply(packages, library, character.only = TRUE)

#unload a package
detach("package:babynames", unload=TRUE)
#help
help(vioplot, package="vioplot")
#see what is inside a package
library(babynames)
ls("package:babynames")
#vignettes are documentation of functionalities of package in a more detailed way
vignette(package="ggplot2")
vignette("ggplot2-specs")
