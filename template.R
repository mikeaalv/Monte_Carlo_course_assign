rm(list=ls())
options(warn=1)
options(stringsAsFactors=FALSE)
options(digits=15)
require(stringr)
require(magrittr)
require(R.matlab)
require(ggplot2)
require(xml2)
require(cowplot)
dir="/Users/mikeaalv/Library/Mobile Documents/com~apple~CloudDocs/working/course/phys8601/assignment/assign2/result/"
setwd(dir)
##plot infinite cluster probability with error bar
