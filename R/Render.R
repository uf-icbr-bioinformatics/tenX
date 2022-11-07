library("rmarkdown")

args = commandArgs(trailingOnly=TRUE)
rmarkdown::render(args[1])
