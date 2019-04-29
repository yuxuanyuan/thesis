#!/usr/bin/env Rscript

############################################################
# Written by Andy Yuan                                     #
# Email: yuxuan.yuan@outlook.com                           #
# Last modified date: 14/03/2018                           #
############################################################

##====================== description ======================#
#  calculate a genome size based on the kmer histo table   #
#    the histo table can be obtained using jellyfish       #
#       or bbmap. You may try different kmer size.         #
#==========================================================#

##Jellyfish 
# jellyfish count -t 8 -C -m 19 -s 5G -o 19mer_out --min-qual-char=? test_R1.fq test_R2.fq
# jellyfish histo -o 19mer_out.histo 19mer_out

##bbmap 
#kmercountexact.sh

suppressMessages(library(argparse))

args <- commandArgs(trailingOnly = TRUE)
hh <- paste(unlist(args),collapse=' ')
listoptions <- unlist(strsplit(hh,'--'))[-1]
options.args <- sapply(listoptions,function(x){
    unlist(strsplit(x, ' '))[-1]
  }, simplify=FALSE
)
options.names <- sapply(listoptions,function(x){
    option <- unlist(strsplit(x, ' '))[1]
})

names(options.args) <- unlist(options.names)

parser <- ArgumentParser(description='Description: estimate a genome size based on the kmer histo.')
parser$add_argument('--table', help='Kmer-based histo table')

if (length(options.args$table) <1 ){
  message()
  parser$print_help()
  message()
  stop()
  quit(save = "default", status = 1, runLast = TRUE)
}

calGenomeSize <- function (table) {
  dataFrame <- read.table(table)
  nPoints <- nrow(dataFrame)
  contamination <- dataFrame[1:10,]
  minDensity <- min(dataFrame[1:10, 2])
  startPoint <- contamination$V1[contamination$V2 == minDensity]
  newDataFrame <- dataFrame[startPoint: 200, ]
  maxDensity <- max(dataFrame[startPoint: 200, 2])
  maxPoint <- newDataFrame$V1[newDataFrame$V2 == maxDensity]
  genomesize <- sum(as.numeric(dataFrame[startPoint:nPoints,1]*dataFrame[startPoint:nPoints,2]))/maxPoint
  return (genomesize)
}

calGenomeSize(options.args$table)
