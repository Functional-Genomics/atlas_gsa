#!/usr/bin/env Rscript
#; -*- mode: R;-*-
# =========================================================
#
# This is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with iRAP.  If not, see <http://www.gnu.org/licenses/>.
#
#
# =========================================================
# use the bit package
#install.packages("bit")
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("parallel"))
library(bit)

myParseArgs <- function(usage,option_list,filenames.exist=NULL,multiple.options=NULL,mandatory=NULL,...) {

  # get command line options, if help option encountered print help and exit,
  # otherwise if options not found on command line then set defaults,
  parser <- OptionParser(usage = usage, option_list=option_list)
  opt <- parse_args(parser,...)
  
  for ( m in mandatory ) {
    if ( is.null(opt[[m]]) ) {
        perror("Parameter ",m," needs to be defined")
        quit(status=1)
    }
  }
  
  for ( p in filenames.exist ) {
    if (! is.null(opt[[p]]) ) {
      if (! file.exists(opt[[p]]) ) {
        perror("File ",opt[[p]]," not found")
        quit(status=1)
      }
    }
  }

  for ( op in names(multiple.options) ) {
    if ( ! opt[[op]] %in% multiple.options[[op]] ) {
      perror("Invalid value ",opt[[op]]," for option ",op)
      quit(status=1)
    }
  }
  return(opt)
}


readList <- function(filename) {
  df <- NULL
  df <- tryCatch(read.table(filename,header=F,colClasses=c("character")),error=function(x) NULL)
  if ( is.null(df) ) {
    stop("Error loading ",filename)
  }
  return(df$V1)
}

load.tsv <- function(filename) {

  df <- NULL
  df <- tryCatch(read.table(filename,header=T,sep="\t",check.names=F,comment.char="",quote=""))
  if ( is.null(df) ) {
    perror("Error loading ",filename)
    quit(status=1)
  }
  return(df)
}


############################################################
# Useful functions
pinfo <- function(...) {
  cat(paste("[INFO] ",...,"\n",sep=""))
}

pwarning <- function(...) {
  cat("[WARNING] ",...,"\n",file=stderr())
}

perror <- function(...) {
  cat("[ERROR] ",...,"\n",file=stderr())
}

debug <- FALSE
pdebug <- function(...) {
  if (debug) 
    cat("[ERROR] ",...,"\n",file=stderr())
}
###############################################################################
usage <- "prepare_data.R --tsv-list file --out prefix"
option_list <- list(
  make_option(c("-c", "--cores"), type="character",default="2",dest="num_cores",help="Number of cores to use ([default %default])"),
  make_option(c("-i", "--tsv-list"), type="character", dest="tsv.list.file", default=NULL,help="File with the list of TSV file names  (gene ids should appear in the first column)"),
  make_option(c("-o", "--out"), type="character",default=NULL,dest="out.file",help="Output file name prefix."),
  make_option(c("-d", "--debug"), ,action="store_true",default=FALSE,dest="debug",help="Enable debug mode (much more verbose).")
)
multiple.options <- NULL
filenames <- c("tsv.list.file") ;#filenames that must exist (if defined)

# check multiple options values
mandatory <- c("tsv.list.file","out.file")
opt <- myParseArgs(usage = usage, option_list=option_list,filenames.exist=filenames,multiple.options=multiple.options,mandatory=mandatory)

debug <- opt$debug

pinfo("Parameters parsed.")
for  ( n in names(opt)) {
      cat(paste(n,"=",opt[n]),"\n",sep="")
}

# cores to use
tryCatch(num.cores <- as.integer(as.numeric(opt$num_cores)),warning=
         function(w) {
           stop("Invalid number of cores ",opt$num_cores)
       }
)
if (num.cores<1) {
  stop("Invalid number of cores ",opt$num_cores)
}

if ( num.cores>detectCores()) {
  num.cores <- detectCores()
  pwarning("The number of cores to use exceeds the cores available. Reducing the limit to ",detectCores())
}

options(mc.cores=num.cores)

##############################
pinfo("Reading the list of tsv file names",opt$tsv.list.file)
tsv.files <- readList(opt$tsv.list.file)
pinfo("Number of TSV files to process: ",length(tsv.files))

# TODO: check for duplicated files
################################
# data
pvals <- c(0.01,0.05,0.1)
gene.name.col <- 1
# number of significant genes and non-significant genes
nSigGenes <- list()
sigGenes <- list()
# Total number of genes considered in an experiment:::contrast
exp2ngenes <- list()
# Set of genes considered on each experiment
expgenes <- list()

nfiles <- length(tsv.files)
sort.pvals <- sort(pvals)
for ( i in seq(1,nfiles)) {
  #cat(i,"\n")
  pinfo("File:",tsv.files[i])
  tsv.df <- load.tsv(tsv.files[i])
  pinfo("File ",tsv.files[i]," loaded.")
  # Each file may contain different comparisons
  colsOfInterest <- grep("p-value",colnames(tsv.df),value=T)
  pinfo("Found ",length(colsOfInterest)," comparisons by looking at p-value columns")
  exp <- gsub("-analytics.*","",basename(tsv.files[[i]]))
  for ( df.col in colsOfInterest) {
    cont <- gsub(".p-value","",df.col)
    ## TODO: replace tsv.file by the name/accession of the experiment
    ## TODO: how to deal with the NAs?? -> pvalue=1
    
    key=paste(exp,":::",cont,sep="")
    #pinfo("Contrast: ", cont)
    pinfo("key: ", key)
    # create a matrix for each contrast
    sigGenes[[key]] <- list()
    nSigGenes[[key]] <- rep(NA,length(pvals))
    names(nSigGenes[[key]]) <- as.character(pvals)

    exp2ngenes[[key]] <- nrow(tsv.df)
    expgenes[[key]] <- as.character(tsv.df[,gene.name.col])
    pinfo("Number of genes: ",exp2ngenes[[key]])
    #cat("PVALS:",as.character(pvals))
    #print(sigGenes[[key]])
    for (pval in sort.pvals) {
      # TODO: optimize memory usage
      pdebug("pval:",pval)
      #pinfo(df.col)
      #sg <- sum(tsv.df[,df.col]<=pval,na.rm=TRUE)
      sg <- tsv.df[tsv.df[,df.col]<=pval,gene.name.col]
      # exclude NA
      sg <- as.character(sg[!is.na(sg)])
      nsigGenes <- length(sg)
      pdebug("sigGenes: ",sg)
      #save.image("test.Rdata")
      #load("test.Rdata")
      pdebug("nsigGenes: ",nsigGenes)
      #pinfo(names(sigGenes[[key]]))
      sigGenes[[key]][[as.character(pval)]] <- sg
      nSigGenes[[key]][as.character(pval)] <- nsigGenes
    }
    names(sigGenes[[key]]) <- as.character(pvals)    
  }
  pinfo(round(i/nfiles*100,2),"% done")
}

#print(sigGenes)
#print(exp2ngenes)
pinfo("Saving to file ",opt$out.file)
save(pvals,nSigGenes,sigGenes,exp2ngenes,expgenes,file=opt$out.file)
# Make a few plots
# Distribution of the number of genes per experiment
#
pinfo("That's all folks!")
q(status=0)

