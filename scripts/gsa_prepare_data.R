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

initial.options <- commandArgs(trailingOnly = FALSE)
script.dir <- dirname(sub("--file=", "", initial.options[grep("--file", initial.options)]))
source(paste(script.dir,"/../R/atlasGSA.R",sep=""))

#library(atlasGSA)

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("parallel"))

###############################################################################
usage <- "prepare_data.R --tsv-list file  --out prefix"
option_list <- list(
  make_option(c("-c", "--cores"), type="character",default="2",dest="num_cores",help="Number of cores to use ([default %default])"),
  make_option(c("-i", "--tsv-list"), type="character", dest="tsv.list.file", default=NULL,help="File with the list of TSV file names  (gene ids should appear in the first column)"),
  make_option(c("-x", "--xml-list"), type="character", dest="xml.list.file", default=NULL,help="File with the list of XML configuration file names (deprecated)."),
  make_option(c("-o", "--out"), type="character",default=NULL,dest="out.file",help="Output file name prefix."),
  make_option(c("-d", "--debug"), ,action="store_true",default=FALSE,dest="debug",help="Enable debug mode (much more verbose).")
)
multiple.options <- NULL
filenames <- c("tsv.list.file","xml.list.file") ;#filenames that must exist (if defined)

# check multiple options values
mandatory <- c("tsv.list.file","xml.list.file","out.file")
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
if ( ! is.null(opt$xml.list.file) ) {
  pinfo("Reading the configuration file names",opt$xml.list.file)
  xml.files <- readList(opt$xml.list.file)
  pinfo("Number of XML files to process: ",length(xml.files))
}

#if ( ! is.null(opt$xml.list.file) ) {
# Sanity check
#  if ( length(xml.files) != length(tsv.files) ) {
#    perror("Number of xml files (",length(xml.files),") is different from the number of tsv files (",length(tsv.files),")")
#    quit(status=1)
#  }
#}

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
    # some exps (MA) have the same gene multiple times E-GEOD-3076_A-AFFY-27:::g13_g18
    expgenes[[key]] <- unique(as.character(tsv.df[,gene.name.col]))
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
      sg <- unique(as.character(sg[!is.na(sg)]))
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
#pinfo("Saving to file ",opt$out.file)
#save(pvals,nSigGenes,sigGenes,exp2ngenes,expgenes,file=opt$out.file)
# Make a few plots
# Distribution of the number of genes per experiment
#

#############################################
# optimization attempt
# install.packages(bit)
library(bit)

gene.list2bit <- function(genes,ref.genes) {  
  tf <- ref.genes %in% genes
  return(as.bit(tf))
}
# how many and which genes do we have?
genes <- c()
for ( n in names(expgenes) ) {
  if ( !is.null(expgenes[[n]]) ) {
    genes <- unique(append(genes,expgenes[[n]]))
  }
}
#genes <- sort(unique(unlist(expgenes)))
pinfo("#genes:",length(genes))

exp.index <- list()
# all genes for a given species
exp.index$genes <- genes
# number of genes
exp.index$ngenes <- length(genes)
# number of genes per contrast
exp.index$exp2ngenes <- exp2ngenes
# set of genes in each contrast
exp.index$expgenes <- lapply(expgenes,gene.list2bit,ref.genes=exp.index$genes)
# p-values considered
exp.index$pvals <- pvals
# number of significant genes per contrast/pvalue
exp.index$nSigGenes <- nSigGenes
# set of sig. genes per contrast/pvalue
exp.index$sigGenes <- list()
for ( c in names(sigGenes) ) {
  exp.index$sigGenes[[c]] <- list()
  for ( p in names(sigGenes[[c]])) {
    exp.index$sigGenes[[c]][[p]] <- gene.list2bit(sigGenes[[c]][[p]],ref.genes=exp.index$genes)
  }
}
pinfo("Saving index to file ",paste(opt$out.file,sep=""))
save(exp.index,file=paste(opt$out.file,sep=""))
pinfo("That's all folks!")
q(status=0)
