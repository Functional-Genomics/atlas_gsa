#!/usr/bin/env Rscript
#; -*- mode: R;-*-
# =========================================================
#
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

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library("bit"))


###############################################################################
usage <- "gsa_run.R --db db.po --gs 'gene1 gene2 ...'  [--pvalue 0.01] --out prefix [-c cores]"
# note: the number of genes should be smaller than ~800 (shell constraints)
option_list <- list(
  make_option(c("-c", "--cores"), type="character",default="2",dest="num_cores",help="Number of cores to use ([default %default])"),
  make_option(c("-p", "--pvalue"), type="numeric",default=0.01,dest="pvalue",help="Pvalue threshold ([default %default])."),
  make_option(c("-f", "--fdr"), type="numeric",default=0.01,dest="fdr",help="FDR ([default %default])."),
  make_option(c("-d", "--db"), type="character", dest="db.file", default=NULL,help="Path to a file created with prepare_data.R."),
  make_option(c("-g", "--gs"), type="character",default=NULL,dest="gs",help="Set of genes."),
  make_option(c("-v", "--verbose"), ,action="store_true",default=FALSE,dest="verbose",help="Increase verbosity."),
  make_option(c("-s", "--server"), ,action="store_true",default=FALSE,dest="server_mode",help="Server mode (listen in port 9999)."),
  make_option(c("-o", "--out"), type="character",default=NULL,dest="out.file",help="Output file name prefix.")
)
multiple.options <- NULL
filenames <- c("db.file") ;#filenames that must exist (if defined)

# check multiple options values
mandatory <- c("db.file","out.file","gs")
opt <- myParseArgs(usage = usage, option_list=option_list,filenames.exist=filenames,multiple.options=multiple.options,mandatory=mandatory)

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

###############################################################

pinfo("p-value=",opt$pvalue)
genes <- strsplit(opt$gs,split="[ \t;,]+")[[1]]
pinfo("|Set of genes provided|=",length(genes))
pinfo("|Set of unique genes provided|=",length(unique(genes)))
genes <- unique(genes)

pinfo("Loading ",opt$db,"...")
system.time(load(opt$db,verbose=T))
pinfo("Loading ",opt$db,"...done.")

pinfo("Number of contrasts: ",length(exp.index$exp2ngenes))
# TODO: check if all variables are present?

# Number of possible pvals defined when preprocessing the data
#pvals <- c(0.001,0.01,0.05,0.1)
if ( sum(opt$pvalue %in% exp.index$pvals)!=1 ) {
  perror("Invalid pvalue ",opt$pvalue)
  q(status=1)
}

# convert genes to bitmap
genes.b <- gene.list2bit(genes,exp.index$genes)
pinfo("genes.b:",sum(genes.b))
pinfo("ready to go!")

if ( opt$server ) {
  server.mode()
}
####################################################################
final.table <- run.GSA(names(exp.index$expgenes),FUN=compute.stats.bit,genes=genes.b,pval=opt$pvalue,fdr=opt$fdr,verbose=opt$verbose)

# nothing found
if (is.null(final.table) ) {
  q(status=3)
}
# write to a file
write.table(final.table,file=paste(opt$out.file,".tsv",sep=""),sep="\t",quote=F,row.names=F)

q(status=0)

