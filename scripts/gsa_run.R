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

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("parallel"))

# duplicated code from prepare_data.R
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

###############################################################
#
# Assumptions: expgenes, sigGenes, expngenes,nSigGenes exist
compute.stats <- function(exp,genes,pval,fdr,verbose=FALSE) {
  #
  pval.chr <- as.character(pval)
  genes2consider <- genes[genes %in% expgenes[[exp]]]
  if ( length(genes2consider) == 0  ) {
    pwarning("No overlap between the provided gene list and ",exp)
    return(NULL)
  } 
  ngenes <- length(genes2consider)
  overlap <- genes2consider[genes2consider %in% sigGenes[[exp]][[pval.chr]]]
  n.overlap <- length(overlap)
  if ( verbose ) {
    pinfo(exp)
    pinfo("|Genes to consider|:",ngenes)
    pinfo("|Genes DE|:",nSigGenes[[exp]][[pval.chr]])
                                        #pinfo("Genes DE:",sigGenes[[exp]][[as.character(opt$pvalue)]])
                                        #pinfo("GS:",genes2consider)
    pinfo("|DE and GS|:",n.overlap)
    pinfo("|genes|:",exp2ngenes[[exp]])
  }
                                        # DE&GS; DE&not GS; not DE & GS ; not DE & not GS
  go <- matrix(c(n.overlap,nSigGenes[[exp]][[pval.chr]]-n.overlap,
                 ngenes-n.overlap,exp2ngenes[[exp]]-nSigGenes[[exp]][[pval.chr]]-ngenes+n.overlap),
               byrow=T,
               nrow=2)
  #print(go)
  
  ft <- fisher.test(go,alternative="greater")
  gidx <- ""
  # be lazy...only keep the list of genes for the cases that
  # we really may need them
  if ( ft$p.value <= fdr ) {    
    gidx <- paste(which(genes %in% overlap,arr.ind=TRUE),collapse=",")
  }
  res <- data.frame(exp,ft$p.value,n.overlap,ngenes*nSigGenes[[exp]][[pval.chr]]/exp2ngenes[[exp]],gidx)  
  colnames(res) <- c("exp","pval","observed","expected","genes")
  return(res)
}


server.mode <- function(){

  while(TRUE){
    writeLines("Listening...")
    con <- socketConnection(host="localhost", port = 9999, blocking=TRUE,
                            server=TRUE, open="r+")
    on.exit(close.socket(con))
    data <- readLines(con, 1)
    # first word - species
    # second word fdr
    # 3rd word pvalue
    # remaining words - genes
    print(data)
    l <- strsplit(data,split="[ \t,;]+")[[1]]
    print(l)
    species <- l[1]
    opt$fdr <- as.numeric(l[2])
    opt$pvalue <- as.numeric(l[3])
    genes <- l[-c(1,2,3)]

    # uggly code ... duplicated!
    expRes <- list()
    pinfo("Performing tests...")
    expRes <- mclapply(names(expgenes),compute.stats,genes=genes,pval=opt$pvalue,fdr=opt$fdr,verbose=opt$verbose,mc.allow.recursive = FALSE)
#expRes <- lapply(names(expgenes),compute.stats)
    pinfo("Performing tests...done.")
#remove NULLs
#expRes <- expRes[!sapply(expRes,is.null)]
    pinfo("Preparing summary table...")
    final.table <- do.call(rbind.data.frame, expRes)
    final.table$expected <- round(final.table$expected,2)
    final.table$adj.pvalue  <- p.adjust(final.table$pval,method="BH")
    final.table$effect.size <- round(final.table$observed/final.table$expected,2)

# Table complete... 
    pinfo("Table with ",nrow(final.table)," contrasts")
    final.table <- final.table[final.table$adj.pvalue<=opt$fdr,]
    pinfo("Final table with ",nrow(final.table)," significant contrasts where fdr=",opt$fdr)

# return with an exit status of 3 if no overlap is found
    if ( nrow(final.table) ==0 ) {
      writeLines("0", con) 
    } else {

# order by effect size and adj.pvalue
      final.table <- final.table[order(-1*final.table$effect.size,final.table$adj.pvalue,decreasing=F),]

# write to a file
      write.table(final.table,file=con,sep="\t",quote=F,row.names=F)
    }
    close(con)
  }
}

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

pinfo("Loading ",opt$db,"...")
load(opt$db,verbose=T)
pinfo("Loading ",opt$db,"...done.")

pinfo("Number of contrasts: ",length(exp2ngenes))
# TODO: check if all variables are present?

# Number of possible pvals defined when preprocessing the data
#pvals <- c(0.001,0.01,0.05,0.1)
if ( sum(opt$pvalue %in% pvals)!=1 ) {
  perror("Invalid pvalue ",opt$pvalue)
  q(status=1)
}
pinfo("ready to go!")

if ( opt$server ) {
  server.mode()
}
####################################################################
expRes <- list()
pinfo("Performing tests...")
expRes <- mclapply(names(expgenes),compute.stats,genes=genes,pval=opt$pvalue,fdr=opt$fdr,verbose=opt$verbose,mc.allow.recursive = FALSE)
#expRes <- lapply(names(expgenes),compute.stats)
pinfo("Performing tests...done.")
#remove NULLs
#expRes <- expRes[!sapply(expRes,is.null)]
pinfo("Preparing summary table...")
final.table <- do.call(rbind.data.frame, expRes)
final.table$expected <- round(final.table$expected,2)
final.table$adj.pvalue  <- p.adjust(final.table$pval,method="BH")
final.table$effect.size <- round(final.table$observed/final.table$expected,2)

# Table complete... 
pinfo("Table with ",nrow(final.table)," contrasts")
final.table <- final.table[final.table$adj.pvalue<=opt$fdr,]
pinfo("Final table with ",nrow(final.table)," significant contrasts where fdr=",opt$fdr)

# return with an exit status of 3 if no overlap is found
if ( nrow(final.table) ==0 ) {
  cat("No experiments found.\n")
  q(status=3)
}

# order by effect size and adj.pvalue
final.table <- final.table[order(-1*final.table$effect.size,final.table$adj.pvalue,decreasing=F),]

# write to a file
write.table(final.table,file=paste(opt$out.file,".tsv",sep=""),sep="\t",quote=F,row.names=F)

q(status=0)
# not implemented
# make a plot
if ( nrow(final.table) > 1 ) {
  png(file=paste(opt$out.file,".png",sep="")width=550,height=550,res=150)
  dev.off()      
} else {
  pwarning("Insufficient experiments to make a plot")
}
q(status=0)








