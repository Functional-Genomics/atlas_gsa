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
suppressPackageStartupMessages(library("bit"))

###############################################################
# 
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
# Useful IO functions
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
  #if ( ft$p.value <= fdr ) {    
  #  gidx <- paste(which(genes %in% overlap,arr.ind=TRUE),collapse=",")
  #}
  res <- data.frame(exp,ft$p.value,n.overlap,ngenes*nSigGenes[[exp]][[pval.chr]]/exp2ngenes[[exp]])  
  colnames(res) <- c("exp","pval","observed","expected")
  return(res)
}

###########################################
# Get a bit vector from a set of gene names
gene.list2bit <- function(genes,ref.genes) {  
  tf <- ref.genes %in% genes
  return(as.bit(tf))
}

# global var:  exp.index 
# user's "genes" is already a bit vector
compute.stats.bit <- function(exp,genes.b,pval,fdr,verbose=FALSE) {
                                        #
  stopifnot(is.bit(genes.b))
  pval.chr <- as.character(pval)
  # overlap between the user's set of genes and the background set of genes
  genes2consider <- genes.b &  exp.index$expgenes[[exp]]
  ngenes <- sum(genes2consider)
  if ( ngenes == 0  ) {
    pwarning("No overlap between the provided gene list and ",exp)
    return(NULL)
  }
  # 
  overlap <- genes2consider & exp.index$sigGenes[[exp]][[pval.chr]]
  n.overlap <- sum(overlap)
  if ( verbose ) {
    pinfo(exp)
    pinfo("|Genes to consider|:",ngenes)
    pinfo("|Genes DE|:",exp.index$nSigGenes[[exp]][[pval.chr]])
                                        #pinfo("Genes DE:",sigGenes[[exp]][[as.character(opt$pvalue)]])
                                        #pinfo("GS:",genes2consider)
    pinfo("|DE and GS|:",n.overlap)
    pinfo("|genes|:",exp.index.$exp2ngenes[[exp]])
  }
                                        # DE&GS; DE&not GS; not DE & GS ; not DE & not GS
  go <- matrix(c(n.overlap,exp.index$nSigGenes[[exp]][[pval.chr]]-n.overlap,
                 ngenes-n.overlap,exp.index$exp2ngenes[[exp]]-exp.index$nSigGenes[[exp]][[pval.chr]]-ngenes+n.overlap),
               byrow=T,
               nrow=2)
  #print(go)  
  ft <- fisher.test(go,alternative="greater")
  gidx <- ""
  # be lazy...only keep the list of genes for the cases that
  # we may need them
  #if ( ft$p.value <= fdr ) {    
  #  gidx <- paste(which(genes %in% overlap,arr.ind=TRUE),collapse=",")
  #}
  res <- data.frame(exp,ft$p.value,n.overlap,ngenes*exp.index$nSigGenes[[exp]][[pval.chr]]/exp.index$exp2ngenes[[exp]])  
  colnames(res) <- c("exp","pval","observed","expected")
  return(res)
}

###################################################
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

###########################################################


run.GSA <- function(exps,FUN,pval,genes,fdr,...) {

  expRes <- list()
  pinfo("Performing tests...")
  expRes <- mclapply(exps,FUN,genes=genes,pval=pval,fdr=fdr,...,mc.allow.recursive = FALSE)

  pinfo("Preparing summary table...")
  final.table <- do.call(rbind.data.frame, expRes)
  #print(final.table)
  final.table$expected <- round(final.table$expected,2)
  final.table$adj.pvalue  <- p.adjust(final.table$pval,method="BH")
  final.table$effect.size <- round(final.table$observed/final.table$expected,2)
  
  # Table complete... 
  pinfo("Table with ",nrow(final.table)," contrasts")
  final.table <- final.table[final.table$adj.pvalue<=opt$fdr,]
  pinfo("Final table with ",nrow(final.table)," significant contrasts for fdr=",opt$fdr)

  # return NULL
  if ( nrow(final.table) ==0 ) {
    cat("No experiments found.\n")
    return(NULL)
  }

# order by effect size and adj.pvalue
  final.table <- final.table[order(-1*final.table$effect.size,final.table$adj.pvalue,decreasing=F),]
  return(final.table)
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
    pinfo("Performing tests...done.")
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
