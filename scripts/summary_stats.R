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
#library(atlasGSA)
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("parallel"))

#args <- c("v2_Anas_platyrhynchos.po","v2_Bos_taurus.po")
args <- c("Anas_platyrhynchos.po","Arabidopsis_thaliana.po","Bos_taurus.po","Brassica_oleracea.po","Caenorhabditis_elegans.po","Danio_rerio.po","Drosophila_melanogaster.po","Gallus_gallus.po","Glycine_max.po","Homo_sapiens.po","Hordeum_vulgare.po","Macaca_mulatta.po","Medicago_truncatula.po","Mus_musculus.po","Oryza_sativa.po","Ovis_aries.po","Physcomitrella_patens.po","Populus_trichocarpa.po","Rattus_norvegicus.po","Saccharomyces_cerevisiae.po","Solanum_tuberosum.po","Sorghum_bicolor.po","Sus_scrofa.po","Triticum_aestivum.po","v2_Anas_platyrhynchos.po","v2_Arabidopsis_thaliana.po","v2_Bos_taurus.po","v2_Brassica_oleracea.po","v2_Caenorhabditis_elegans.po","v2_Danio_rerio.po","v2_Drosophila_melanogaster.po","v2_Gallus_gallus.po","v2_Glycine_max.po","v2_Homo_sapiens.po","v2_Hordeum_vulgare.po","v2_Macaca_mulatta.po","v2_Medicago_truncatula.po","v2_Mus_musculus.po","v2_Oryza_sativa.po","v2_Ovis_aries.po","v2_Physcomitrella_patens.po","v2_Populus_trichocarpa.po","v2_Rattus_norvegicus.po","v2_Saccharomyces_cerevisiae.po","v2_Solanum_tuberosum.po","v2_Sorghum_bicolor.po","v2_Sus_scrofa.po","v2_Triticum_aestivum.po","v2_Vitis_vinifera.po","v2_Zea_mays.po","Vitis_vinifera.po","Zea_mays.po")
#args <- commandArgs(trailingOnly = TRUE)

ngenesbg <- list()
nde <- list()
contrasts <- list()
de.pval <- "0.05"

for ( db.file in args ) {                                     #
  #db.file <- "v2_Anas_platyrhynchos.po"
  cat("Loading ",db.file,"...")
  load(db.file)
  ppfilename <- gsub("_"," ",gsub("v2_","",gsub(".po","",db.file)))
  ngenesbg[[ppfilename]] <- unlist(exp.index[["exp2ngenes"]])
  contrasts[[ppfilename]] <- names(exp.index[["exp2ngenes"]])
  nde[[ppfilename]] <- as.numeric(data.frame(exp.index$nSigGenes)[de.pval,])
  cat("done.\n")
}

#
save.image("sum_data.RData")
#setwd("/home/nf/Research/Projects/WIP/Atlas/BBR_Software/atlas_gsa/tests/stats")
#load(".RData")


####################
# Create the plots

get.norm.nde <- function(nde,ngenesbg) {
  do.comp <- function(l,v1,v2) {
    return(round(v1[[l]]/v2[[l]]*100,2))
  }
  x <- lapply(names(ngenesbg),do.comp,nde,ngenesbg)
  names(x) <- names(ngenesbg)
  return(x)
}
norm.nde <- get.norm.nde(nde,ngenesbg)


png("summary_de_genes_per_specie.png",width=1000,height=1000,res=150)
s.order <- sort(names(nde),decreasing=T)
par(bty="l",mar=c(5,13,5,3))
boxplot(nde[s.order],horizontal=T,las=2,main="Differentially Expressed Genes\n(a=0.05)")
dev.off()


png("summary_norm_de_genes_per_specie.png",width=1000,height=1000,res=150)
s.order <- sort(names(nde),decreasing=T)
par(bty="l",mar=c(5,13,5,3))
boxplot(norm.nde[s.order],horizontal=T,las=2,main="Percentage of Genes\n Differentially Expressed Genes\n(a=0.05)",xlab="%")
dev.off()


png("summary_genes_per_specie.png",width=1000,height=1000,res=150)
s.order <- sort(names(nde),decreasing=T)
par(bty="l",mar=c(5,13,5,3))
boxplot(ngenesbg[s.order],horizontal=T,las=2,main="Number of Genes Considered",log="x")
dev.off()


png("summary_contrasts_per_specie.png",width=1000,height=1000,res=150)
par(bty="l",mar=c(5,13,5,3))
vals <- as.numeric(lapply(contrasts,length))
names(vals) <- names(contrasts)
vals <- sort(vals,decreasing=TRUE)
barplot(vals,las=2,horiz=T,main="Contrasts per Specie")
dev.off()

q(status=0)








