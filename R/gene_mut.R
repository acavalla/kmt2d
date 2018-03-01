### freads mafs and greps lines containing user-defined transcript. In MAF files, there is a main gene affected, and then other possible
## effects in a single field all bunched together. If the transcr of interest is not the main gene affected, function searches within
## the single field and retrieves the correct transcript and effect. Also fixes the issue of thymine (T) being interpreted as "TRUE" and 
## removes "chr" from before chromosomes, as found in GRCh38-annotated files.

library(data.table)
library(plyr)
process.maf <- function(x){
  a <- file.names[i]
  a <- paste0("/projects/rmorin/projects/nhl_meta_analysis/results/strelka_pipeline-1.0.0/", a)
  x <- tryCatch(fread(paste("grep", transcr[j], file.names[i]), header=FALSE, sep="\t", stringsAsFactors=FALSE), error=function(e) NULL)
  if(!is.null(x)){
    x <- as.data.frame(as.matrix(x))
    x <- x[,which(unlist(lapply(x, function(x)!all(is.na(x))))),]
    x$file <- file.names[i]
    x$i <- i
    y <- x[!grepl(transcr[j], x$V38),]
    x <- x[grepl(transcr[j], x$V38),]
    if(nrow(y) != 0){
      y$V46 <- as.character(y$V46)
      aa <- strsplit(y$V46, ";", fixed = FALSE, perl = FALSE, useBytes = FALSE)
      aa <- lapply(aa, function(ch) grep(transcr[j], ch, value = TRUE))
      aa <- strsplit(as.character(aa), ",", fixed = FALSE, perl = FALSE, useBytes = FALSE)
      y$V1 <- as.character(y$V1)
      y$V9 <- as.character(y$V9)
      if(!is.null(y$V37)){
        y$V37 <- as.character(y$V37)
      }
      y$V38 <- as.character(y$V38)
      for(k in 1:length(aa)){
        y$V1[k] <- aa[[k]][1]
        y$V9[k] <- aa[[k]][2]
        y$V37[k] <- aa[[k]][3]
        y$V38[k] <- aa[[k]][4]
      }
      if(length(x$V37)!=length(y$V37)){
        x$V37 <- rep("NA", nrow(x))
      }
      x <- rbind(y, x)
      x <- x[order(x$file, x$V6),] 
      x$V9 <- gsub("_variant$", "", x$V9)
    }
    for(k in 1:nrow(x)){
      if(!is.null(x$V37[k])){
        x$V37 <- as.character(x$V37)
        x$V37[k] <- ifelse(x$V37[k]=="", "NA", x$V37[k])
      } else {
        x$V37[k] <- "NA"
        x$V37 <- as.character(x$V37)
      }
    }
    
    x <- cbind(x$file, x[1], x[4:7], x[11:13], x$V37, x$V9, x$i)
    x$V11 <- as.factor(revalue(warn_missing=F, x$V11, c("TRUE" = "T")))
    x$V12 <- as.factor(revalue(warn_missing=F, x$V12, c("TRUE" = "T")))
    x$V13 <- as.factor(revalue(warn_missing=F, x$V13, c("TRUE" = "T")))
    x$V5 <- gsub("^chr", "", x$V5)
  }
  return(x)
}

##transcript of your choice 
transcr <- "ENST00000301067" #KMT2D

##file names within /projects/rmorin/projects/nhl_meta_analysis/results/strelka_pipeline-1.0.0/ listed in a file created by /projects/acavalla_prj/kmt2d/301067/patient_data/files/file.extract.R
file.names <- read.table("/projects/acavalla_prj/kmt2d/301067/patient_data/files/filenames.maf.txt")
setwd("/projects/rmorin/projects/nhl_meta_analysis/results/strelka_pipeline-1.0.0/")

datalist <- list()
for(j in 1:length(transcr)){
  for(i in 1:length(file.names)){
    x <- process.maf(file.names[i])
    datalist[[i]] <- x ##consecutively add it to list 
    if(i%%10==0){
      print(i)
    }
  }
}

datalist <- datalist[-which(sapply(datalist, is.null))]
file <- do.call(rbind, datalist)
colnames(file) <- c("file", "gene", "assembly", "chr", "startpos", "endpos", "ref", "alt", "alt2", "aachange", "mut.type", "i")

write.table(file, file=paste("/projects/acavalla_prj/kmt2d/301067/patient_data/gene_muts/", transcr[j], "_mutations.txt", sep=""), quote=FALSE, sep="\t")
