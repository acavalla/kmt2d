### Deconstruct Sigs to identify mutational processes
# Input required:
# sample.id Column name in the mutation file corresponding to the Sample ID
# chr Column name in the mutation file corresponding to the chromosome
# pos Column name in the mutation file corresponding to the mutation position
# ref Column name in the mutation file corresponding to the reference base
# alt Column name in the mutation file corresponding to the alternate base

install.packages("deconstructSigs")
source("https://bioconductor.org/biocLite.R")
biocLite("BSgenome.Hsapiens.UCSC.hg18")
biocLite("BSgenome.Hsapiens.NCBI.GRCh38")
library(deconstructSigs)
library('VariantAnnotation')
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Hsapiens.UCSC.hg18)
library(dplyr)
install.packages("dendsort")

setwd("/projects/rmorin/projects/nhl_meta_analysis/results/strelka_pipeline-1.0.0/")
post.process <- function(x) {
  x$pos <- as.integer(x$pos)
  x <- x[,1:5]
  # remove factors
  x$sample.id <- factor(x$sample.id)
  x$samples <- as.numeric(x$sample.id)
  return(x)
}
v2df = function(x) {
  df = as.data.frame(rowRanges(x))
  if (nrow(df) == 0) {
    print(snv.vcf[i])
  } else {
    df = cbind(df, as.data.frame(info(x)))
    if ( any(c('ANN', 'EFF') %in% names(info(x))) ) {
      ann = c('ANN', 'EFF')[ c('ANN', 'EFF') %in% names(info(x)) ][1]
      dfann = .anncols(df$ANN, info(header(x))[ann, ]$Description)
      df = df[, colnames(df) != ann]
      df = cbind(df, dfann)
    }
    n  = names(geno(x))
    tmp = lapply(n, function(col) {
      return(as.data.frame(geno(x)[[col]]))
    })
    ncols = sapply(tmp, ncol)
    tmp = do.call(cbind, tmp)
    colnames(tmp) = paste(rep(n, times = ncols), colnames(tmp), sep =
                            "_")
    df = cbind(df, tmp)
    df[sapply(df, is.list)] <- apply(df[sapply(df, is.list)], 2, 
                                     function(x) { 
                                       unlist(lapply(x, paste,  
                                                     sep=",", collapse=";"))
                                     } )
    print(dim(df))
  }
  x <- df
  return(x)
}

sig <- list()
for(i in 1:length(snv.vcf)){
  x <- readVcf(snv.vcf[i])
  x <- v2df(x)
  if (nrow(x) == 0) {
    print(snv.vcf[i])
  } else {
    x <- as.data.frame(cbind("sample.id" = paste0(snv.vcf[i]), "chr" = x$seqnames, "pos" = x$start, "ref" = x$REF, "ALT" = x$ALT))
    x$chr <- as.numeric(as.character(x$chr))
    x<-x[(x$chr<23),]
    x$chr <- sub("^", "chr", x$chr)
    x$ref <- as.character(x$ref)
    x$pos <- as.numeric(as.character(x$pos))
    x$alt <- as.character(x$ALT)
    x <- x[complete.cases(x), ]
    x$i <- i # maybe you want to keep track of which iteration produced it?
    sig[[i]] <- x # add it to your list
    print(i)
  }
}

sig <- do.call(rbind, sig)
sig.done <- post.process(sig)

write.table(sig.done, "/projects/acavalla_prj/kmt2d/301067/patient_data/mutsigs/v2dfrunvcfsnvs.txt")

sigall <- sig.done
sigall$ave <- as.numeric(as.factor(as.character(sigall$sample.id)))
sigall$ALT <- as.character(sigall$ALT)
s <- strsplit(sigall$ALT, split = ";")
sigall <- data.frame(sample.id = rep(sigall$sample.id, sapply(s, length)), chr = rep(sigall$chr, sapply(s, length)), pos = rep(sigall$pos, sapply(s, length)), ref = rep(sigall$ref, sapply(s, length)), ALT = unlist(s), samples = rep(sigall$samples, sapply(s, length)), ave = rep(sigall$ave, sapply(s, length)))
sigall$chr <- as.character(sigall$chr)
sigall$ref <- as.character(sigall$ref)
sigall$ALT <- as.character(sigall$ALT)
names(sigall)[5] <- "alt"

trimers <- list()
for(i in 1:nrow(sigall)){
  sigalli <- sigall[sigall$samples==i,]
  sigalli$ave <- 1
  a <- if(grepl("BL_Love|MCL", sigalli[1,1])){
    BSgenome.Hsapiens.UCSC.hg38
  } else {
    BSgenome.Hsapiens.UCSC.hg19
  }
  
  trimers[[i]] <-  mut.to.sigs.input(mut.ref = sigalli,
                                     sample.id = "ave",
                                     chr = "chr",
                                     pos = "pos",
                                     ref = "ref", 
                                     alt = "alt",
                                     bsg = a)
  print(paste(i, "/", max(sigall$samples)))
}

trimers <- do.call(rbind, trimers)
rownames(trimers) <- 1:nrow(trimers)
write.table(trimers, "/projects/acavalla_prj/kmt2d/301067/patient_data/mutsigs/trimers.tsv")

percents <- list()
for (i in 1:nrow(trimers)) {
  sample <- whichSignatures(tumor.ref = trimers, 
                            signatures.ref = signatures.cosmic, 
                            sample.id = i, 
                            contexts.needed = TRUE,
                            tri.counts.method = 'genome')
  sample$i <- i  # maybe you want to keep track of which iteration produced it?
  percents[[i]] <- sample 
  print(i)
}

barchart <- percents[[1]]$weights
for(i in 2:length(percents)) {
  barchart <- rbind(barchart, percents[[i]]$weights)
}

write.table(barchart, "/projects/acavalla_prj/kmt2d/301067/patient_data/mutsigs/sigproportions")


#smry.mutsigs order:
# 1:6 = BLIC (6)
# 7:19 = BLLD (13) (PT1092 removed)
# 20:163 = CLL (144)
# 164:275 = DLBCL (112)
# 276:397 = FL (not inc Fitzgibbon) (122)
# 398:597 = MCL (two Zhangs and M006 removed) (200) -- arranged alphabetically

smry.mutsigs <- read.table("/projects/acavalla_prj/kmt2d/301067/patient_data/summarybysample.txt")
smry.mutsigs <- smry.mutsigs[smry.mutsigs$snvs!=0,]
rownames(smry.mutsigs) <- 1:nrow(smry.mutsigs)

trimers <- read.table("/projects/acavalla_prj/kmt2d/301067/patient_data/mutsigs/trimers.tsv")
sigs <- read.table("/projects/acavalla_prj/kmt2d/301067/patient_data/mutsigs/sigproportions")
types <- list(trimers,sigs)


library(pheatmap)
library(viridis)
library(dendsort)
library(RColorBrewer)
sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

for(i in 1:length(types)){
  pheat <- types[i]
  pheat <- as.matrix(types[[i]])
  #normalise pheat to 1 
  pheatnorm <- prop.table(pheat,1) 
  row.names(pheatnorm) <- 1:nrow(pheatnorm)
  
  #split into types
  a <- rownames(smry.mutsigs[grepl("BL", smry.mutsigs$mafindel),])
  b <- rownames(smry.mutsigs[grepl("CLL", smry.mutsigs$mafindel),])
  c <- rownames(smry.mutsigs[grepl("DLBCL", smry.mutsigs$mafindel),])
  d <- rownames(smry.mutsigs[grepl("FL", smry.mutsigs$mafindel),])
  e <- rownames(smry.mutsigs[grepl("MCL", smry.mutsigs$mafindel),])
  
  BLpheat <- pheatnorm[a,]
  CLLpheat <- pheatnorm[b,]
  DLBCLpheat <- pheatnorm[c,]
  FLpheat <- pheatnorm[d,]
  MCLpheat <- pheatnorm[e,]
  
  pheats <- list(BLpheat, CLLpheat, DLBCLpheat, FLpheat, MCLpheat)
  ids <- list(a, b, c, d, e)
  for(j in 1:length(pheats)){
    mat_breaks <- quantile_breaks(pheats[[j]], n = 11)
    col_groups <- c(as.character(smry.mutsigs[ids[[j]][1]:ids[[j]][length(ids[[j]])],10]))
    mat_col <- data.frame(group = col_groups)
    rownames(mat_col) <- colnames(t(pheats[[j]]))
    mat_colors <- list(group = brewer.pal(3, "Set3"))
    names(mat_colors$group) <- unique(col_groups)
    mat_cluster_cols <- sort_hclust(hclust(dist(pheats[[j]])))
    mat_cluster_rows <- sort_hclust(hclust(dist(t(pheats[[j]]))))
    
    jpeg(paste0("/projects/acavalla_prj/kmt2d/301067/patient_data/mutsigs/", types[i], "_", pheat[j], ".jpg"), height = 1800, width = 2500)
    pheatmap(
      mat               = t(pheats[[j]]),
      color             = viridis(length(mat_breaks) - 1),
      breaks            = mat_breaks,
      border_color      = NA,
      cluster_cols      = mat_cluster_cols,
      cluster_rows      = mat_cluster_rows,
      show_rownames     = TRUE,
      show_colnames = FALSE,
      fontsize_col = 2,
      annotation_col    = mat_col,
      annotation_colors = mat_colors,
      drop_levels       = TRUE,
      fontsize          = 8,
      main              = x
    )
    dev.off()
    print(paste("j =", j))
  }
  print(paste("i =", i))
}
