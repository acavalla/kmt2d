library(data.table)

##load files of KMT2D mutation status
kmt2d.all <- read.table("/projects/acavalla_prj/kmt2d/301067/patient_data/kmt2d_muts.tsv")
kmt2d.all <- kmt2d.all[!grepl("MCL",kmt2d.all$folder),]

kmt2c.all <- read.table("/projects/acavalla_prj/kmt2d/301067/patient_data/kmt2c_muts.tsv")
kmt2c.all <- kmt2c.all[!grepl("MCL",kmt2c.all$folder),]

parp1.all <- read.table("/projects/acavalla_prj/kmt2d/301067/patient_data/parp1_muts.tsv")
parp1.all <- parp1.all[!grepl("MCL",parp1.all$folder),]

effect.files <- c("/projects/acavalla_prj/kmt2d/301067/patient_data/provean/kmt2d.result.tsv", "/projects/acavalla_prj/kmt2d/301067/patient_data/provean/kmt2c.result.tsv",
                  "/projects/acavalla_prj/kmt2d/301067/patient_data/provean/parp1.result.tsv")
effect.list <- list()
clean.effects <- function(x) {
  x <- read.csv(effect.files[i], sep = "\t")
  x <- x[!duplicated(x[1]),]
  x <- x[,c(2,3,8:12)]
  colnames(x)[7] <- "PREDICTION"
  x$PREDICTION <- as.character(x$PREDICTION)
  return(x)
}
for (i in 1:length(effect.files)){
  x <- clean.effects(x)
  effect.list[[i]] <- x
}

kmt2d.effects <- cbind(kmt2d.all, effect.list[1])
kmt2c.effects <- cbind(kmt2c.all, effect.list[2])
parp1.effects <- cbind(parp1.all, effect.list[3])

x <- kmt2d.effects
clean <- function(x){
  x <- x[x$PREDICTION != "Neutral",]
  x <- x[x$mut.type != "utr_3",]
  x <- x[x$mut.type != "utr_5",]
  x <- x[x$mut.type != "intron",]
  x <- x[x$mut.type != "downstream_variant",]
  x <- x[x$mut.type != "upstream_variant",]
  x <- x[x$gene != NA,]
}

x <- clean(x)
x <- x[complete.cases(x),]
kmt2d.effects <- clean(kmt2d.effects)
kmt2c.effects <- clean(kmt2c.effects)
parp1.effects <- clean(parp1.effects)


snv <- read.table("/projects/acavalla_prj/kmt2d/301067/patient_data/all_snvs/all_snvs.tsv", sep="\t")

samples <- data.frame("iteration" = 1:16, "samples" = c("BL_ICGC.vcf", "BL_Love_Dave_2012.vcf", "CLL_Connors.vcf", "CLL_ICGC.vcf", "DLBCL_Connors.vcf", "DLBCL_ICGC.vcf", "DLBCL_refractory.vcf", "DLBCL_TCGA.vcf", "FL_Connors.vcf",
                                                        "FL_Shah.vcf", "MCL_Bea_Campo.vcf", "MCL_Morin_2015.vcf", "MCL_mutsigcv.vcf", "MCL_oncoDriveFM.vcf", "MCL_Sweden.vcf", "MCL_Zhang_Dave_2014.vcf"))
snv <- merge(snv, samples, by.x = "i", by.y = "iteration")

make <- function(x) {
  snv.kmt <- merge(snv, x, by = "filename", all.x = TRUE, all.y = FALSE)
  snv.kmt <- snv.kmt[!duplicated(snv.kmt[1]),]
  snv.kmt$samples <- as.character(snv.kmt$samples)
  snv.kmt$type <- sapply(strsplit(snv.kmt$samples, "_", fixed = FALSE), function(x) (x[1]))
  snv.kmt$mut.type <- as.character(snv.kmt$mut.type)
  snv.kmt$filename <- as.character(snv.kmt$filename) 
  # Combine snvs and indels into one overall count
  #snv.kmt.singles <- snv.kmt[grepl("all", snv.kmt$filename),]
  x <- snv.kmt[!grepl("all", snv.kmt$filename),]
  x <- x[!grepl("DLBCL_TCGA", x$filename),]
  x$type <- factor(x$type)
  x <- x[order(x$type),]
  rownames(x) <- NULL
  x$file1 <- sapply(strsplit(x$filename, "/", fixed = FALSE), function(x) (x[1]))
  x$file2 <- sapply(strsplit(x$filename, "/", fixed = FALSE), function(x) (x[2]))
  x$file3 <- sapply(strsplit(x$filename, "/", fixed = FALSE), function(x) (x[3]))
  x$file4 <- sapply(strsplit(x$filename, "/", fixed = FALSE), function(x) (x[4]))
  x$file7 <- sapply(strsplit(x$filename, "_", fixed = FALSE), function(x) (x[1]))
  x$file8 <- sapply(strsplit(x$filename, "\\.", fixed = FALSE), function(x) (x[1]))
  x$file8 <- sapply(strsplit(x$file8, "_", fixed = FALSE), function(x) (x[4]))
  x$match <- ifelse(x$samples == "BL_ICGC.vcf", paste(x$file1, x$file2, x$file3, sep="/"),
                    ifelse(x$samples == "BL_Love_Dave_2012.vcf", paste(x$file1, x$file2, x$file3, x$file4, sep="/"),
                           ifelse(x$samples == "CLL_ICGC.vcf", x$file7,
                                  ifelse(x$samples == "CLL_Connors.vcf", paste(x$file1, x$file2, sep="/"),
                                         ifelse(x$samples == "DLBCL_Connors.vcf", paste(x$file1, x$file2, x$file3, sep="/"),
                                                ifelse(x$samples == "DLBCL_ICGC.vcf", paste(x$file1, x$file2, sep="/"),
                                                       ifelse(x$samples == "DLBCL_refractory.vcf", x$file8,
                                                              ifelse(x$samples == "DLBCL_TCGA.vcf", x$file2,
                                                                     ifelse(x$samples == "FL_Connors.vcf", x$file1,
                                                                            ifelse(x$samples == "FL_Shah.vcf", x$file1,
                                                                                   ifelse(x$samples == "MCL_Bea_Campo.vcf", paste(x$file1, x$file2, x$file3, sep="/"),
                                                                                          ifelse(x$samples == "MCL_Morin_2015.vcf", paste(x$file1, x$file2, sep="/"),
                                                                                                 ifelse(x$samples == "MCL_mutsigcv.vcf", paste(x$file1, x$file2, sep="/"),
                                                                                                        ifelse(x$samples == "MCL_oncoDriveFM.vcf", paste(x$file1, x$file2, sep="/"),
                                                                                                               ifelse(x$samples == "MCL_Sweden.vcf", paste(x$file1, x$file2, sep="/"),
                                                                                                                      ifelse(x$samples == "MCL_Zhang_Dave_2014.vcf", paste(x$file1, x$file2, sep="/"),
                                                                                                                             NA))))))))))))))))
  
  x$ind <- with(x, ave(as.character(match), match, FUN = seq_along))
  x <- reshape(x, direction = "wide", idvar = "match", timevar = "ind")
  x$total.events <- x$'events.1' + x$'events.2'
  x$mut <- ifelse(grepl("KMT2D|KMT2C|PARP1", x$'gene.1') | grepl("KMT2D|KMT2C|PARP1", x$'gene.2'), "MUT", "WT")
  x$mut[is.na(x$mut)] <- "WT"
  #x$mut <- factor(x$mut, levels=c("WT", "MUT"))
  x$V1 <- paste(x$'type.1', x$mut, sep="_")
  x$V1 <- factor(x$V1, levels = c("BL_WT", "BL_MUT", "CLL_WT", "CLL_MUT", "DLBCL_WT", "DLBCL_MUT", 
                                  "FL_WT", "FL_MUT", "MCL_WT", "MCL_MUT"))
  x$match <- as.character(x$match)
  x$mut <- as.character(x$mut)
  x$V1 <- as.character(x$V1)
  x$total.events <- as.numeric(as.character(x$total.events))
  x <- x[,c(1:3,5,6,8,11:15,17,24,27,30,34:38,46:48)]
  row.names(x) <- NULL
  colnames(x) <- c("match", "indel.file", "i", "indels", "sample.1", "gene.1", "pos.1", "ref.1", "alt.1", "alt2.1", "mut.type.1", "type.1", "snv.file", 
                    "snvs", "gene.2", "pos.2", "ref.2", "alt.2", "alt2.2", "mut.type.2", "total.events", "mut", "type_mut")
  return(x)
}

snv.kmt2d <- make(kmt2d.effects)
row.names(snv.kmt2d) <- NULL

snv.kmt2d$ if (snv.kmt2d$mut.type.1 == "frameshift") {snv.kmt2d$PREDICTION..cutoff..2.5..1 <- "deleterious"
})











snv.kmt2d$PREDICTION <- as.character(snv.kmt2d$PREDICTION)
snv.kmt2d$PREDICTION[is.na(snv.kmt2d$PREDICTION)] <- "WT"
snv.kmt2d[snv.kmt2d$PREDICTION == "",] <- "WT"

snv.kmt2d$PREDICTION..cutoff..2.5..2 <- as.character(snv.kmt2d$PREDICTION)
snv.kmt2d$PREDICTION..cutoff..2.5..2[is.na(snv.kmt2d$PREDICTION)] <- "WT"
snv.kmt2d[snv.kmt2d$PREDICTION..cutoff..2.5..2 == "",] <- "WT"




events.factor.b <- boxplot(total.events~mut, data=snv.kmt2d, plot=0)
events.factor <- boxplot(total.events~PREDICTION..cutoff..2.5..1,data=snv.kmt2d, las = 2,
                         #col=(c("red","blue")),
                         #names=paste(events.factor.b$names, "(n=", events.factor.b$n, ")"),
                         main = "KMT2D status versus number of SNV and indels in genome",
                         #notch = TRUE,
                         par(mar = c(12, 5, 4, 2)+ 0.1) )
events.factor.ylim <- boxplot(total.events~mut,data=snv.ind.comb, las = 2,
                              col=(c("red","blue")),
                              ylim=c(0,35000),
                              main = "KMT2D status versus number of SNV and indels in genome",
                              names=paste(events.factor.b$names, "(n=", events.factor.b$n, ")"),
                              notch = TRUE)
