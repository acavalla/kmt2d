### loads in list of KMT2D mutational status and list of all samples and their snv and indel counts, both from MAF files. Finds 
## unique string identifiers and uses them to match snv and indel files per sample, reshaping the table and 
## aggregating the two files to give a total somatic event count and overall kmt2d mutational status. Also
## merged in vcf files, to indicate which samples  we can analyse mutational signatures from.
## (deconstructSigs requires vcfs, and where there is no matching maf, the kmt2d stat cannot be determined)

match <- function(x){
  x$file <- as.character(x$file)
  x$file1 <- sapply(strsplit(x$file, "/", fixed = FALSE), function(x) (x[3]))
  x$file2 <- sapply(strsplit(x$file, "/", fixed = FALSE), function(x) (x[4]))
  x$file3 <- sapply(strsplit(x$file, "/", fixed = FALSE), function(x) (x[5]))
  x$file4 <- sapply(strsplit(x$file1, "_", fixed = FALSE), function(x) (x[1]))
  x$file5 <- sapply(strsplit(x$file2, "_", fixed = FALSE), function(x) (x[4]))
  x$file5 <- sapply(strsplit(x$file5, "\\.", fixed = FALSE), function(x) (x[1]))
  x$file6 <- sapply(strsplit(x$file1, "_", fixed = FALSE), function(x) (x[4]))
  x$file6 <- sapply(strsplit(x$file6, "\\.", fixed = FALSE), function(x) (x[1]))
  x$file7 <- sapply(strsplit(x$file2, "_", fixed = FALSE), function(x) (x[3]))
  x$match <- ifelse(grepl("BL_ICGC", x$file), x$file3,
                    ifelse(grepl("BL_Love", x$file), x$file3,
                           ifelse(grepl("CLL_Connors", x$file), x$file1,
                                  ifelse(grepl("CLL_ICGC", x$file), x$file4,
                                         ifelse(grepl("batch7", x$file), x$file2,
                                                ifelse(grepl("first30.*TASK_STRELKA", x$file), x$file7,
                                                       ifelse(grepl("first30.*VCF2MAF", x$file), x$file5,
                                                              ifelse(grepl("DLBCL_ICGC", x$file), x$file1,
                                                                     ifelse(grepl("DLBCL_TCGA", x$file), x$file2,
                                                                            ifelse(grepl("DLBCL_refrac", x$file), x$file6,
                                                                                   ifelse(grepl("FL_Connors", x$file), x$file1,
                                                                                          ifelse(grepl("FL_Shah", x$file), x$file1,
                                                                                                 ifelse(grepl("FL_Fitz", x$file), x$file1,
                                                                                                        ifelse(grepl("MCL", x$file), x$file3,
                                                                                                               NA))))))))))))))
  x <- x[, c(1,2,3,ncol(x))]
  x <- cbind(x[,1], x[,4], x[,2], x[,3])
  x <- as.data.frame(x)
  names(x) <- c("file", "match", "events", "mut")
  return(x)
}

kmt2d <- read.table("/projects/acavalla_prj/kmt2d/301067/patient_data/gene_muts/kmt2d_mafmuts.tsv")
snv <- read.table("/projects/acavalla_prj/kmt2d/301067/patient_data/all_snvs/all_snvs.maf.txt", sep="\t")

snv.x <- merge(snv, kmt2d, all.x = TRUE)
snv.x <- snv.x[!duplicated(snv.x$file),]
snv.x$mut <- ifelse(is.na(snv.x$gene), "WT", "MUT")
snv.x <- snv.x[,c(1,2, ncol(snv.x))]
x <- snv.x

y <- read.table("/projects/acavalla_prj/kmt2d/301067/patient_data/files/filenames.vcf.txt")
names(y) <- "file"
y <- match(y)
y <- y[,1:2]
names(y) <- c("file", "match")

x <- match(x)
x$snv <- ifelse(grepl("snv|SNV", x$file), "snv", "indel")
x <- x[order(x$snv),]
x$matchsnv <- paste(x$match, x$snv)

y$snv <- ifelse(grepl("snv|SNV", y$file), "snv", "indel")
y <- y[order(y$snv),]
y$matchsnv <- paste(y$match, y$snv)


x <- merge(x, y, by = "matchsnv")
x <- x[!grepl("Fitz", x$file.x),]

x$ind <- with(x, ave(as.character(match.x), match.x, FUN = seq_along))
x <- reshape(x, direction = "wide", idvar = "match.x", timevar = "ind")
x <- x[order(x$file.x.1),]
rownames(x) <- 1:nrow(x)

x$mut <- ifelse(x$mut.1=="MUT" | x$mut.2=="MUT", "MUT", "WT")
x$total.somatic <- as.numeric(as.character(x$events.1)) + as.numeric(as.character(x$events.2))

x <- as.data.frame(cbind("match" = as.character(x$match.x), "match.snv" = x$matchsnv.1,
                         "mafindel" = as.character(x$file.x.1), "vcfindel" = as.character(x$file.y.1), "indels" = as.numeric(as.character(x$events.1)),
                         "mafsnv" = as.character(x$file.x.2), "vcfsnv" = as.character(x$file.y.2), "snvs" = as.numeric(as.character(x$events.2)),
                         "total.somatic" = x$total.somatic, "mut" = x$mut))
write.table(x, "/projects/acavalla_prj/kmt2d/301067/patient_data/summarybysample.txt", row.names = F, quote = F, col.names=T, sep = "\t")
