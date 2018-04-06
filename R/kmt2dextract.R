## extracting mutations of a transcript of interest to enter into eg lolliplot.R


library(ggplot2)
##EFF= Effect ( Effect_Impact | Functional_Class | Codon_Change | Amino_Acid_Change| Amino_Acid_Length | Gene_Name | Transcript_BioType | Gene_Coding | Transcript_ID | Exon_Rank  | Genotype_Number [ | ERRORS | WARNINGS ] )

process.vcf <- function(x){
  x <- tryCatch(fread(paste("grep", transcr, file.names.celllines[i]), header=FALSE, sep="\t", stringsAsFactors=FALSE), error=function(e) NULL) 
  #file <- file[file$V1 == chr,]
  if(!is.null(x)){
    x <- x[,c(-3, -7, -9)]
    x <- as.data.frame(as.matrix(x))
    x$file <- file.names.celllines[i]
    x$i <- i
    x$V8 <- as.character(x$V8)
    aa <- strsplit(x$V8, ",", fixed = FALSE, perl = FALSE, useBytes = FALSE)
    aa <- lapply(aa, function(ch) grep(transcr, ch, value = TRUE))
    aa <- strsplit(as.character(aa), "=", fixed = FALSE, perl = FALSE, useBytes = FALSE)
    aa <- lapply(aa, function(ch) grep(transcr, ch, value = TRUE))
    aa <- strsplit(as.character(aa), "\\|", fixed = FALSE, perl = FALSE, useBytes = FALSE)
    for(j in 1:length(aa)){
      x[j,6] <- sapply(strsplit(as.character(aa[[j]][1]), "\\(", fixed = FALSE, perl = FALSE, useBytes = FALSE), function(x) x[1])
      x[j,10] <- sapply(strsplit(as.character(aa[[j]][1]), "\\(", fixed = FALSE, perl = FALSE, useBytes = FALSE), function(x) x[2])
      x[j,11] <- aa[[j]][2]
      x[j,12] <- aa[[j]][3]
      x[j,13] <- aa[[j]][4]
    }
    
    x$V1 <- gsub("^chr", "", x$V1)
    #x$V1 <- as.character(x$V1)
    x <- x[!grepl("PT1092",x$file),]
    print(i)
  } else {
    print(i)
  }
  return(x)
}


file.names.celllines <- c("/projects/analysis/analysis19/A12437/merge_bwa-0.5.7/100nt/hg19a/A12437_3_lanes_dupsFlagged.varFilter.eff.snvs.vcf",
                "/projects/analysis/analysis19/A12438/merge_bwa-0.5.7/100nt/hg19a/A12438_3_lanes_dupsFlagged.varFilter.eff.snvs.vcf",
                "/projects/analysis/analysis26/A51671/HNV23CCXX_2/A51671/150nt/hg19a/bwa-0.5.7/A51671_1_lane_dupsFlagged.varFilter.eff.vcf",
                "/projects/analysis/analysis25/A51679/HNV23CCXX_4/A51679/150nt/hg19a/bwa-0.5.7/A51679_1_lane_dupsFlagged.varFilter.eff.vcf",
                "/projects/analysis/analysis25/A51620/HNTHYCCXX_2/A51620/150nt/hg19a/bwa-0.5.7/A51620_1_lane_dupsFlagged.varFilter.eff.vcf",
                "/projects/analysis/analysis19/A12439/merge_bwa-0.5.7/100nt/hg19a/A12439_3_lanes_dupsFlagged.varFilter.eff.snvs.vcf",
                "/projects/analysis/analysis19/A12440/merge_bwa-0.5.7/100nt/hg19a/A12440_3_lanes_dupsFlagged.varFilter.eff.snvs.vcf",
                "/projects/analysis/analysis19/A12441/merge_bwa-0.5.7/100nt/hg19a/A12441_3_lanes_dupsFlagged.varFilter.eff.snvs.vcf",
                "/projects/analysis/analysis19/A12442/merge_bwa-0.5.7/100nt/hg19a/A12442_3_lanes_dupsFlagged.varFilter.eff.snvs.vcf",
                "/projects/analysis/analysis19/A12443/merge_bwa-0.5.7/100nt/hg19a/A12443_3_lanes_dupsFlagged.varFilter.eff.snvs.vcf",
                "/projects/analysis/analysis19/A12449/merge_bwa-0.5.7/100nt/hg19a/A12449_3_lanes_dupsFlagged.varFilter.eff.snvs.vcf",
                "/projects/analysis/analysis19/A12444/merge_bwa-0.5.7/100nt/hg19a/A12444_3_lanes_dupsFlagged.varFilter.eff.snvs.vcf",
                "/projects/analysis/analysis25/A51629/HNTHLCCXX_2/A51629/150nt/hg19a/bwa-0.5.7/A51629_1_lane_dupsFlagged.varFilter.eff.vcf",
                "/projects/analysis/analysis19/A12445/merge_bwa-0.5.7/100nt/hg19a/A12445_3_lanes_dupsFlagged.varFilter.eff.snvs.vcf",
                "/projects/analysis/analysis19/A10088/merge_bwa-0.5.7/100nt/hg19a/A10088_2_lanes_dupsFlagged.varFilter.eff.snvs.vcf",
                "/projects/analysis/analysis25/A51638/HNTHLCCXX_6/A51638/150nt/hg19a/bwa-0.5.7/A51638_1_lane_dupsFlagged.varFilter.eff.vcf",
                "/projects/analysis/analysis25/A51647/HNJJ5CCXX_2/A51647/150nt/hg19a/bwa-0.5.7/A51647_1_lane_dupsFlagged.varFilter.eff.vcf",
                "/projects/analysis/analysis19/A12446/merge_bwa-0.5.7/100nt/hg19a/A12446_3_lanes_dupsFlagged.varFilter.eff.snvs.vcf",
                "/projects/analysis/analysis19/A12447/merge_bwa-0.5.7/100nt/hg19a/A12447_3_lanes_dupsFlagged.varFilter.eff.snvs.vcf",
                "/projects/analysis/analysis25/A51656/HNJJ5CCXX_5/A51656/150nt/hg19a/bwa-0.5.7/A51656_1_lane_dupsFlagged.varFilter.eff.vcf",
                "/projects/analysis/analysis25/A51664/HNJJ5CCXX_8/A51664/150nt/hg19a/bwa-0.5.7/A51664_1_lane_dupsFlagged.varFilter.eff.vcf",
                "/projects/analysis/analysis21/A14217/merge_bwa-0.5.7/100nt/hg19a/A14217_3_lanes_dupsFlagged.varFilter.eff.snvs.vcf",
                "/projects/analysis/analysis26/A51672/HTTVJCCXX_8/A51672/150nt/hg19a/bwa-0.5.7/A51672_1_lane_dupsFlagged.varFilter.eff.vcf",
                "/projects/analysis/analysis25/A51680/HNT57CCXX_1/A51680/150nt/hg19a/bwa-0.5.7/A51680_1_lane_dupsFlagged.varFilter.eff.vcf",
                "/projects/analysis/analysis25/A51621/HNTHYCCXX_3/A51621/150nt/hg19a/bwa-0.5.7/A51621_1_lane_dupsFlagged.varFilter.eff.vcf",
                "/projects/analysis/analysis25/A51630/HNTHLCCXX_3/A51630/150nt/hg19a/bwa-0.5.7/A51630_1_lane_dupsFlagged.varFilter.eff.vcf",
                "/projects/analysis/analysis25/A51639/HNTHLCCXX_7/A51639/150nt/hg19a/bwa-0.5.7/A51639_1_lane_dupsFlagged.varFilter.eff.vcf",
                "/projects/analysis/analysis21/A51648/merge_bwa-0.5.7/125nt/hg19a/A51648_2_lanes_dupsFlagged.varFilter.eff.snvs.vcf")

transcr <- "ENST00000301067" ##transcript of your choice 

datalist <- list()
for(i in 1:length(file.names.celllines)){
  file <- process.vcf(file.names.celllines[i])
  datalist[[i]] <- file # add it to your list
}

file <- do.call(rbind, datalist)

colnames(file) <- c("filename", "i", "CHROM", "POS","ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "P8")

metadata <- data.frame("cell.line" = c("DB", "DoHH2", "Farage", "HBL-1", "HT", "Karpas422", "MD903", "Nu-DHL-1", "Nu-DUL-1", "OCI-Ly1", "OCI-Ly3", "OCI-Ly7", "OCI-Ly10", "OCI-Ly19", "Pfeiffer", "Su-DHL-4",
                                       "Su-DHL-5", "Su-DHL-6", "Su-DHL-9", "SU-DHL-10", "Toledo", "WSU_DLCL2", "WSU-NHL", "GM12878", "GM12891", "GM17010", "GM18507",  "HEK293_D320"),
                       "i" = 1:28, c.type =  c("GCB-DLBCL", "GCB-DLBCL", "GCB-DLBCL", "GCB-DLBCL", "DLBCL", "GCB-DLBCL", "ABC-DLBCL", "GCB-DLBCL", "ABC-DLBCL", "GCB-DLBCL", "ABC-DLBCL", "GCB-DLBCL", 
                                               "ABC-DLBCL", "GCB-DLBCL", "GCB-DLBCL", "GCB-DLBCL", "GCB-DLBCL", "GCB-DLBCL", "ABC-DLBCL", "GCB-DLBCL", "GCB-DLBCL", "GCB-DLBCL", "GCB-DLBCL",
                                               "B-Lymphoblastoid", "B-Lymphoblastoid", "B-Lymphoblastoid", "B-Lymphoblastoid", "Human Embryonic Kidney - MLL2 KO"), stringsAsFactors = FALSE)
file <- merge(metadata, file, by = "i", all.y=TRUE)



write.table(file, file="/projects/acavalla_prj/kmt2d/301067/cell_lines/kmt2d.tsv")

snvlist <- list()
for(i in 1:length(file.names.celllines)){
  x <- tryCatch(fread(paste("wc -l", file.names.celllines[i]), header=FALSE, sep="\t", stringsAsFactors=FALSE), error=function(e) NULL) 
  x$file <- sapply(strsplit(x$V1, " ", fixed = FALSE), function(x) (x[2]))
  x$V1 <- sapply(strsplit(x$V1, " ", fixed = FALSE), function(x) (as.numeric(x[1])))
  x$V1 <- x$V1-1
  snvlist[[i]] <- x # add it to your list
  print(i)
}
snv <- do.call(rbind, snvlist)
write.table(snv, "/projects/acavalla_prj/kmt2d/301067/patient_data/all_snvs/all_snvs.maf.txt", sep="\t")
