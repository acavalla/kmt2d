# makes files appropriate for upload to pecan st jude protein paint (requirements and format underneath code)
# The arrange function output, when uploaded to https://pecan.stjude.cloud/proteinpaint, produces an infographic labelled with the alternative (new & mutated) nucleotide
# arrange.amino.acid creates an infographic where mutations are labelled by the previous amino acid, position, and new amino acid.

##load tsv created by gene_muts.R script
transcr <- "ENST00000460911" #EZH2
x <- read.table(paste0("XXX/xxx", transcr, "_mutations.txt"))

gene.name <- "EZH2"
refseq <- "NM_001203247"
arrange.nucleotide <- function(x, gene.name, refseq) {
  x$mut.type <- as.character(x$mut.type)
  x$pos <- x$startpos
  x$alt <- as.character(x$alt)
  x$alt2 <- as.character(x$alt2)
  x$alt <- x$alt2
  x <- as.data.frame(cbind(gene.name, refseq, x$chr, x$pos, x$alt, x$mut.type))
  colnames(x) <- c("gene", "refseq", "chromosome", "start", "aachange", "class")
  x <- as.data.frame(x)
  remove(gene.name)
  remove(refseq)
  return(x)
  }

x <- arrange(x, gene.name, refseq)

gene.name <- "EZH2"
refseq <- "NM_001203247"
code3 <- c("Ala", "Arg", "Asn", "Asp", "Cys", "Glu", "Gln", "Gly", "His", 
           "Ile", "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", 
           "Tyr", "Val")
code1 <- c("A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", 
           "M", "F", "P", "S", "T", "W", "Y", "V")
arrange.amino.acid <- function(x, gene.name, refseq) {
  x$mut.type <- as.character(x$mut.type)
  x$pos <- x$startpos
  x$alt <- as.character(x$alt)
  x$alt2 <- as.character(x$alt2)
  x$alt <- x$alt2
  x$mut.type <- ifelse(x$mut.type == "synonymous", "silent", x$mut.type)
  # a <- ifelse(x$mut.type == "missense", "M",
  #             ifelse(x$mut.type =="exon", "E",
  #                    ifelse(x$mut.type =="frameshift", "F",
  #                           ifelse(x$mut.type =="nonsense", "N",
  #                                  ifelse(x$mut.type =="silent", "S",
  #                                         ifelse(x$mut.type =="proteinDel", "D",
  #                                                ifelse(x$mut.type =="proteinIns", "I",
  #                                                       ifelse(grepl("splice", x$mut.type), "P",
  #                                                              ifelse(x$mut.type =="intron", "Intron",
  #                                                                     ifelse(x$mut.type =="utr_3", "Utr3",
  #                                                                            ifelse(x$mut.type =="utr_5", "Utr5",
  #                                                                                   NA)))))))))))
  

  x$aachange <- gsub("p.", "", x$aachange)
  b <- substr(gsub("[[:digit:]]","",x$aachange), 1,3)
  c <- substr(gsub("[[:digit:]]","",x$aachange), 4,6)
  for (i in 1:length(code3)){
    b <- gsub(code3[i],code1[i],b,ignore.case=TRUE)
    c <- gsub(code3[i],code1[i],c,ignore.case=TRUE)
  }
  x$aachange <- paste0(b,gsub("[^[:digit:]]", "", x$aachange), c)
  x$aachange <- ifelse(grepl("NANANA", x$aachange), as.character(x$alt), as.character(x$aachange))
  
  x <- as.data.frame(cbind(gene.name, refseq, x$chr, x$pos, x$aachange, x$mut.type))
  colnames(x) <- c("gene", "refseq", "chromosome", "start", "aachange", "class")
  x <- as.data.frame(x)
  remove(gene.name)
  remove(refseq)
  return(x)
}

x <- arrange.amino.acid(x, gene.name, refseq)

  
write.table(x, paste0("/projects/acavalla_prj/kmt2d/301067/patient_data/gene_muts/pecan/chop_", gene.name, ".tsv"), sep="\t", row.names = FALSE, quote = FALSE)

# Required columns
# gene: The gene symbol, e.g. "TP53"
# refseq: The RefSeq accession of the gene isoform, e.g. "NM_000546" (but not version name NM_000546.5)
# chromosome: Can be either "chr1" or "1", case-insensitive
# start: Genomic coordinate, 1-based
# aachange: A string describing the amino acid change caused by the mutation, can be any string, will be displayed as-is as the label on the right of the discs
# class: Mutation class of the amino acid change, will allow ProteinPaint to render the mutation using predefined colors. Names of supported classes:
# missense; nonsense (representing stoploss, nonstop, and stopgain); frameshift; silent;
# proteinDel (in-frame deletion); proteinIns (in-frame insertion); splice; splice_region; utr_3; utr_5
# intron; noncoding.

# To indicate amino acid change
# 
# Example: R267W:M
# 
# "R" is reference residue, and "W" is the mutant residue. "267" is the amino acid position of "R", and is used to position this variant on protein. "M" is mutation class, joined by colon with the amino acid change. The class can allow rendering the mutation in a predefined color. Full class list:
#   
#   M ------ MISSENSE
# E ------ EXON
# F ------ FRAMESHIFT
# N ------ NONSENSE  (covers nonstop, stoploss, stopgain)
# S ------ SILENT
# D ------ PROTEINDEL  (in-frame insertion)
# I ------ PROTEININS    (in-frame deletion)
# P ------ SPLICE_REGION  (within 10 bp range from splice site)
# L ------ SPLICE
# Intron - INTRON
# Utr3 --- UTR_3
# Utr5 --- UTR_5
# noncoding --- (either intergenic or intronic)
