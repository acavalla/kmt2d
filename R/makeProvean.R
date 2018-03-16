# PROVEAN Genome Variants - Input Format
# 
# Input: a list of genome variants
# 
# Each genome variant is represented in comma-separated values as the following:
#   <chromosome>,<position>,<reference allele>,<variant allele>,<comment(optional)>
#   
#   Chromosome:
#   Chromosome name (1, 2, 3, ..., 22, X, or Y).
# Position:
#   Reference position, with the 1st base having position 1.
# Reference allele:
#   One or more nucleotides in the reference genome. The value in the "Position" field refers to the position of the first base in this field. At least one base is required.
# Variant allele:
#   One or more nucleotides for non-reference allele. The bases in the "Reference allele" field are replaced by the bases in this field. At least one base is required.
# Example 1 (SNP): 1,100382265,C,G,some comments
# Example 2 (Deletion): 7,117199646,CTTT,C
# Example 3 (Insertion): 1,43217995,G,GCCA
# Example 4 (MNP): 10,102762471,AG,CC

##input: file created by gene_muts.R

transcr <- "ENST00000301067" ##transcript of your choice
x <- read.table(paste0("XXX/xxx", transcr, "_mutations.txt"))
x <- x[x$assembly=="GRCh37",]

x$alt <- ifelse(as.character(x$ref) == as.character(x$alt), as.character(x$alt), as.character(x$alt2))

x <- cbind(x$chr, x$pos, x$ref, x$alt)
write.table(x, paste0("XXX/xxx, transcr, ".csv"), sep=",",  col.names=FALSE, row.names = FALSE, quote=FALSE)
