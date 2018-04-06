library(Gviz)
library(rtracklayer)
library(trackViewer) ##lolliplot
library(VariantAnnotation)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(data.table) ##fread
library(org.Hs.eg.db)

##Read kmt2dextract product, select single transcript
data <- read.table("/projects/acavalla_prj/kmt2d/ALLmuts.tsv")
MLL2.1 <- with(data, data[ grepl( "ENST00000301067", data$INFO),])

##split info column, remove other transcript annotations, remove EFF= at beginning
spl.info.col <- function(x) {
  colnames(x) <- c("CHR","POS", "ID", "REF", "ALT", "QUAL", "V6", "INFO", "FORMAT", "P8")
  x <- with(x, x[ -grep ( "INDEL", x$INFO),])
  a <- as.data.frame(gsub(".*;","",x$INFO))
  colnames(a) <- c("b")
  a$b <- as.character(a$b)
  aa <- strsplit(a$b, ",", fixed = FALSE, perl = FALSE, useBytes = FALSE)
  aa <- lapply(aa, function(ch) grep("ENST00000301067", ch, value = TRUE))
  aa <- gsub("^EFF=","",aa)
  x$INFO <- as.character(x$INFO)
  x$DP <- sapply(strsplit(x$INFO, ";", fixed = FALSE), function(x) (x[1]))
  x$DP <- gsub("DP=","",x$DP)
  x$AF1 <- sapply(strsplit(x$INFO, ";", fixed = FALSE), function(x) (x[2]))
  x$AF1 <- gsub("AF1=","",x$AF1)
  x$AC1 <- sapply(strsplit(x$INFO, ";", fixed = FALSE), function(x) (x[3]))
  x$AC1 <- gsub("AC1=","",x$AC1)
  x$DP4 <- sapply(strsplit(x$INFO, ";", fixed = FALSE), function(x) (x[4]))
  x$DP4 <- gsub("DP4=","",x$DP4)
  x$MQ <- sapply(strsplit(x$INFO, ";", fixed = FALSE), function(x) (x[5]))
  x$MQ <- gsub("MQ=","",x$MQ)
  x$FQ <- sapply(strsplit(x$INFO, ";", fixed = FALSE), function(x) (x[6]))
  x$FQ <- gsub("FQ=","",x$FQ)
  x$PV4 <- sapply(strsplit(x$INFO, ";", fixed = FALSE), function(x) (x[7]))
  x$PV4 <- gsub("PV4=","",x$PV4)
  x$INFO <- as.character(aa)
  x$LOCATION <- sapply(strsplit(x$INFO, "\\(", fixed = FALSE), function(x) (x[1]))
  x$INFO <- sapply(strsplit(x$INFO, "\\(", fixed = FALSE), function(x) (x[2]))
  x <- cbind(x[1:6], x[10:18], x[7:9])
  return(x)
}
#MLL2.1$PV4[which (MLL2.1$PV4 ^"EFF") ] <- ""

##"Features" - ref exon coords + largest and smallest labelled SNVs.
grch37 <- fread(paste('grep MLL2-001', "/projects/acavalla_prj/kmt2d/Homo_sapiens.GRCh37.69.gtf"), sep="\t")
grch37 <- subset(grch37, grch37$V3=="exon")
grch37 <- grch37[order(grch37$V4),] 
grch37$V10 <- grch37$V5-grch37$V4
exons <- c(min(MLL2.1$POS), grch37$V4, max(MLL2.1$POS))
lengths <- c(1, grch37$V10, 1)

###colour lollipops with reference to type of mut
colour.lol <- function(x) {
  x$location.colour <- ifelse(x$LOCATION=="DOWNSTREAM", 1,
                              ifelse(x$LOCATION=="INTRON", 2,
                                     ifelse(x$LOCATION=="NON_SYNONYMOUS_CODING", 3,
                                            ifelse(x$LOCATION=="STOP_GAINED", 4,
                                                   ifelse(x$LOCATION=="SYNONYMOUS_CODING", 5,
                                                          ifelse(x$LOCATION=="UPSTREAM", 6,
                                                                 ifelse(x$LOCATION=="UTR_3_PRIME", 7,
                                                                        ifelse(x$LOCATION=="FRAME_SHIFT", 8,
                                                                               NA  ))))))))
  return(x)
}

## remove multiple labels for repeat SNPs, make variable to be called later
a <- as.data.frame(MLL2.1$POS)
b <- as.data.frame(table(MLL2.1$POS))
b$Mid <- ceiling(b$Freq/2)

a$Index = 0
for(j in 1:length(MLL2.1$POS)){
  snpinds <- length(a$`MLL2.1$POS`==a$`MLL2.1$POS`[j])
  a$Index[snpinds] <- 1:length(snpinds)
  }

a$no<-rownames(a)
a$no <- as.numeric(a$no)
ab <- merge(a, b, by.x="MLL2.1$POS", by.y="Var1")
ab <- ab[order(ab$no),]

## lay down base sample gr and features
sample.gr <- GRanges("chr12", IRanges(c(MLL2.1$POS), width=1, names=NULL))
features <- GRanges("chr12", IRanges(exons, 
                                     width=lengths,
                                     names=NULL))


# customising of features and sample.gr
features$height <- 0.04
names(sample.gr)[ab$Mid!=ab$Index] <- " "
sample.gr$color <- MLL2.1$location.colour
##axis fix
xaxis <- c(min(MLL2.1$POS), 49410000, 49415000, 49420000, 49425000, 49430000, 49435000, 49440000, 49445000, 49450000, max(MLL2.1$POS)) ## define the position
yaxis <- c(40, 50, 100, 150, 200, 225, 350)
names(yaxis) <- yaxis
names(yaxis)[7] <- ""
sample.gr$score <- MLL2.1$QUAL #Lack of this with y-axis label fs up the diagram

##make legend
cols <- c("DOWNSTREAM",
          "FRAME_SHIFT",
          "INTRON",
          "NON_SYNONYMOUS_CODING",
          "STOP_GAINED",
          "SYNONYMOUS_CODING",
          "UPSTREAM",
          "UTR_3_PRIME")
sample.gr$score <- MLL2.1$i
legend <- 1:8 ## legend fill color
names(legend) <- cols
legend <- list(labels=cols, 
               col=palette()[1:8], 
               fill=palette()[legend])

#plot & title
lolliplot(sample.gr, features, xaxis = xaxis, ylab=" ", yaxis = FALSE, cex=.4, legend=legend)
grid.text("KMT2D mutations from WGS", x=.5, y=.5,
          gp=gpar(cex=1.5, fontface="bold"))

##for running individual plots for a sequence of filenames
#lolliplot(sample.gr, features, xaxis = xaxis, yaxis = yaxis, ylab=" ")
#grid.text(MLL2.1$cell.types[MLL2.1[,1]==i], x=.5, y=.98, just="top", 
#          gp=gpar(cex=5, fontface="bold"))


##add track of protein domains - incomplete
uni <- read.csv("/projects/acavalla_prj/kmt2d/uniprot.csv")
uni$relstart <- uni$Start + exons[2]
features_domain <- GRanges("chr12", IRanges(uni$relstart, 
                                            width=uni$Length_bp,
                                            names=NULL))
features.mul <- rep(features_domain, 2)
lolliplot(sample.gr, features.mul)

##attempt with results migrated to grch38 (from ensembl tool)
#new VCF SNPs
datamig <- read.table("/projects/acavalla_prj/kmt2d/301067/output_301067_to_migrate.vcf")
datamig <- cbind(MLL2.1$i, MLL2.1$cell.types, datamig)


#new exons
grch38 <- fread("/projects/acavalla_prj/kmt2d/ExonsSpreadsheet-Homo_sapiens_Transcript_Exons_ENST00000301067grch38.csv")
toDelete <- seq(2, nrow(grch38), 2)
grch38 <- grch38[ toDelete ,]
grch38$Start <- as.numeric(gsub(",","",grch38$Start))
grch38$End <- as.numeric(gsub(",","",grch38$End))
grch38$Length <- as.numeric(gsub(",","",grch38$Length))
grch38 <- grch38[order(grch38$Start),]
exons38 <- c(min(datamig$POS), grch38$End)
lengths <- c(1, grch38$Length)
sample.gr <- GRanges("chr12", IRanges(c(datamig$POS), width=1, names=NULL))
features <- GRanges("chr12", IRanges(exons38, 
                                     width=lengths,
                                     names=NULL))


datamig <- spl.info.col(datamig)
datamig <- colour.lol(datamig)

features$height <- 0.04
sample.gr$color <- datamig$location.colour
lolliplot(sample.gr, features)
