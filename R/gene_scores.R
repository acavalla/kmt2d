#Creating a scatterplot of likelihood of genes being comutated based on Provean scores for deleterious mutations. Currently only including scores given
#May later assign scores to nonsense or frameshift mutations.
#Use cell lines' mutation extractor to obtain mutations from maf & vcf sample files, then makeProvean to make the csv for input into provean's online tool, then 
#loadProvean 'clean' function to deal with cleaning the output then extract functionality from make function to align snv & indel files from the same sample
#then run scatterplot!

#KMT2D vs ASH2L vs RBBP5 as these have been shown to affect the functionality & enzyme efficiency of the KMT2 family complexes (NRC, Yali & Dou)

#Gene mutational status

process.maf <- function(x){
  x <- tryCatch(fread(paste("grep", transcr, file.names.maf[i]), header=FALSE, sep="\t", stringsAsFactors=FALSE), error=function(e) NULL) 
  if(!is.null(x)){
    x <- as.data.frame(as.matrix(x))
    x <- x[,which(unlist(lapply(x, function(x)!all(is.na(x))))),]
    x$file <- file.names.maf[i]
    x$i <- i
    y <- x[!grepl(transcr, x$V38),]
    x <- x[grepl(transcr, x$V38),]
    if(nrow(y) != 0){
      y$V46 <- as.character(y$V46)
      aa <- strsplit(y$V46, ";", fixed = FALSE, perl = FALSE, useBytes = FALSE)
      aa <- lapply(aa, function(ch) grep(transcr, ch, value = TRUE))
      aa <- strsplit(as.character(aa), ",", fixed = FALSE, perl = FALSE, useBytes = FALSE)
      y$V1 <- aa[[1]][1]
      y$V9 <- aa[[1]][2]
      y$V38 <- aa[[1]][4]
      x <- rbind(y, x)
      x <- x[order(x$file, x$V6),] 
    }
    
    #new col of mutation type
    x$mut.type <- ifelse(grepl("Frame_Shift_Del", x$V9), "frameshift",
                         ifelse(grepl("Frame_Shift_Ins", x$V9), "frameshift",
                                ifelse(grepl("Intron", x$V9), "intron",
                                       ifelse(grepl("Missense_Mutation", x$V9), "missense",
                                              ifelse(grepl("Nonsense_Mutation",x$V9), "nonsense",
                                                     ifelse(grepl("In_Frame_Del", x$V9), "proteinDel",
                                                            ifelse(grepl("In_Frame_Ins", x$V9), "proteinIns",
                                                                   ifelse(grepl("Silent", x$V9), "silent",
                                                                          ifelse(grepl("Splice", x$V9), "splice",
                                                                                 ifelse(grepl("3'UTR", x$V9), "utr_3",
                                                                                        ifelse(grepl("3'Flank", x$V9), "utr_3",
                                                                                               ifelse(grepl("5'Flank", x$V9), "utr_5",
                                                                                                      ifelse(grepl("5'UTR", x$V9), "utr_5",
                                                                                                             x$V9)))))))))))))
    
    x <- cbind(x$file, x[1], x[4:7], x[11:13], x$mut.type, x$i)
    x$V11 <- as.factor(revalue(warn_missing=F, x$V11, c("TRUE" = "T")))
    x$V12 <- as.factor(revalue(warn_missing=F, x$V12, c("TRUE" = "T")))
    x$V13 <- as.factor(revalue(warn_missing=F, x$V13, c("TRUE" = "T")))
    x$V5 <- gsub("^chr", "", x$V5)
    #x$V1 <- as.character(x$V1)
    x <- x[!grepl("PT1092",x$`x$file`),]
    print(i)
  } else {
    print(i)
  }
  return(x)
}

##Read in the filenames
file.names <- read.table("/projects/acavalla_prj/kmt2d/301067/patient_data/mutfilenames.txt", sep = "\n", fill = TRUE, stringsAsFactors = FALSE)
file.names <- file.names$V1


transcr <- "ENST00000301067" ##transcript of your choice 

datalist <- list()
setwd("/projects/rmorin/projects/nhl_meta_analysis/results/strelka_pipeline-1.0.0/")
for(i in 1:length(file.names)){
  x <- process.maf(x)
  datalist[[i]] <- x # add it to your list
}

datalist <- datalist[-which(sapply(datalist, is.null))]

file <- do.call(rbind, datalist)
kmt2d <- file



#ASH2L mutational status
transcr <- "ENST00000343823" ##transcript of your choice 
datalist = list()
setwd("/projects/rmorin/projects/nhl_meta_analysis/results/strelka_pipeline-1.0.0/")
for(i in 1:length(file.names)){
  x <- process.maf(x)
  datalist[[i]] <- x # add it to your list
}
datalist <- datalist[-which(sapply(datalist, is.null))]

file <- do.call(rbind, datalist)

ash2l <- file

#RBBP5 mutational status
transcr <- "ENST00000264515" ##transcript of your choice 
datalist = list()

setwd("/projects/rmorin/projects/nhl_meta_analysis/results/strelka_pipeline-1.0.0/")
for(i in 1:length(file.names)){
  x <- process.maf(x)
  datalist[[i]] <- x # add it to your list
}

datalist <- datalist[-which(sapply(datalist, is.null))]
file <- do.call(rbind, datalist)

rbbp5 <- file

write.table(kmt2d, "/projects/acavalla_prj/kmt2d/301067/patient_data/kmt2d_muts.tsv", sep=",",  col.names=FALSE, row.names = FALSE, quote=FALSE)
#### make PROVEAN input csvs 
x <- ash2l
x$V12 <- ifelse(x$V13 != "-", as.character(x$V13), as.character(x$V12))
x.prov <- cbind(x[4], x[5], x[7], x[8])
write.table(x.prov, "/projects/acavalla_prj/kmt2d/301067/patient_data/provean/ash2l.csv", sep=",",  col.names=FALSE, row.names = FALSE, quote=FALSE)

x <- kmt2d
x$V12 <- ifelse(x$V13 != "-", as.character(x$V13), as.character(x$V12))
x.prov <- cbind(x[4], x[5], x[7], x[8])
write.table(x.prov, "/projects/acavalla_prj/kmt2d/301067/patient_data/provean/kmt2d.gs.csv", sep=",",  col.names=FALSE, row.names = FALSE, quote=FALSE)

x <- rbbp5
x$V12 <- ifelse(x$V13 != "-", as.character(x$V13), as.character(x$V12))
x.prov <- cbind(x[4], x[5], x[7], x[8])
write.table(x.prov, "/projects/acavalla_prj/kmt2d/301067/patient_data/provean/rbbp5.csv", sep=",",  col.names=FALSE, row.names = FALSE, quote=FALSE)




####load effects
effect.files <- c("/projects/acavalla_prj/kmt2d/301067/patient_data/provean/ash2l.tsv", "/projects/acavalla_prj/kmt2d/301067/patient_data/provean/kmt2dgs.tsv",
                  "/projects/acavalla_prj/kmt2d/301067/patient_data/provean/rbbp5.tsv")
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

ash2l.effects <- cbind(ash2l, effect.list[1])
kmt2d.effects <- cbind(kmt2d, effect.list[2])
rbbp5.effects <- cbind(rbbp5, effect.list[3])

x <- kmt2d.effects
clean <- function(x){
  x <- x[x$PREDICTION != "Neutral",]
  x <- x[x$`x$mut.type` != "utr_3",]
  x <- x[x$`x$mut.type` != "utr_5",]
  x <- x[x$`x$mut.type` != "intron",]
  x <- x[x$`x$mut.type` != "downstream_variant",]
  x <- x[x$`x$mut.type` != "upstream_variant",]
  x <- x[!is.na(x$`x$file`),]
  return(x)
}

x <- clean(x)
x <- x[complete.cases(x),]
kmt2d.effects <- clean(kmt2d.effects)
ash2l.effects <- clean(ash2l.effects)
rbbp5.effects <- clean(rbbp5.effects)

all <- rbind(kmt2d.effects, ash2l.effects, rbbp5.effects)


all$ind <- with(all, ave(as.character(`x$file`), `x$file`, FUN = seq_along))
colnames(all)[1] <- "file" 
all <- reshape(all, direction = "wide", idvar = "file", timevar = "ind")



as.vector(all[1])
str(as.vector(all$`x$file`))


######There are further details in colV46 of the maf files of muts. This unfinished script extracts them (i can't get lapply to work with the if statement)

transcr <- "ENST00000301067" ##transcript of your choice 
datalist = list()
file.names <- read.table("/projects/acavalla_prj/kmt2d/301067/patient_data/mutfilenames.txt", sep = "\n", fill = TRUE, stringsAsFactors = FALSE)
file.names <- file.names$V1
setwd("/projects/rmorin/projects/nhl_meta_analysis/results/strelka_pipeline-1.0.0/")
for(i in 1:length(file.names)){
  x <- tryCatch(fread(paste("grep", transcr, file.names[714]), header=FALSE, sep="\t", stringsAsFactors=FALSE), error=function(e) NULL) 
  #file <- file[file$V1 == chr,]
    if(!is.null(x)){
      x <- as.data.frame(as.matrix(x))
      x <- x[,which(unlist(lapply(x, function(x)!all(is.na(x))))),]
      x$file <- file.names[i]
      x$temp <- x$V46
      tmp <- sapply(x, function(x)
        if(x[1] != "KMT2D"){
          a <- as.data.frame(x[ncol(x)])
          colnames(a) <- c("b")
          a$b <- as.character(a$b)
          aa <- strsplit(a$b, ";", fixed = FALSE, perl = FALSE, useBytes = FALSE)
          x$tmp <- as.character(lapply(aa, function(ch) grep("ENST00000301067", ch, value = TRUE)))
        } else {
          x$tmp <- "PH"
        }
        )
      ##lapply(dflist, function(x) ifelse(nrow(dflist[[x]])%%2==!0, "ODD", dflist[x]))
      x$V46 <- apply(x, 1, function(x){
        if(temp != "PH"){
            x[21] <- temp
          }
      })
      
      if(temp != "PH"){
        x$V1 <- apply(x, 1, function(x) {
          x[1] <- strsplit(temp, ",", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][1]
        })
      }
      if(temp != "PH"){
        x$V9 <- apply(x, 1, function(x) {
        x[9] <- strsplit(temp, ",", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][2]
        })
      }
      #new col of mutation type
      x$mut.type <- ifelse(grepl("Frame_Shift_Del", x$V9), "frameshift",
                           ifelse(grepl("Frame_Shift_Ins", x$V9), "frameshift",
                                  ifelse(grepl("Intron", x$V9), "intron",
                                         ifelse(grepl("Missense_Mutation", x$V9), "missense",
                                                ifelse(grepl("Nonsense_Mutation",x$V9), "nonsense",
                                                       ifelse(grepl("In_Frame_Del", x$V9), "proteinDel",
                                                              ifelse(grepl("In_Frame_Ins", x$V9), "proteinIns",
                                                                     ifelse(grepl("Silent", x$V9), "silent",
                                                                            ifelse(grepl("Splice_Site", x$V9), "splice",
                                                                                   ifelse(grepl("3'UTR", x$V9), "utr_3",
                                                                                          ifelse(grepl("3'Flank", x$V9), "utr_3",
                                                                                                 ifelse(grepl("5'Flank", x$V9), "utr_5",
                                                                                                        ifelse(grepl("5'UTR", x$V9), "utr_5",
                                                                                                               NA)))))))))))))
      x <- cbind(x$file, x[1], x[4:7], x[11:13], x$mut.type)
      x$i <- i
      x$V11 <- as.factor(revalue(x$V11, c("TRUE" = "T")))
      x$V12 <- as.factor(revalue(x$V12, c("TRUE" = "T")))
      x$V13 <- as.factor(revalue(x$V13, c("TRUE" = "T")))
      datalist[[i]] <- x # add it to your list
      print(i)
    } else {
      print(i)
       }
}


ids.to.remove <- sapply(datalist, function(i) length(as.data.frame(i)) < 3)
datalist <- datalist[!ids.to.remove]

file <- do.call(rbind, datalist)
file$V1 <- as.character(file$V1)
file <- file[!grepl("PT1092",file$`x$file`),]



# remove found elements

# create data.frame
df <- do.call(rbind, datalist)



colnames(file) <- c("filename", "i", "CHROM", "POS","ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "P8")

grep("1011", file.names)
datalist[[1]] <- "PH"

tryCatch(fread(paste("grep SNW", file.names[1])), error=function(e) NULL) 

file[1,1]
fread(paste("grep SNW", file.names[1])) 
fread(file.names[1])


#ASH2L mutational status





#RBBP5 mutational status







#### make PROVEAN input csvs






x <- read.maf(file.names[714], removeSilent = TRUE, useAll = FALSE)

##CHEXCHEX
unique <- unique(kmt2d.all$folder)
length(unique)
for(i in 1:length(unique)){
  y <- file[grepl("CLL_ICGC",file$`x$file`),]
  z <- kmt2d.all[grepl("CLL_ICGC",kmt2d.all$folder),]
  length(y) == length(z)
}
y <- file[grepl("CLL_Connors",file$`x$file`),]
z <- kmt2d.all[grepl("CLL_Connors",kmt2d.all$folder),]
