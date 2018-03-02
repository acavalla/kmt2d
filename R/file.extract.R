### This script pulls out the filenames for the varied filenames of the maf, vcf, and vep.vcf files on our server and saves them to
### a txt file to be called in analyses.


source("https://bioconductor.org/biocLite.R")
biocLite("maftools")
library(maftools)


strings <- c("*.passed.maf$", "*assed.somatic.snvs.maf$", "*assed.somatic.indels.maf$", "*_snv_mnv.maf$", "*_indel.maf$", "*assed_indels.maf$",
             "*assed_snvs.maf")
file.names.maf <- list()
for(i in 1:length(strings)){
  x <- list.files(path = "/XXX/xxx/", pattern = strings[i], all.files = FALSE,
                  full.names = TRUE, recursive = TRUE,
                  ignore.case = TRUE, include.dirs = TRUE, no.. = FALSE)
  file.names.maf[[i]] <- x
  print(i)
}

file.names.maf <- unlist(file.names.maf)
file.names.maf <- file.names.maf[!grepl("redundant", file.names.maf)]
file.names.maf <- file.names.maf[!grepl("1092", file.names.maf)]
file.names.maf <- file.names.maf[order(file.names.maf)]
file.names.maf <- droplevels(file.names.maf)

#find out which maf files are unreadable by read.maf (because empty or corrupt) and exclude these. pay attention to condition removeSilent=FALSE
library(futile.logger)
for(i in 1:length(file.names.maf)){
  result <- tryCatch({
    a <- file.names.maf[i]
    a <- paste0("/XXX/xxx/", a)
    x <- read.maf(a, removeSilent=FALSE, isTCGA = FALSE, verbose = FALSE)
  }, error = function(err) {
    # error handler picks up where error was generated
    print(paste("MY_ERROR:  ",err))
    print(i)
    write.table(file.names.maf[i], "/XXX/xxx/", append = TRUE)
  })
}

broken.files <- read.table("/XXX/xxx/", sep = "\t")
broken.files <- broken.files[c(seq(0,nrow(broken.files), 2)),]
broken.files <- gsub("^1 ", "", broken.files)
for(i in 1:length(broken.files)){
  file.names.maf <- file.names.maf[!grepl(broken.files[i], file.names.maf)]
}

write.table(file.names.maf, "/XXX/xxx/", row.names = F, quote = T, col.names=F, sep = "\t")

##load files of gene mutation status
x <- read.table("/XXX/xxx/")
x <- x[!duplicated(x$file),]
x <- x$file
write.table(x, "/XXX/xxx/", row.names = F, quote = F, col.names=F, sep = "\t")




strings <- c("*assed_indels.vcf",
             "*assed_snvs.vcf",
             "*assed.somatic.snvs.vcf",
             "*assed.somatic.indels.vcf",
             "*_indel.vcf",
             "*_snv_mnv.vcf")
file.names.vcf <- list()
for(i in 1:length(strings)){
  x <- list.files(path = "/XXX/xxx/", pattern = strings[i], all.files = FALSE,
                  full.names = TRUE, recursive = TRUE,
                  ignore.case = TRUE, include.dirs = TRUE, no.. = FALSE)
  file.names.vcf[[i]] <- x
  print(i)
}

file.names.vcf <- unlist(file.names.vcf)
file.names.vcf <- file.names.vcf[!grepl("TASK_STRELKA_TCGA", file.names.vcf)] #checked
file.names.vcf <- file.names.vcf[!grepl("MCL_.*TASK_STRELKA.*resul", file.names.vcf)] #checked
file.names.vcf <- file.names.vcf[!grepl("BL_ICGC.*pt_.*pt_", file.names.vcf)] #checked
file.names.vcf <- file.names.vcf[!grepl("BL_L.*PT.*PT", file.names.vcf)] #checked
file.names.vcf <- file.names.vcf[order(file.names.vcf)]
file.names.vcf <- file.names.vcf[!grepl("1092", file.names.vcf)]
file.names.vcf <- file.names.vcf[!grepl("M006", file.names.vcf)]
file.names.vcf <- droplevels(file.names.vcf)

write.table(file.names.vcf, "/XXX/xxx/", row.names = F, quote = F, col.names=F, sep = "\t")

strings <- c("*assed.vep.vcf$", "*assed.somatic.indels.vep.vcf$", "*assed.somatic.snvs.vep.vcf$", "*_indel.vep.vcf$", "*_snv_mnv.vep.vcf$")
file.names.vep.vcf <- list()
for(i in 1:length(strings)){
  x <- list.files(path = "/XXX/xxx/", pattern = strings[i], all.files = FALSE,
                  full.names = TRUE, recursive = TRUE,
                  ignore.case = TRUE, include.dirs = TRUE, no.. = FALSE)
  file.names.vep.vcf[[i]] <- x
  print(i)
}
x <- file.names.vep.vcf
file.names.vep.vcf <- unlist(file.names.vep.vcf)
file.names.vep.vcf <- file.names.vep.vcf[order(file.names.vep.vcf)] 
write.table(file.names.vep.vcf, "/XXX/xxx/", row.names = F, quote = F, col.names=F, sep = "\t")
