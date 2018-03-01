### Count SNVs

setwd("/projects/rmorin/projects/nhl_meta_analysis/results/strelka_pipeline-1.0.0/")
snvlist <- list()

for(i in 1:length(file.names.maf)){
  x <- tryCatch(fread(paste("wc -l", file.names.maf[i]), header=FALSE, sep="\t", stringsAsFactors=FALSE), error=function(e) NULL) 
  x$file <- sapply(strsplit(x$V1, " ", fixed = FALSE), function(x) (x[2]))
  x$V1 <- sapply(strsplit(x$V1, " ", fixed = FALSE), function(x) (as.numeric(x[1])))
  x$V1 <- x$V1-2
  snvlist[[i]] <- x # add it to your list
  print(i)
}

snv <- do.call(rbind, snvlist)
write.table(snv, "/projects/acavalla_prj/kmt2d/301067/patient_data/all_snvs/all_snvs.maf.txt", sep="\t")
