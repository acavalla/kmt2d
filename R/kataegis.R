##Kataegis analysis
source("https://bioconductor.org/biocLite.R")
biocLite("maftools")
library(maftools)

file.names <- read.table("/projects/acavalla_prj/kmt2d/301067/patient_data/files/filenames.maf.txt", sep = "\n", fill = TRUE, stringsAsFactors = FALSE)
file.names <- file.names[!grepl("indel|INDEL", file.names$V1),]

fl <- list()

a <- list("BL", "CLL", "DLBCL", "FL", "MCL")
for(i in 1:length(a)){
  fl[[i]] <- file.names[grepl(a[[i]], file.names)]
}

for(i in 1:length(fl)){
  y <- fl[[i]]
  for(j in 1:length(y)){
    x <- tryCatch(read.maf(y[j], removeSilent = FALSE), error=function(e) NULL)
    if(length(x) != 0){
      name <- as.data.frame(strsplit(y[j], "/")) 
      type <- as.data.frame(strsplit(as.character(name[2,1]), "_"))
      type <- type[1,1]
      
      dir.create(paste0('/projects/acavalla_prj/kmt2d/301067/patient_data/kataegis', "/", type), showWarnings = TRUE, recursive = FALSE, mode = "0777")
      name <- ifelse(grepl("FL|CLL_Connors", y[j]), as.character(name[3, 1]),
                                  as.character(name[nrow(name), 1]))
      jpeg(paste0('/projects/acavalla_prj/kmt2d/301067/patient_data/kataegis', "/", type, "/", name, ".jpg"), height = 400, width = 2000)
      b <- ifelse(grepl("MCL|BL_Love", name), "hg38", "hg19")
      rainfallPlot(x,
                   tsb = NULL, detectChangePoints = FALSE,
                   ref.build = b, color = NULL, savePlot = FALSE, width = 6,
                   height = 3, fontSize = 22, pointSize = 1)
      dev.off()
    }
    print(j)
    print(y[j])
  }
  print(i)
}
