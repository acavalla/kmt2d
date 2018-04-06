## This script relies on maf files concatenated from the server by cat_maf.R. It creates a variety of summary visuals 
## with concatenated maf files that have already been separated into samples in which KMT2D was mutated and those in which it was wildtype.


## Load maftools toolbox and ggplot2 
source("https://bioconductor.org/biocLite.R")
biocLite("maftools")
library(maftools)
library(ggplot2)


## Load concatenated mafs
mutwtcat <- list.files(path = "/projects/acavalla_prj/kmt2d/301067/patient_data/cat_maf/mut", pattern = ".maf$", all.files = FALSE,
                       full.names = TRUE, recursive = FALSE,
                       ignore.case = TRUE, include.dirs = FALSE, no.. = FALSE)
a <- "XXX/xxx/" ##maftools directory path
b <- "XXX/xxx/" ##mafsummary directory path
c <- "XXX/xxx/" ##oncoplot directory path
d <- "XXX/xxx/" ##lolliplot directory path
e <- "XXX/xxx/" ##oncostrip directory path
f <- "XXX/xxx/" ##synthetic lethal oncoplot directory path

dir.create(a, showWarnings = TRUE, recursive = FALSE, mode = "0777")
dir.create(paste0(a, "/", b), showWarnings = TRUE, recursive = FALSE, mode = "0777")
dir.create(paste0(a, "/", c), showWarnings = TRUE, recursive = FALSE, mode = "0777")
dir.create(paste0(a, "/", d), showWarnings = TRUE, recursive = FALSE, mode = "0777")
dir.create(paste0(a, "/", e), showWarnings = TRUE, recursive = FALSE, mode = "0777")
dir.create(paste0(a, "/", f), showWarnings = TRUE, recursive = FALSE, mode = "0777")

x <- list()
name <- list()
type <- list()
for(i in 1:length(mutwtcat)){
  x[[i]] <- read.maf(mutwtcat[i], removeSilent=TRUE, isTCGA = FALSE, verbose = FALSE)
  name <- as.data.frame(strsplit(mutwtcat[i], "/"))
  name <- as.character(name[nrow(name),1])
  type <- ifelse(grepl("wt", name), "WT", "MUT")
  name <- gsub("wt.maf$", "", name)
  name <- gsub("mut.maf$", "", name)
  call <- paste0(name, "_", type)
  basename <- paste0("/projects/acavalla_prj/kmt2d/301067/patient_data/cat_maf/mut/mafsummary/text", "/", call)
  write.mafSummary(maf = x[[i]], basename = basename)
  
  jpeg(paste0(b, "/", call, ".jpg"), height=800, width=1200)
  plot.new()
  mtext(paste(call, "oncoplot"), side=1, line=2, adj=0.5, padj=0.5, cex=1.25, col="black")
  tryCatch(plotmafSummary(maf = x[[i]], rmOutlier = TRUE, addStat = 'median', dashboard = TRUE), error=function(e) NULL)
  dev.off()

  jpeg(paste0(c, "/", call, ".jpg"), height=800, width=1200)
  plot.new()
  mtext(paste(call, "oncoplot"), side=1, line=2, adj=0.5, padj=0.5, cex=1.25, col="black")
  tryCatch(oncoplot(maf = x[[i]], top = 10, fontSize = 12), error=function(e) NULL)
  dev.off()

  jpeg(paste0(d, "/", name, ".jpg"), height = 200, width = 2400)
  tryCatch(lollipopPlot(maf = x[[i]], gene = "KMT2D", AACol = "HGVSp_Short", showMutationRate = TRUE, domainLabelSize = 3, defaultYaxis = FALSE), error=function(e) NULL)
  dev.off()
  
  jpeg(paste0(e, "/", name, ".jpg"), height = 800, width = 1200)
  tryCatch(oncostrip(x[[i]], genes = c("KMT2D", "CDC23", "CHEK1", "MCM7", "PCNA", "POLR2J", "PSMB4", "PSMD2", "RPA3", "STAG2", "ATM")), error=function(e) NULL)
  dev.off()

  jpeg(paste0(f, "/", name, ".jpg"), height = 800, width = 1200)
  tryCatch(oncoplot(x[[i]], genes = c("KMT2D", "CDC23", "CHEK1", "MCM7", "PCNA", "POLR2J", "PSMB4", "PSMD2", "RPA3", "STAG2", "ATM"), fontSize = 12), error=function(e) NULL)
  dev.off()

  
  print(i)
}
