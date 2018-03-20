##maftools toolbox: drawing visuals for KMT2D mutational status on MAF files split:
##BL: mut/wt
##CLL: mut/wt
##DLBCL: mut/wt
##FL: mut/wt
##MCL: mut/wt
##MAF files were concatenated using a shell script and 

source("https://bioconductor.org/biocLite.R")
biocLite("maftools")
library(maftools)
library(ggplot2)

mutwtcat <- list.files(path = "/projects/acavalla_prj/kmt2d/301067/patient_data/cat_maf/mut", pattern = ".maf$", all.files = FALSE,
                       full.names = TRUE, recursive = FALSE,
                       ignore.case = TRUE, include.dirs = FALSE, no.. = FALSE)
a <- "XXX/xxx" ##maftools directory path
b <- "XXX/xxx" ##mafsummary directory path
c <- "XXX/xxx" ##oncoplot directory path
d <- "XXX/xxx" ##lolliplot directory path
e <- "XXX/xxx" ##oncostrip directory path
f <- "XXX/xxx" ##synthetic lethal oncoplot directory path

dir.create(paste0(a, "/"), showWarnings = TRUE, recursive = FALSE, mode = "0777")
dir.create(paste0(b, "/"), showWarnings = TRUE, recursive = FALSE, mode = "0777")
dir.create(paste0(c, "/"), showWarnings = TRUE, recursive = FALSE, mode = "0777")
dir.create(paste0(d, "/"), showWarnings = TRUE, recursive = FALSE, mode = "0777")
dir.create(paste0(e, "/"), showWarnings = TRUE, recursive = FALSE, mode = "0777")

x <- list()
name <- list()
type <- list()
for(i in 1:length(mutwtcat)){
  x[[i]] <- read.maf(mutwtcat[i], removeSilent=TRUE, isTCGA = FALSE, verbose = FALSE)
  name[[i]] <- as.data.frame(strsplit(mutwtcat[i], "/"))
  name[[i]] <- as.character(name[[i]][nrow(name[[i]]),1])
  type[[i]] <- ifelse(grepl("wt", name[[i]]), "WT", "MUT")
  name[[i]] <- gsub("wt.maf$", "", name[[i]])
  name[[i]] <- gsub("mut.maf$", "", name[[i]])
  call[[i]] <- paste0(name[[i]], "_", type[[i]])
  basename <- paste0("/projects/acavalla_prj/kmt2d/301067/patient_data/cat_maf/mut/mafsummary/text", "/", call[[i]])
  write.mafSummary(maf = x[[i]], basename = basename)
  
  jpeg(paste0(b, "/", call[[i]], ".jpg"), height=800, width=1200)
  plot.new()
  mtext(paste(call[[i]], "oncoplot"), side=1, line=2, adj=0.5, padj=0.5, cex=1.25, col="black")
  tryCatch(plotmafSummary(maf = x[[i]], rmOutlier = TRUE, addStat = 'median', dashboard = TRUE), error=function(e) NULL)
  dev.off()

  jpeg(paste0(c, "/", call[[i]], ".jpg"), height=800, width=1200)
  plot.new()
  mtext(paste(call[[i]], "oncoplot"), side=1, line=2, adj=0.5, padj=0.5, cex=1.25, col="black")
  tryCatch(oncoplot(maf = x[[i]], top = 10, fontSize = 12), error=function(e) NULL)
  dev.off()

  jpeg(paste0(d, "/", name, ".jpg"), height = 200, width = 2400)
  tryCatch(lollipopPlot(maf = x, gene = "KMT2D", AACol = "HGVSp_Short", showMutationRate = TRUE, domainLabelSize = 3, defaultYaxis = FALSE), error=function(e) NULL)
  dev.off()
  
  jpeg(paste0(e, "/", name, ".jpg"), height = 800, width = 1200)
  tryCatch(oncostrip(x, genes = c("KMT2D", "CDC23", "CHEK1", "MCM7", "PCNA", "POLR2J", "PSMB4", "PSMD2", "RPA3", "STAG2", "ATM")), error=function(e) NULL)
  dev.off()

  jpeg(paste0(f, "/", name, ".jpg"), height = 800, width = 1200)
  tryCatch(oncoplot(x, genes = c("KMT2D", "CDC23", "CHEK1", "MCM7", "PCNA", "POLR2J", "PSMB4", "PSMD2", "RPA3", "STAG2", "ATM"), fontSize = 12), error=function(e) NULL)
  dev.off()

  
  print(i)
}
