### The result printed to the console should be copied to the shell, preceded by "cat " and followed by 
## "| awk '/Hugo/&&c++>0 {next} 1' | awk '/version/&&c++>0 {next} 1' > outfile". It segregates files into whether they come from a sample
## with mutated KMT2D or wildtype (or any other gene, depending on the gene_mut.R script) and concatenates those mafs to 
## be entered into the maf_tools.R script


x <- read.table("/projects/acavalla_prj/kmt2d/301067/patient_data/summarybysample.txt", sep = "\t", fill = TRUE, stringsAsFactors = FALSE, header = TRUE)

a <- c("BL_ICGC_MALY-DE_orig", "BL_Love_Dave_2012", "CLL", "batch7",
       "DLBCL_ICGC", "DLBCL_TCGA", "FL_Connors",
       "FL_Shah", "MCL_Bea_Campo_2013", "MCL_Sweden", "MCL_Zhang_Dave_2014")

ff <- list()
for(i in 1:length(a)){
    ff[[i]] <- x[grepl(a[[i]], x[,4]),]
  }

wt <- list()
mut <- list()
for(i in 1:length(ff)){
  mut[[i]] <- ff[[i]][ff[[i]]$mut=="MUT",]
  wt[[i]] <- ff[[i]][ff[[i]]$mut=="WT",]
}

z <- list(wt,mut)

for(h in 1:length(z)){
  y <- z[[h]]
  for(i in 1:length(y)){
    y[[i]] <- paste(y[[i]]$mafsnv, y[[i]]$mafindel)
  }
  if(h==1){
    wt <- y
  } else {
    mut <- y
  }
}

z <- list(wt,mut)

for(h in 1:length(z)){
  y <- z[[h]]
  for(i in 1:length(y)){
    for (j in 1:length(y[[i]])) {cat(y[[i]][j], "")}
    cat("\n")
  }
}



