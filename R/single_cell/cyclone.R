##Analysis of where cells are in cell cycle based on sc transcriptomics

library(repmis)
library(scran)
library(preprocessCore)
library(tools)
source_data("https://github.com/PMBio/cyclone/blob/c982e7388d8e49e1459055504313f87bb3eb0ceb/R/pairs_method/core/pairs_functions.RData?raw=true")

gene.ids <- read.table("/home/elim/single_cell_projects/TRANSCRIPTOME/hodgkin_lymphoma_cell_line/SI-3A-B2/outs/filtered_gene_bc_matrices/hg19/genes.tsv")
mouse.ids <- read.table("/projects/acavalla_prj/summaries/normalised/one2onemouse.txt", sep = "\t")

mouse.ids <- mouse.ids[-1,]
row.names(mouse.ids) <- NULL
onefactor <- mouse.ids$V2 == "ortholog_one2one"
oneone <- mouse.ids[onefactor,]

names(oneone)[1] <- "human"
names(oneone)[3] <- "mouse"
oneone <- cbind(oneone[1], oneone[3])

G1.marker.pairs <- G1.marker.pairs[(G1.marker.pairs$Gene.1 %in% oneone$mouse) & (G1.marker.pairs$Gene.2 %in% oneone$mouse),]
S.marker.pairs <- S.marker.pairs[(S.marker.pairs$Gene.1 %in% oneone$mouse) & (S.marker.pairs$Gene.2 %in% oneone$mouse),]
G2M.marker.pairs <- G2M.marker.pairs[(G2M.marker.pairs$Gene.1 %in% oneone$mouse) & (G2M.marker.pairs$Gene.2 %in% oneone$mouse),]
training.data <- training.data[row.names(training.data) %in% oneone$mouse,]

pairs <- list(G1.marker.pairs, S.marker.pairs, G2M.marker.pairs)

sandbag <- sandbag(training.data, is.G1 = id.G1, is.S = id.S, is.G2M = id.G2M, gene.names=rownames(training.data), 
                   fraction=0.5, subset.row=NULL)


#link normalised gene matrix to mouse IDs, run machine learning algorithm to predict cell cycle
path = "/projects/acavalla_prj/summaries/normalised"
out.file<-""
file.names <- dir(path, pattern =".normalised.txt")

for(i in 2:length(file.names)){
  setwd("/projects/acavalla_prj/summaries/normalised")
  data <- read.table(file.names[i], sep = "\t")
  #sandbag1 <- read.table("/projects/acavalla_prj/summaries/getcell.cycle.stages/sandbag.G1.txt")
  #sandbag2 <- read.table("/projects/acavalla_prj/summaries/cell.cycle.stages/sandbag.S.txt")
  #sandbag3 <- read.table("/projects/acavalla_prj/summaries/cell.cycle.stages/sandbag.G2M.txt")
  #sandbag <- list(sandbag1, sandbag2, sandbag3)
  #data <- read.table("/projects/acavalla_prj/summaries/normalised/SI-3A-A1.normalised.txt")
  
  data <- as.matrix(data)
  normalize.quantiles(data,copy=FALSE)
  row.names(data) <- gene.ids$V1
  
  #Use authors' approach: only 1:1 orthologues
  
  data <- data[row.names(data) %in% oneone$human,]
  data <- merge(oneone, data, by.x = "human", by.y = "row.names")
  row.names(data) <- data$mouse
  data <- data[,-1]
  data <- data[,-1]
  
  data <- as.matrix(data)
  cell.cycle.stages <- cyclone(data, sandbag, gene.names=rownames(data), iter=1000, min.iter=100, min.pairs=50,
                               BPPARAM=SerialParam(), verbose=FALSE, subset.row=NULL)
  
  setwd("/projects/acavalla_prj/summaries/cell.cycle.stages/")
  
  filename.ccs <- paste(file_path_sans_ext(file.names, compression = FALSE)[i], ".ccstages.txt", sep="")
  write.table(cell.cycle.stages, filename.ccs, sep="\t",quote=FALSE)
}

write.table(sandbag$G1, "/projects/acavalla_prj/summaries/cell.cycle.stages/sandbag.G1.txt", sep="\t",quote=FALSE)
write.table(sandbag$S, "/projects/acavalla_prj/summaries/cell.cycle.stages/sandbag.S.txt", sep="\t",quote=FALSE)
write.table(sandbag$G2M, "/projects/acavalla_prj/summaries/cell.cycle.stages/sandbag.G2M.txt", sep="\t",quote=FALSE)

path = "/projects/acavalla_prj/summaries/cell.cycle.stages"
out.file<-""
file.names <- dir(path, pattern =".ccstages.txt")

#scatterplots
for(i in 1:length(file.names)){
  cell.cycle.stages <- read.table(file.names[i], sep = "\t")
  col <- character(length(cell.cycle.stages$phases))
  col[cell.cycle.stages$phases=="G1"] <- "red"
  col[cell.cycle.stages$phases=="G2M"] <- "blue"
  col[cell.cycle.stages$phases=="S"] <- "darkgreen"
  plot(cell.cycle.stages$scores.G1, cell.cycle.stages$scores.G2M, col=col, pch=16,
       main="Cell cycle stages", sub=file.names[i],
       xlab="G1 score", ylab="G2M score")
}

#percentage of cells per stage by run
library(data.table)

path = "/projects/acavalla_prj/summaries/cell.cycle.stages"
out.file<-""
file.names <- dir(path, pattern ="ccstages.txt")

for(i in 1:length(file.names)){
  setwd("/projects/acavalla_prj/summaries/cell.cycle.stages")
  file <- read.table(file.names[i], sep = "\t")
  
  G1 <- file[file$phases=="G1",]
  G1 <- nrow(G1)
  G2M <- file[file$phases=="G2M",]
  G2M <- nrow(G2M)
  S <- file[file$phases=="S",]
  S <- nrow(S)
  
  all.phases <- cbind(G1, S, G2M)
  names(all.phases)[1] <- i
  
  out.file <- rbind(out.file, all.phases)
}

out.file <- out.file[-1,]
row.names(out.file) <- NULL

out.file <- t(out.file)
out.file <- as.data.table(out.file)

out.file$G1 <- as.numeric(out.file$G1)
out.file$G2M <- as.numeric(out.file$G2M)
out.file$S <- as.numeric(out.file$S)

rownames(out.file) <-  c("HEK293", "OCI-LY10 1", "Hodgkin lymphoma 1", "OCI-LY10 2", "Hodgkin lymphoma 2", "OCI-LY10 3",
           "Hodgkin lymphoma 3", "WSU-NHL 1", "Hodgkin lymphoma 4", "WSU-NHL 2", "Hodgkin lymphoma 5", "WSU-NHL 3", 
           "Hodgkin lymphoma 6")

library(ggplot2)
library(scales)

names <- rownames(out.file)
rownames(out.file) <- NULL
out.file <- cbind(names,out.file)

dfm = melt(out.file, id.vars = "names")
ggplot(dfm,aes(names, y = value), las=2) + 
  geom_bar(aes(fill=variable),position = "fill",stat = "identity") + 
  scale_y_continuous(labels = percent_format())

