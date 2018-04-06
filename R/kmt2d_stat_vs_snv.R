##needs tidying - visualisations of status of KMT2D versus number of SNVs


library(ggplot2)
library(vioplot)
install.packages("ggpubr")
library(ggpubr)
library(dplyr)

##load summary file of sample, kmt2d status and snv/indels
snv.x <- read.table("/projects/acavalla_prj/kmt2d/301067/patient_data/summarybysample.txt")

snv.x <- snv.x[order(snv.x$`file.1`),]
# row.names(snv.x) <- 1:nrow(snv.x)
# kd <- as.data.frame(cbind("file" = snv.x$`file.1`, "type" = as.character(snv.x$type), "total.events" = snv.x$total.events, "mut" = snv.x$mut, "V1" = snv.x$V1))
# kd$total.events <- as.numeric(as.character(kd$total.events))

slim <- snv.kmt2d[,c(1,2,12,20:23)]
atm1 <- make(atm1)
atm1 <- as.data.frame(cbind("match" = kc$match, "mut" = kc$mut))

atm2 <- make(atm2)
atm2 <- as.data.frame(cbind("match" = atm2$match, "mut" = atm2$mut))


atm <- merge(atm1, atm2, by = "match")
ka <- merge(slim, atm, by = "match")
colnames(kcd) <- c("match","kmt2c.mut",  "type", "total.events", "kmt2d.mut", "V1")
kcd$V1 <- factor(kcd$V1, levels = c("BL_WT", "BL_MUT", "CLL_WT", "CLL_MUT", "DLBCL_WT", "DLBCL_MUT", 
                                    "FL_WT", "FL_MUT", "MCL_WT", "MCL_MUT"))
kcd$cdmut <- paste(kcd$kmt2c.mut, kcd$kmt2d.mut, sep="_")
kcd$total.events <- as.numeric(as.character(kcd$total.events))

kdp <- merge(p1, kd, by = "V1")
colnames(kdp) <- c("match", "p1.mut", "total.events", "kmt2d.mut", "V1")
kdp$V1 <- factor(kdp$V1, levels = c("BL_WT", "BL_MUT", "CLL_WT", "CLL_MUT", "DLBCL_WT", "DLBCL_MUT", 
                                    "FL_WT", "FL_MUT", "MCL_WT", "MCL_MUT"))
kdp$kpmut <- paste(kdp$p1.mut, kdp$kmt2d.mut, sep="_")



kcd.tumor <- merge(kc.tumor, kd.tumor, by = "V1")
colnames(kcd.tumor) <- c("match", "kmt2c.mut", "total.events", "kmt2d.mut", "V1")
kcd.tumor$V1 <- factor(kcd.tumor$V1, levels = c("BL_WT", "CLL_WT", "CLL_MUT", "DLBCL_WT", "DLBCL_MUT", 
                                                "FL_WT", "FL_MUT", "MCL_WT", "MCL_MUT"))
kcd.tumor$cdmut <- paste(kcd.tumor$kmt2c.mut, kcd.tumor$kmt2d.mut, sep="_")

kcd.normal <- merge(kc.normal, kd.normal, by = "V1")
colnames(kcd.normal) <- c("match", "kmt2c.mut", "total.events", "kmt2d.mut", "V1")
kcd.normal$V1 <- factor(kcd.normal$V1, levels = c("BL_WT", "CLL_WT", "CLL_MUT", "DLBCL_WT", "DLBCL_MUT", 
                                                  "FL_WT", "FL_MUT", "MCL_WT", "MCL_MUT"))
kcd.normal$cdmut <- paste(kcd.normal$kmt2c.mut, kcd.normal$kmt2d.mut, sep="_")


sic.pval <- k[k$'type.1' != "BL",]
sic.pval$'type.1' <- factor(sic.pval$'type.1', levels = c("CLL", "DLBCL", "FL", "MCL"))
sic.pval$mut <- factor(sic.pval$mut, level = c("MUT", "WT"))
dat <- data.frame(
  x=1.5, y=55,
  pval = sapply(split(sic.pval, sic.pval$mut), function(x) t.test(total.events ~ x$'type.1')$p.value),
  variable=factor(sic.pval, levels= c("CLL", "DLBCL", "FL", "MCL")
  ))


pval <- sapply(split(sic.pval, sic.pval$mut), function(x) t.test(total.events ~ mut, x)$p.value)
variable <- c("CLL", "DLBCL", "FL", "MCL")

snv.kmt$mut <- ifelse(is.na(snv.kmt$mut.type), "WT", "MUT")
snv.kmt$mut <- factor(snv.kmt$mut, levels = c("WT", "MUT")) 
snv.kmt$V1 <- paste(snv.kmt$type, snv.kmt$mut, sep="_")

snv.kmt.snv <- snv.kmt[grep("snv", snv.kmt$filename),]
snv.kmt.ind <- snv.kmt[grep("ind", snv.kmt$filename),]
snv.kmt.snv$non <- ifelse(grepl("nonsense", snv.kmt.snv$mut.type), "Nonsense",
                          ifelse(is.na(snv.kmt.snv$mut.type), "WT",
                                 "MUT"))
snv.kmt.snv$factor <- factor(snv.kmt.snv$V1, levels=c("BL_WT", "CLL_WT", "CLL_MUT", "DLBCL_WT", "DLBCL_MUT", 
                                                      "FL_WT", "FL_MUT", "MCL_WT", "MCL_MUT"))

snv.kmt.all <- merge(kmt2d.all, snv, by = "filename", all.x = TRUE)
snv.kmt.all$pos <- as.numeric(snv.kmt.all$pos)
snv.kmt.all.snv <- snv.kmt.all[grep("snv", snv.kmt.all$filename),]

snv.kmt.all$type <- sapply(strsplit(snv.kmt.all$folder, "_", fixed = FALSE), function(x) (x[1]))
snv.kmt.all$type <- as.factor(snv.kmt.all$type)

snv.x$type <- factor(snv.x$type, levels=c("DLBCL", "FL", "CLL", "BL", "MCL"), ordered=TRUE)
counts <- table(snv.x$mut, snv.x$type)
jpeg("/projects/acavalla_prj/kmt2d/301067/patient_data/all_snvs/stack_bar.jpg", height = 3500, width = 2500)
barplot(counts, main="Proportion of mutated to wildtype KMT2D in samples by lymphoma type",
        xlab="NHL type", ylab = "Samples", col=c("darkblue","red"),
        legend = rownames(counts), cex.names = 5, cex.axis = 5)
dev.off()

jpeg("/projects/acavalla_prj/kmt2d/301067/patient_data/all_snvs/stack_bar.jpg", height = 3500, width = 2500)
ggplot(snv.x, aes(type)) + geom_bar(aes(fill = mut), position="stack") + theme_minimal() +
  theme(text = element_text(size = 65)) + 
  ggtitle("Proportion of mutated KMT2D in samples by lymphoma type")
dev.off()


ggplot(snv.x, aes(x = snv.x$type, y = snv.x$mut)) + geom_col(position="stack")



tmp <- rbind(snv.x[snv.x$type == "BL",], snv.x[snv.x$type == "CLL",], snv.x[snv.x$type == "DLBCL",], 
             snv.x[snv.x$type == "FL" & snv.x$total.events < 30000,], snv.x[snv.x$type == "MCL" & snv.x$total.events < 1000,])

events.total.1 <- boxplot(as.numeric(as.character(total.events))~mut, data=snv.x, plot=0)
ggplot(tmp, aes(fill = mut, x = total.events)) +
  geom_density(alpha=0.3) + theme_minimal() + 
  theme(text = element_text(size = 20)) + 
  ggtitle("Density distribution of somatic alterations versus KMT2D mutational status") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x = "Total somatic alterations") + facet_wrap(~type, scales = "free")

jpeg("/projects/acavalla_prj/kmt2d/301067/patient_data/all_snvs/boxplot_facet.jpg", height=1200, width=2000)
ggplot(tmp, aes(x=mut, y=total.events, fill=mut)) + 
  geom_boxplot(alpha =0.3) + 
  theme_minimal() +
  theme(axis.title=element_text(size=32), strip.text = element_text(size=60), axis.text=element_text(size=20), legend.position = "none") +
  ylab("Total somatic events") +
  xlab("KMT2D mutational status") +
  facet_wrap(~type, scales = "free")
dev.off()

CLLind <- wilcox.test(V1.2 ~ mut, data=snv.x[snv.x$type=="CLL",], paired = FALSE) 
DLBCLind <- wilcox.test(V1.2 ~ mut, data=snv.x[snv.x$type=="DLBCL",], paired = FALSE) 
FLind <- wilcox.test(V1.2 ~ mut, data=snv.x[snv.x$type=="FL",], paired = FALSE) 
MCLind <- wilcox.test(V1.2 ~ mut, data=snv.x[snv.x$type=="MCL",], paired = FALSE)


ind <- rbind(snv.x[snv.x$type == "BL" & snv.x$V1.2 < 15000,], snv.x[snv.x$type == "CLL" & snv.x$V1.2 < 1000,], snv.x[snv.x$type == "DLBCL" & snv.x$V1.2 < 4000,], 
            snv.x[snv.x$type == "FL" & snv.x$V1.2 < 1500,], snv.x[snv.x$type == "MCL" & snv.x$V1.2 < 20,])

jpeg("/projects/acavalla_prj/kmt2d/301067/patient_data/all_snvs/boxplot_indel_facet.jpg", height=1500, width=1500)
ggplot(ind, aes(x=mut, y=V1.2, fill=mut)) + 
  geom_boxplot(alpha =0.3) + 
  theme_minimal() +
  theme(axis.title=element_text(size=32), strip.text = element_text(size=60), axis.text=element_text(size=20), legend.position = "none") +
  ylab("Total indels") +
  xlab("KMT2D mutational status") +
  facet_wrap(~type, scales = "free")
dev.off()

CLLsnv <- wilcox.test(V1.1 ~ mut, data=snv.x[snv.x$type=="CLL",], paired = FALSE) 
DLBCLsnv <- wilcox.test(V1.1 ~ mut, data=snv.x[snv.x$type=="DLBCL",], paired = FALSE) 
FLsnv <- wilcox.test(V1.1 ~ mut, data=snv.x[snv.x$type=="FL",], paired = FALSE) 
MCLsnv <- wilcox.test(V1.1 ~ mut, data=snv.x[snv.x$type=="MCL",], paired = FALSE)
#fake <- data.frame(file = "fake", events = 99999999, gene = "fake", assembly = "fake", chr = "fake", startpos = "fake", endpos = "fake", ref = "fake", alt = "fake", alt2 = "fake", 
#                                             mut.type = "fake", i = "fake", type = "BL", mut = "MUT", V1 = "fake")
snv <- rbind(snv.x[snv.x$type == "BL" & snv.x$V1.1 < 15000,], snv.x[snv.x$type == "CLL" & snv.x$V1.1 < 8700,], snv.x[snv.x$type == "DLBCL" & snv.x$V1.1 < 40000,], 
             snv.x[snv.x$type == "FL" & snv.x$V1.1 < 25000,], snv.x[snv.x$type == "MCL" & snv.x$V1.1 < 400,])
jpeg("/projects/acavalla_prj/kmt2d/301067/patient_data/all_snvs/boxplot_snv_facet.jpg", height=1500, width=1500)
ggplot(snv, aes(x=mut, y=V1.1, fill=mut)) +
  geom_boxplot(alpha =0.3) + 
  theme_minimal() +
  theme(axis.title=element_text(size=32), strip.text = element_text(size=60), axis.text=element_text(size=20), legend.position = "none") +
  ylab("Total SNVs") +
  xlab("KMT2D mutational status") +
  facet_wrap(~type, scales = "free")
dev.off()

ind <- wilcox.test(V1.2 ~ mut, data=snv.x, paired = FALSE) 
jpeg("/projects/acavalla_prj/kmt2d/301067/patient_data/all_snvs/boxplot_indel.jpg", height=1500, width=1000)
ggplot(snv.x, aes(x=mut, y=V1.2, fill=mut)) +
  geom_boxplot(alpha =0.3) + 
  theme_minimal() +
  theme(axis.title=element_text(size=32), strip.text = element_text(size=60), axis.text=element_text(size=20), legend.position = "none") +
  coord_cartesian(ylim = c(0, 2000)) + 
  ylab("Total indels") +
  xlab("KMT2D mutational status")
dev.off()


jpeg("/projects/acavalla_prj/kmt2d/301067/patient_data/all_snvs/boxplot_sorted.jpg", height=500, width=500)
snv.x$type <- with(snv.x, reorder(type, -total.events, median))
ggplot(snv.x, aes(x = type, y = total.events)) +
  geom_boxplot(fill = c("yellow", "green", "cyan", "pink", "red"), alpha = 0.3) + theme_minimal() +
  scale_x_discrete(name = "NHL type") +
  scale_y_continuous(name = "Total somatic mutations") +
  coord_cartesian(ylim = c(0, 45000)) +
  ggtitle("Boxplot of total somatic mutation events by lymphoma type")
dev.off()

ggplot(snv.xlog, aes(x = reorder(type, total.events, median), y = (total.events/3234.83))) +
  geom_boxplot(fill = c("red", "pink", "cyan", "green", "yellow"), alpha = 0.3) +
  coord_flip(ylim = c(0,45)) +
  theme_minimal() +
  scale_x_discrete(name = "NHL type") +
  scale_y_continuous(name = "Total somatic mutations") +
  ggtitle("Boxplot of total somatic mutation events by lymphoma type")
dev.off()
3234.83
MBp





library(beeswarm)
beeswarm(total.events~mutdel, data = snv.kmt2d,
         pch = 16,
         xlab = "", ylab = "number of SNVs+indels",
         labels = c("MUT_Deleterious", "MUT_Neutral", "WT"),
         ylim = c(0,70000),
         corral = "gutter",
         col = c("red", "green", "green"))
bxplot(total.events~ mutdel, data = snv.kmt2d, add = TRUE)
title("SNVs and indels in genome by KMT2D status")







events.total.b <- boxplot(as.numeric(as.character(total.events))~type, data=kd, plot=0)
events.total <- boxplot(as.numeric(as.character(kd$total.events))~type,data=kd, las = 2,
                        col=(c("purple", "yellow", "red", "green", "blue")),
                        ylim=c(0,75000),
                        main = "SNVs and indels in genome by cancer type",
                        names=paste(events.total.b$names, "(n=", events.total.b$n, ")"),
                        notch = TRUE,
                        par(mar = c(12, 5, 4, 2)+ 0.1) )
abline(h = median(as.numeric(as.character(kd$total.events))))

library(ggplot2)
countwt <- nrow(kd[kd$mut=="WT",])
countmut <- nrow(kd[kd$mut=="MUT",])


geom_violin(mapping = x1, x2, names=c(paste0("WT (n=", countwt, ")"), paste0("MUT (n=", countmut, ")")),
        col=c("red", "green"))
abline(h = median(as.numeric(as.character(kd$total.events))))
title("SNVs and indels in genome")
p <- ggplot(kd, aes(factor(kd$mut), kd$total.events))
p + geom_violin(scale = "count",
                aes(fill = factor(mut)))
                + geom_jitter(height = 0, width = 0.1)



events.factor.b <- boxplot(total.events~mut, data=kd, plot=0)
events.factor <- boxplot(total.events~mut,data=snv.ind.comb, las = 2,
                         col=(c("red","blue")),
                         names=paste(events.factor.b$names, "(n=", events.factor.b$n, ")"),
                         main = "KMT2D status versus number of SNV and indels in genome",
                         notch = TRUE)
events.factor.ylim <- boxplot(total.events~mut,data=snv.ind.comb, las = 2,
                              col=(c("red","blue")),
                              ylim=c(0,35000),
                              main = "KMT2D status versus number of SNV and indels in genome",
                              names=paste(events.factor.b$names, "(n=", events.factor.b$n, ")"),
                              notch = TRUE)

snvs.type <- boxplot(events~type,data=snv.kmt.snv, las = 2,
                     col=(c("red","blue")),
                     ylim=c(0,50000),
                     main = "SNVs versus cancer type",
                     notch=TRUE)
inds.type <- boxplot(events~type,data=snv.kmt.ind, las = 2,
                     col=c("red","blue"),
                     ylim=c(0,8000),
                     main = "Indels versus cancer type",
                     notch=TRUE)

snvs.inds.type.b <- boxplot(total.events~mut, data=k, plot=0)
snvs.inds.type <- boxplot(total.events~mut,data=k, las = 2,
                          col=c("red","blue"),
                          ylim=c(0,30000),
                          main = "SNV + indel events by cancer type and KMT2D status",
                          names=paste(snvs.inds.type.b$names, "(n=", snvs.inds.type.b$n, ")"),
                          notch=TRUE)


snvs.inds.type.split.b <- boxplot(total.events~V1, data=k, plot=0)
snvs.inds.type.split <- boxplot(total.events~V1,data=k, las = 2,
                                col=c("red","blue"),
                                ylim=c(0,30000),
                                main = "SNV + indel events by cancer type and KMT2D status",
                                par(mar = c(12, 5, 4, 2)+ 0.1),
                                names=paste(snvs.inds.type.split.b$names, "(n=", snvs.inds.type.split.b$n, ")"),
                                at =c(1,3,4,6,7,9,10,12,13))
mut.wt.b <- boxplot(events~mut, data=snv.kmt.snv, plot=0)
mut.wt <- boxplot(events~mut, data=snv.kmt.snv, las = 2,
                  col=c("red", "blue"),
                  names=paste(mut.wt.b$names, "(n=", mut.wt.b$n, ")"),
                  main = "SNVs versus KMT2D status")
mut.wt.lim <- boxplot(events~mut, data=snv.kmt.snv, las =2,
                      col=c("red", "blue"),
                      ylim=c(0,29000),
                      names=paste(mut.wt.b$names, "(n=", mut.wt.b$n, ")"),
                      main = "SNVs versus KMT2D status",
                      notch = TRUE)
df2<-melt(k,id.var=c("type.1","Name"))
ggplot(k, aes(x=V1,y=total.events,fill=factor(V1))+
         geom_boxplot() + labs(title="SNV + indel events by cancer type and KMT2D status") +facet_wrap(~variable))

library(sm)
attach(mtcars)
snv.kmt2d$total.events <- as.integer(snv.kmt2d$total.events)
plot((density(snv.kmt2d$total.events)) +
       density(snv.kmt2d$total.events[snv.kmt2d$mut == "MUT"]) +
       density(snv.kmt2d$total.events[snv.kmt2d$mut == "WT"]))

# create value labels 
x <- as.data.frame(cbind(snv.kmt2d$mutdel, as.factor(snv.kmt2d$mutdel)))
cyl.f <- factor(x, levels= c("MUT_Del", "MUT_Neutral", "WT_Neutral"),
                labels = c("Deleterious mutation", "Neutral mutation", "Wild-type")) 

# plot densities 
sm.density.compare(snv.kmt2d$total.events, snv.kmt2d$V1, xlab="total genomic somatic alterations")
title(main="somatic alterations by kmt2d mutational status")

# add legend via mouse click
colfill<-c(2:(2+length(levels(cyl.f)))) 
legend(locator(1), levels(cyl.f), fill=colfill)


par(mfrow=c(3, 3))
colnames <- dimnames(snv.kmt2d)[[2]]
for (i in 2:8) {
  hist(crime[,i], xlim=c(0, 3500), breaks=seq(0, 3500, 100), main=colnames[i], probability=TRUE, col="gray", border="white")
  d <- density(crime[,i])
  lines(d, col="red")
}par(mfrow=c(3, 3))

colnames <- dimnames(crime.new)[[2]]
for (i in 2:8) {
  d <- density(crime[,i])
  plot(d, type="n", main=colnames[i])
  polygon(d, col="red", border="gray")
}




nomcl <- snv.kmt2d[snv.kmt2d$type.1 != "MCL",]

x1 <- nomcl$total.events[nomcl$mutdel == "WT_Neutral"]
x2 <- nomcl$total.events[nomcl$mutdel == "MUT_Neutral"]
x3 <- nomcl$total.events[nomcl$mutdel == "MUT_Deleterious"]
vioplot(x3, x2, x1, names=c("MUT_Deleterious (n=92)", "MUT_Neutral (n=40)", "WT (n=488)"), 
        col=c("red", "green", "green"))
boxplot(total.events~mutdel,data=snv.kmt2d[snv.kmt2d$type.1 != "MCL",], ylim=c(0,60000), las = 2, notch = TRUE)
boxplot(total.events~mut,data=snv.kmt2d[snv.kmt2d$type.1 != "MCL",], ylim=c(0,60000), las = 2, notch = TRUE)



plot(density(nomcl$total.events[nomcl$mutdel=="MUT_Deleterious"]))

d <- density(kd$mpg)
plot(d, main="Kernel Density of Miles Per Gallon")
polygon(d, col="red", border="blue")


# create value labels 
attach(snv.kmt2d)
cyl.f <- factor(mutdel, levels= c("MUT_Del", "MUT_Neutral", "WT_Neutral"),
                labels = c("Deleterious mutation", "Neutral mutation", "Wild-type")) 
x <- as.factor(mutdel)
# plot densities 
sm.density.compare(snv.kmt2d$total.events, x)

title(main="distributions of total somatic alterations by kmt2d status")

colfill<-c(2:(2+length(levels(cyl.f)))) 
legend(locator(1), levels(cyl.f), fill=colfill)


b <- boxplot(total.events~mut, data=k, plot=0)
boxplot(total.events ~ mut, data=k,
        col=(c("red","blue")),
        names=paste(b$names, "(n=", b$n, ")"), las =2, notch = TRUE,
        ylim=c(0,30000),
        main = "KMT2D status versus number of SNVs in genome \n pval=0.002655671")

ggplot(snv.kmt.all[snv.kmt.all$pos > 49200000,], aes(pos, events)) +
  geom_point(aes(colour = type)) +
  #geom_smooth(method ="lm") +
  coord_cartesian() +
  theme_bw()


types <- c("BL", "CLL", "DLBCL", "FL", "MCL")
snv.x$mut <- as.factor(snv.x$mut)
for(i in 1:length(types)){
  mytable <- table(snv.x$mut[snv.x$type==types[i]])
  lbls <- paste(names(mytable), "\n", mytable, sep="")
  jpeg(paste0("/projects/acavalla_prj/kmt2d/301067/patient_data/gene_muts/kmt2dpie_", types[i], ".jpg"), height=500, width=750)
  pie(mytable, labels = lbls, 
      main=paste("KMT2D mutation status:", types[i], "samples"))
  dev.off()
}


library(plotrix)
mytable <- table(kcd$cdmut)
lbls <- paste(names(mytable), "\n", mytable, sep="")
pie(mytable, labels = lbls, 
    main="Distribution of KMT2C and KMT2D mutation status, respectively")

mytable <- table(kcd$cdmut)
lbls <- paste(names(mytable), "\n", mytable, sep="")
pie(mytable, labels = lbls, 
    main="Distribution of KMT2C and KMT2D mutation status within tumours, respectively")

mytable <- table(kcd.tumor$cdmut)
lbls <- paste(names(mytable), "\n", mytable, sep="")
pie(mytable, labels = lbls, 
    main="Distribution of KMT2C and KMT2D mutation status within tumours, respectively")

mytable <- table(kcd.normal$cdmut)
lbls <- paste(names(mytable), "\n", mytable, sep="")
pie(mytable, labels = lbls, 
    main="Distribution of KMT2C and KMT2D mutation status within normal tissue, respectively")

mytable <- table(kdp$kpmut)
lbls <- paste(names(mytable), "\n", mytable, sep="")
pie(mytable, labels = lbls, 
    main="Distribution of KMT2D and PARP1 mutation status within tumours, respectively")

g <- ggplot(cdtypes, aes(V1))
g + scale_y_continuous(labels = scales::percent) +
  labs ( y = "distribution of kmt2c and kmt2d mut types") +
  geom_bar(aes(fill = V2), color = "black", position = "fill")
#geom_text(aes(y=cumsum, label=label), vjust=-1)

g <- ggplot(k, aes(type.1))
g + scale_y_continuous(labels = scales::percent) +
  labs ( y = "percentage of cases with at least one KMT2D mutation") +
  geom_bar(aes(fill = mut), color = "black", position = "fill")
#geom_text(aes(y=cumsum, label=label), vjust=-1)

g <- ggplot(order, aes(num,log10(total.events))) + 
  labs ( title = "individual cases and number of genomic mutation events") +
  geom_point(aes(colour = factor(kmt2d.mut))) + 
  facet_wrap(~type)
g
kcd$match <- as.character(kcd$match)
order <- setorder(kcd, -total.events)
row.names(order) <- NULL
order$num <- rownames(order)
order$num <- as.numeric(order$num)




p <- ggplot(tmp, aes(fill = mut, x = total.events)) +
  geom_density(alpha=0.3) + theme_minimal() + 
  theme(text = element_text(size = 12)) + 
  ggtitle("Density distribution of somatic alterations versus KMT2D mutational status") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x = "Total somatic alterations") + facet_wrap(~type, scales = "free")
p + stat_function(stat_compare_means())

CLL <- wilcox.test(total.events ~ mut, data=kd[kd$type=="CLL",], paired = FALSE) 
DLBCL <- wilcox.test(total.events ~ mut, data=kd[kd$type=="DLBCL",], paired = FALSE) 
FL <- wilcox.test(total.events ~ mut, data=kd[kd$type=="FL",], paired = FALSE) 
MCL <- wilcox.test(total.events ~ mut, data=kd[kd$type=="MCL",], paired = FALSE)


q.p <- wilcox.test(total.events ~ mut, data=kd)
jpeg('/projects/acavalla_prj/kmt2d/301067/patient_data/density_dist.jpg', height = 900, width = 3000)
ggplot(kd, aes(fill = mut, x = total.events)) +
  scale_x_continuous(limits = c(0, 30000)) +
  geom_density(alpha=0.3) + theme_minimal() + 
  theme(text = element_text(size = 20)) + 
  ggtitle(paste("Density distribution of somatic alterations versus KMT2D mutational status\np = ", q.p[[3]])) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x = "Total somatic alterations")
dev.off()



g <- ggplot(kcd, aes(type,kmt2d.mut)) + 
  labs (title = "sample sizes and kmt2d status") +
  geom_col(colour = aes(as.character(kmt2d.mut)), position = "stack")
g
counts <- table(kcd$kmt2d.mut, kcd$type)
barplot(counts, main="sample sizes and kmt2d status",
        xlab="type", col=c("darkblue","red"),
        legend = rownames(counts))

####writing plots to file
jpeg('/projects/acavalla_prj/kmt2d/301067/patient_data/density_dist_facet.jpg', height = 900, width = 3000)
ggplot(tmp, aes(fill = mut, x = total.events)) +
  geom_density(alpha=0.3) + theme_minimal() + 
  theme(text = element_text(size = 20)) + 
  ggtitle("Density distribution of somatic alterations versus KMT2D mutational status") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x = "Total somatic alterations") + facet_wrap(~type, scales = "free")
dev.off()

library(circlize)
circos.par("track.height" = 0.1)
circos.initialize(x = stacked_bar)
