
setwd("~/NGS/HCC.NGS.CNV/rscnvs")
cn.panel <- read_delim("new.hcc.rscnvs.txt", delim = "\t")
cnv.panel <- as.data.frame(cn.panel)
cn.panel <- join(cytobands, cnv.panel, type = "right", by = "Cytoband")
genes = cn.panel[,c(2:4,1)]
#setwd("~/NGS/HCC.NGS.CNV/All.samples.250kbp.fix.iod.masked")
#setwd("~/NGS/HCC.NGS.CNV/clustering/all.var.250kbp.normreadct.masked")
setwd("~/NGS/HCC.NGS.CNV/clustering/all.1MB.fix.normreadct")
segnorm <- read_delim("SegNorm.txt", delim = "\t")
segfix <- read_delim("SegFixed.txt", delim = "\t")
segcopy <- read_delim("SegCopy.txt", delim = "\t")
seg <- segfix #Can use anyone, all commands below are the same.
c <- seg[,c(4:94)] #just samples
c2 <- seg[,c(1:3)]
c <- subset(c, select = hcc.good.1.ctc) #Getting just the columns matching the samples we want
names(c) <- better
segs <- as.data.frame(c(c2, c))
df <- makeGRangesFromDataFrame(segs, keep.extra.columns = F)
df

#Getting gene panel
good <- c(1:22)
#good <- paste("chr", good, sep = "")
#genes <- genes %>% filter(chromosome %in% good)
colnames(genes) <- Cs(chromosome, start, end, cytoband)
genes$loc <- round((as.numeric(genes$start) + as.numeric(genes$end)) / 2, 0)
genes <- genes[, c(1,5,4)]
genes$start <- genes$loc
colnames(genes)[2] <- "end"
genes <- genes[,c(1,4,2,3)]
genes$chromosome <- paste("chr", genes$chromosome, sep = "")
#genes$cytoband = paste("cb_", genes$cytoband, sep = "")
#
genes.gr <- makeGRangesFromDataFrame(genes, keep.extra.columns = T)
genes.gr
#Getting the genes into the format of the Segcopy/fixed/norm, etc so that there are 6194 rows and only saving the gene symbol names
symInCnv <- splitColumnByOverlap(genes.gr,df,column = "cytoband")
y <- unstrsplit(symInCnv, sep = ", ")
y <- as.data.frame(y) #6194 row and 1 character string of gene names in the appropriate places!
colnames(y) <- "cytoband"

seg.genes <- cbind(segs, y) #DOUBLE CHECK THAT THE NUMBER OF COLUMNS IS RIGHT
seg.genes <- seg.genes[,c(1:3,24,4:23)]
genes <- seg.genes %>% filter(cytoband %not in% "") #Keep only CNVs with known cancer genes
genes$cytoband = gsub("\\,.*", "", genes$cytoband)
c <- genes[,c(4:24)] #DOUBLE CHECK THE NUMBER OF COLUMNS
row.names(c) <- c$cytoband
c <- c[,c(2:21)]
#c$H195.CTC = c$H195.CTC - 1.6 #Normalize to CN2
#c$H169.TUMOR = c$H169.TUMOR - 3 #Morn to CN2
c.names <- names(c)
c <- as.matrix(c)
colMeans(c)
#c <- c - 1

heatmap(c, scale = "none")
col<- colorRampPalette(c("blue", "white", "red"))(256)
heatmap(c, scale = "none", col =  col)
colMeans(c)
x <- cor(c, method = "spearman")
corrplot(x, method = "ellipse", diag = T, addrect = 11)
x <- prcomp(t(c))
ggplot2::autoplot(x, label = T)
heatmap.2(c, scale = "row", col = bluered(256), trace = "none", dendrogram = "column", Rowv = T, density.info = "none", key = F, ColSideColors = c(rep("#E64B35FF", 2), rep("#4DBBD5FF", 2), rep("#00A087FF", 2), rep("#3C5488FF", 2), rep("#F39B7FFF", 2), rep("#8491B4FF", 2), rep("#91D1C2FF", 2), rep("#DC0000FF", 2), rep("#7E6148FF", 2), rep("#B09C85FF", 2)),colsep = c(2,4,6,8,10,11,13,15,17,18), sepcolor = "white")
coords <- locator(1)


# * * * * * * Plotting for paper ------------------------------------------
#Plotting for paper
pdf("output.pdf",width=10,height=20)
plot()
heatmap.2(c, scale = "none", col = bluered(256), trace = "none", dendrogram = "column", Rowv = F, density.info = "none", key = F, ColSideColors = c(rep("#E64B35FF", 2), rep("#4DBBD5FF", 2), rep("#00A087FF", 2), rep("#3C5488FF", 2), rep("#F39B7FFF", 2), rep("#8491B4FF", 2), rep("#91D1C2FF", 2), rep("#DC0000FF", 2), rep("#7E6148FF", 2), rep("#B09C85FF", 2)))
plot(1, type="n", axes=F, xlab="", ylab="")
#Older version
legend("center", xpd=TRUE,     
       legend = c.names,
       col = c(rep("#E64B35FF", 2), rep("#4DBBD5FF", 2), rep("#00A087FF", 2), rep("#3C5488FF", 2), rep("#F39B7FFF", 2), rep("#8491B4FF", 2), rep("#91D1C2FF", 2), rep("#DC0000FF", 2), rep("#7E6148FF", 2), rep("#B09C85FF", 2)), 
       lty= 1,             
       lwd = 5,           
       cex=.7
)

#Newer version
plot(1, type="n", axes=F, xlab="", ylab="",family = "arial")
legend("center", xpd=TRUE,     
       legend = Cs(H200, H199,H195,H187,H185,H175,H172,H169,H167,H165),
       col = mypal, 
       lty= 1,
       lwd = 5,           
       cex=.7
)


#with 5 columns
quartz()
plot(1, type="n", axes=F, xlab="", ylab="",family = "arial")
legend("center", xpd=TRUE,     
       legend = Cs(H200, H199,H195,H187,H185,H175,H172,H169,H167,H165),
       col = mypal, 
       lty= 1, ncol = 5, lwd = 5
)
graphics.off() #Turns off quartz windows
legend
dev.off()
x = -6, y = 1

###Getting supplemental table of top genes:
x <- top.bands[,c(1:5,7)]
names(x) <- c("cytoband", "chromosome", "start", "end", "giemsa stain", "cnv type")
write.table(x, "supplemental.table.4.txt", quote = F, row.names = F, sep = '\t') 




# * * CTC VS. TUMOR *****-------------------------------------------------------
setwd("~/NGS/HCC.NGS.CNV/All.samples.250kbp.fix.iod.masked")
setwd("~/NGS/HCC.NGS.CNV/clustering/all.var.250kbp.normreadct.masked")
setwd("~/NGS/HCC.NGS.CNV/clustering/all.1MB.fix.normreadct")
segnorm <- read_delim("SegNorm.txt", delim = "\t")
segfix <- read_delim("SegFixed.txt", delim = "\t")
segcopy <- read_delim("SegCopy.txt", delim = "\t")
seg <- segfix #Can use anyone, all commands below are the same.
c <- seg[,c(4:94)] #just samples
c2 <- seg[,c(1:3)]
c <- subset(c, select = hcc.good.1.ctc) #Getting just the columns matching the samples we want
names(c) <- better
segs <- as.data.frame(c(c2, c))
df <- makeGRangesFromDataFrame(segs, keep.extra.columns = F)
df
genes <- new.genes
#good <- c(1:22)
#good <- paste("chr", good, sep = "")
#genes <- genes %>% filter(chromosome %in% good)
#colnames(genes) <- Cs(chromosome, start, end, cytoband)
#genes$loc <- round((as.numeric(genes$start) + as.numeric(genes$end)) / 2, 0)
#genes <- genes[, c(1,5,4)]
#genes$start <- genes$loc
#colnames(genes)[2] <- "end"
#genes <- genes[,c(1,4,2,3)]
#genes$cytoband = paste("cb_", genes$cytoband, sep = "")
#
genes.gr <- makeGRangesFromDataFrame(genes, keep.extra.columns = T)
genes.gr
#Getting the genes into the format of the Segcopy/fixed/norm, etc so that there are 6194 rows and only saving the gene symbol names
symInCnv <- splitColumnByOverlap(genes.gr,df,column = "cytoband")
y <- unstrsplit(symInCnv, sep = ", ")
y <- as.data.frame(y) #6194 row and 1 character string of gene names in the appropriate places!
colnames(y) <- "cytoband"

seg.genes <- cbind(segs, y) #DOUBLE CHECK THAT THE NUMBER OF COLUMNS IS RIGHT
seg.genes <- seg.genes[,c(1:3,24,4:23)]
genes <- seg.genes %>% filter(cytoband %not in% "") #Keep only CNVs with known cancer genes
genes$cytoband = gsub("\\,.*", "", genes$cytoband)
c <- genes[,c(4:24)] #DOUBLE CHECK THE NUMBER OF COLUMNS
row.names(c) <- c$cytoband
c <- c[,c(2:21)]
#c$H195.CTC = c$H195.CTC - 1.6 #Normalize to CN2
#c$H169.TUMOR = c$H169.TUMOR - 3 #Morn to CN2
c.names <- names(c)
c <- as.matrix(c)
c = c[,c(1:14,17:20)] #this gets rid of H169
colMeans(c)
#c <- c - 1

heatmap(c, scale = "none")
col<- colorRampPalette(c("blue", "white", "red"))(256)
heatmap(c, scale = "none", col =  col)
colMeans(c)
x <- cor(c, method = "spearman")
corrplot(x, method = "ellipse", diag = T, addrect = 11)
x <- prcomp(t(c))
ggplot2::autoplot(x, label = T)
heatmap.2(c, scale = "none", col = bluered(256), trace = "none", dendrogram = "column", Rowv = F, density.info = "none", key = T, ColSideColors = c(rep("#E64B35FF", 2), rep("#4DBBD5FF", 2), rep("#00A087FF", 2), rep("#3C5488FF", 2), rep("#F39B7FFF", 2), rep("#8491B4FF", 2), rep("#91D1C2FF", 2), rep("#DC0000FF", 2), rep("#7E6148FF", 2), rep("#B09C85FF", 2)))
coords <- locator(1)


# * * * * * * Plotting for paper ------------------------------------------

##NEW 20170728
#Plotting for paper
setwd("~/NGS/HCC.NGS.CNV/example")
pdf("output.pdf",width=10,height=20)
plot()
heatmap.2(c, scale = "none", col = bluered(256), trace = "none", dendrogram = "column", Rowv = F, density.info = "none", key = F, ColSideColors = c(rep("#E64B35FF", 2), rep("#4DBBD5FF", 2), rep("#00A087FF", 2), rep("#3C5488FF", 2), rep("#F39B7FFF", 2), rep("#8491B4FF", 2), rep("#91D1C2FF", 2), rep("#7E6148FF", 2), rep("#B09C85FF", 2)))
plot(1, type="n", axes=F, xlab="", ylab="")
#Older version
legend("center", xpd=TRUE,     
       legend = c.names,
       col = c(rep("#E64B35FF", 2), rep("#4DBBD5FF", 2), rep("#00A087FF", 2), rep("#3C5488FF", 2), rep("#F39B7FFF", 2), rep("#8491B4FF", 2), rep("#91D1C2FF", 2), rep("#DC0000FF", 2), rep("#7E6148FF", 2), rep("#B09C85FF", 2)), 
       lty= 1,             
       lwd = 5,           
       cex=.7
)

#Newer version
plot(1, type="n", axes=F, xlab="", ylab="",family = "arial")
legend("center", xpd=TRUE,     
       legend = Cs(H200, H199,H195,H187,H185,H175,H172,H169,H167,H165),
       col = mypal, 
       lty= 1,
       lwd = 5,           
       cex=.7
)

legend
dev.off()
x = -6, y = 1

#Plotting for paper
pdf("output.pdf",width=10,height=20)
plot()
heatmap.2(c, scale = "none", col = bluered(256), trace = "none", dendrogram = "column", Rowv = F, density.info = "none", key = F, ColSideColors = c(rep("#E64B35FF", 2), rep("#4DBBD5FF", 2), rep("#00A087FF", 2), rep("#3C5488FF", 2), rep("#F39B7FFF", 2), rep("#8491B4FF", 2), rep("#91D1C2FF", 2), rep("#DC0000FF", 2), rep("#7E6148FF", 2), rep("#B09C85FF", 2)))
plot(1, type="n", axes=F, xlab="", ylab="")
#Older version
legend("center", xpd=TRUE,     
       legend = c.names,
       col = c(rep("#E64B35FF", 2), rep("#4DBBD5FF", 2), rep("#00A087FF", 2), rep("#3C5488FF", 2), rep("#F39B7FFF", 2), rep("#8491B4FF", 2), rep("#91D1C2FF", 2), rep("#DC0000FF", 2), rep("#7E6148FF", 2), rep("#B09C85FF", 2)), 
       lty= 1,             
       lwd = 5,           
       cex=.7
)

#Newer version
plot(1, type="n", axes=F, xlab="", ylab="",family = "arial")
legend("center", xpd=TRUE,     
       legend = Cs(H200, H199,H195,H187,H185,H175,H172,H169,H167,H165),
       col = mypal, 
       lty= 1,
       lwd = 5,           
       cex=.7
)

legend
dev.off()
x = -6, y = 1

###Getting supplemental table of top genes:
x <- top.bands[,c(1:5,7)]
names(x) <- c("cytoband", "chromosome", "start", "end", "giemsa stain", "cnv type")
write.table(x, "supplemental.table.4.txt", quote = F, row.names = F, sep = '\t') 

#Removing H169
c2 <- c[,c(1:14,17:20)]
heatmap(c2, scale = "none")
col<- colorRampPalette(c("red", "white", "blue"))(256)
heatmap(c2, scale = "none", col =  col)
x <- prcomp(t(c2)) #Helps a little but still fairly seperate
ggplot2::autoplot(x, label = T)
heatmap.2(c2, scale = "none", col = bluered(7), breaks = c(0,1.5,2.5,3.5,4.5,5.5,6.5,7.5), trace = "none", density.info = "none")
#Removing H195 (PRE FIGURING OUT I WAS USING THE WRONG VERSION OF H195 CTC) and H200
c3 = c2[,c(3:4,7:18)]
x <- prcomp(t(c3)) #Helps a little but still fairly seperate
ggplot2::autoplot(x, label = T)
x <- cluster::daisy(t(c))

#Example heatmap:
df <- as.matrix(scale(mtcars))
heatmap(df, scale = "none")
col<- colorRampPalette(c("red", "white", "blue"))(256)
heatmap(scale(as.matrix(mtcars)), scale = "none", col =  col)
#write.table(top.bands, file = "top.hcc.cytobands.115.txt", quote = F, sep = "\t", row.names = F, col.names = T)
#
#



# HCC CTC SCNA Supplemental figures ---------------------------------------


setwd("~/NGS/HCC.NGS.CNV/clustering/all.var.250kbp.normreadct.masked")
segnorm <- read_delim("SegNorm.txt", delim = "\t")
seg <- segnorm #Can use anyone, all commands below are the same.
c <- seg[,c(4:94)] #just samples
c2 <- seg[,c(1:3)]
c <- subset(c, select = h195.samples) #Getting just the columns matching the samples we want
names(c) <- h195.names
c = as.matrix(c)
b = c[,c(1:6,8:10)]
x = prcomp(b)
y = get_eigenvalue(x) #60% and 17%
z = as.data.frame(x$rotation)
h195.names.good <- c("Blood", "Normal Liver", "Tumor", "CTC 1", "CTC 2", "CTC 3", "CTC 5", "CTC 6", "CTC 7")
z$sample = h195.names.good
keep = c("sample", "PC1", "PC3")
pca.vals = subset(z, select = keep)
pca.vals$PC1 = pca.vals$PC1 * 0.6
pca.vals$PC2 = pca.vals$PC2 * 0.17
plot(pca.vals$PC1, pca.vals$PC2)
text(pca.vals$PC1, pca.vals$PC2, labels=pca.vals$sample, cex= 0.7, pos = 3)

p = ggplot(pca.vals, aes(x = PC1, y = PC2)) +
  geom_point() +
  coord_fixed(ratio = 1)
p
