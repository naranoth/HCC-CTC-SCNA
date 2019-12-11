#For HCC CTC SCNA classification paper


# Initial pre-processing ------------------------------------------------------
cb.cnv <- read.delim("cb.cnv.all.samples.20181022.txt", sep = "\t") #11024 obs, 1049 vars (1043 genes + other vars)

all.cbs = read.delim("all_cytobands_unique.txt", sep = "\t", header = T) # From Zack et al, Nature Genetics, DOI: 10.1038/ng.2760)
all.cbs = as.character(all.cbs$all.cbs)

#Getting list of cytobands to be filtered out so that just Zack et al remaining:
cbs = colnames(cb.cnv[7:572])
all.cbs.cb = paste("cb", all.cbs, sep = ".")
top.cbs = intersect(cbs, all.cbs.cb) #116 cytobands
#List of columns to keep, first 6 plus the 116 intersecting SCNAs 
col.keep = c("sample.name", "origin", 'dataset', 'rtsne1', 'rtsne2', 'cin', top.cbs) #274


# Models ------------------------------------------------------------------
###Final DF for model, all cancer subtypes:
top.cnvs = subset(cb.cnv, select = col.keep) 
#Baseline 31 class problem with RF:
#Overall Statistics
#Accuracy : 0.5816          
#95% CI : (0.5601, 0.6028)
#No Information Rate : 0.1059          
#P-Value [Acc > NIR] : < 2.2e-16      

###Final DF model, 21 cancer types in 11 classes:
top.cnvs = subset(cb.cnv, select = col.keep) #11024 obs of 122 vars for global only; 274 vars with disease-specific top cytobands included.
keep.origin = c('ACC', 'CHOL', 'COADREAD', 'GBM', 'KICH', 'KIRC', 'KIRP', 'LAML', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'OV', 'PAAD', 'PCPG', 'PRAD', 'STAD', 'TGCT', 'THCA', 'THYM', 'UVM', 'COADREAD_LP', 'LGG_LP', 'STAD_LP') #2019: added in LP ones
top.cnvs = top.cnvs %>% filter(origin %in% keep.origin) 
top.cnvs$origin = as.factor(top.cnvs$origin)
top.cnvs$origin = plyr::mapvalues(top.cnvs$origin, from = c('ACC', 'CHOL', 'COADREAD', 'GBM', 'KICH', 'KIRC', 'KIRP', 'LAML', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'OV', 'PAAD', 'PCPG', 'PRAD', 'STAD', 'TGCT', 'THCA', 'THYM', 'UVM', 'COADREAD_LP', 'LGG_LP', 'STAD_LP'), to = c('ADRENAL', 'LIVER', 'GI', 'NERVOUS', 'KIDNEY', 'KIDNEY', 'KIDNEY', 'CNLOW', 'NERVOUS', 'LIVER', 'LUNG', 'LUNG', 'OVARY', 'GI', 'PHEO', 'CNLOW', 'GI', 'MALE', 'CNLOW', 'CNLOW', 'MELANOMA', 'GI', 'NERVOUS', 'GI'))
top.cnvs$origin = droplevels(top.cnvs$origin) #7328 Samples in 20190124
#RF 500 trees on 20190131
#Overall Statistics
#Accuracy : 0.8321          
#95% CI : (0.7812, 0.8752)
#No Information Rate : 0.8015          
#P-Value [Acc > NIR] : 0.1214          
#Kappa : 0.6148

#plot(rf.imp) ... it looks like 150 trees is about what i need. After that the OOB error rate basically levels off.

#RF 100 trees new 20190123
#Overall Statistics
#Accuracy : 0.783           
#95% CI : (0.7605, 0.8043)
#No Information Rate : 0.1963          
#P-Value [Acc > NIR] : < 2.2e-16       
#Kappa : 0.7482     



# Random Forest -----------------------------------------------------------

#Using pan-genome markers as well:
top.cnvs2 = top.cnvs
tcga = top.cnvs2[top.cnvs2$dataset == 'tcga',] #10478
tcga$origin = as.factor(tcga$origin)
table(tcga$origin)
ctc = top.cnvs2[top.cnvs2$dataset == 'ctc',] #63
ctc = ctc %>% filter(sample.name %in% keep.ctcs) #44, these are the ones for the paper


#Using just the cbs for classification
top.cnvs2 = top.cnvs[,c(1,2,3,7:274)]
tcga = top.cnvs2[top.cnvs2$dataset == 'tcga',] #10478
tcga$origin = as.factor(tcga$origin)
table(tcga$origin)
ctc = top.cnvs2[top.cnvs2$dataset == 'ctc',] #63
ctc = ctc %>% filter(sample.name %in% keep.ctcs) #45, these are the ones for the paper

#Decide which column to use as classifier: 2 = TCGA origin name
tcga_df = tcga[,c(2,4:length(tcga))]
lp_df = lp[,c(2,4:length(tcga))]
ctc_df = ctc[,c(2,4:length(tcga))]

tcga_df$origin = as.factor(tcga_df$origin)
smp_size = floor(0.8 * nrow(tcga_df))
set.seed(123)
train_ind = sample(seq_len(nrow(tcga_df)), size = smp_size)
train = tcga_df[train_ind,] #3225
train = upSample(x = train, y = train$origin) ###Upsampling training set to get equal partitions
train = train[,-ncol(train)] #removing the extra column that upsample has added.
test = tcga_df[-train_ind,] #807
train$origin = droplevels(train$origin) 
test$origin = droplevels(test$origin) 

#Using RF
rf.imp = randomForest(x = train[,-1], y = train$origin, data = train, importance = T, do.trace = T, ntree = 500)
pred = predict(rf.imp, newdata = test)
confusionMatrix(pred, test$origin)



# CTC classification ------------------------------------------------------

ctc.pred = predict(rf.imp, newdata = ctc_df)
ctc_df$origin = as.factor(ctc_df$origin)
ctc.prediction = as.data.frame.matrix(table(ctc.pred, ctc_df$origin))
samps = as.character(ctc.pred)

#Just the top CTCs: keep.ctcs (the new one) = 44 ctcs from labels section above
sample.pred = data.frame(keep.ctcs, samps)
sample.predictions = data.frame(ctc$sample.name, samps)

#Getting CTC probabilities from RF Importance
preds.ctc = as.data.frame(predict(rf.imp, newdata = ctc_df, "prob"))
preds.ctc$sample = keep.ctcs


#Getting the HCC and Lung CTC Samples in order
hcc.ctc = preds.ctc[c(1:15),]
hcc.ctc = hcc.ctc %>% arrange(-LIVER) 
hcc.ord = hcc.ctc$sample #Get the sample names in descending order of LIVER value
lung.ctc = preds.ctc[16:44,]
lung.ctc = lung.ctc %>% arrange(-LUNG)
lung.ord = lung.ctc$sample

#HCC
hcc.pred = gather(data = hcc.ctc, key = 'origin', value = 'votes', -sample)
hcc.pred = hcc.pred %>% mutate_if(is.character, as.factor)
hcc.pred$origin = relevel(hcc.pred$origin, "LIVER")
hcc.pred$origin = fct_rev(hcc.pred$origin)
hcc.pred$sample = factor(hcc.pred$sample, levels = hcc.ord)
library(ggthemes)
p = ggplot(hcc.pred, aes(x = sample,y = votes, fill = origin)) +
  geom_bar(stat = "identity") +
  theme_calc() + 
  scale_fill_stata()
p


#LUNG
lung.pred = gather(data = lung.ctc, key = 'origin', value = 'votes', -sample)
lung.pred = lung.pred %>% mutate_if(is.character, as.factor)
lung.pred$origin = relevel(lung.pred$origin, "LUNG")
lung.pred$origin = fct_rev(lung.pred$origin)
lung.pred$sample = factor(lung.pred$sample, levels = lung.ord)
library(ggthemes)
p = ggplot(lung.pred, aes(x = sample,y = votes, fill = origin)) +
  geom_bar(stat = "identity") +
  theme_calc() + 
  scale_fill_stata()
p


# Top CTC only cytobands --------------------------------------------------

cb.cnv <- read.delim("cb.cnv.all.samples.20181022.txt", sep = "\t") 
scnas = read.delim("freq_scna.txt", sep = "\t",header = T)
all.cbs = read.delim("all_cytobands_unique.txt", sep = "\t", header = T)
all.cbs = as.character(all.cbs$all.cbs)

#Getting just the CTC and Tumors
samps = cb.cnv[1:35,]
samps = samps %>% separate(sample.name, c("patient.id", "sample"), "\\.",remove = F)

#Getting rid of excess columns
samps = samps[,c(1,2,3,4,9:length(samps))]

#Getting rid of CCA samples
samps = samps %>% filter(origin == "LIHC")

#Getting just 1 CTC per tumor to keep things consistent:
good = Cs(H165.CTC, H167.CTC, H172.CTC, H175.CTC, H185.CTC, H187.CTC, H195.CTC.A, H199.CTC.B, H200.CTC, H165.TUMOR, H167.TUMOR, H172.TUMOR, H175.TUMOR, H185.TUMOR, H187.TUMOR, H195.TUMOR, H199.TUMOR, H200.TUMOR)
samps = samps %>% filter(sample.name %in% good)

#Going to 3 level instead of 5 level
samps = samps %>% mutate_if(is.integer, as.character)
samps[,c(5:length(samps))] = sapply(X = samps[,c(5:length(samps))],function(x) fct_collapse(x, amp = c("1", "2"), neutral = "0", del = c("-2", "-1")))
samps = samps %>% mutate_if(is.character, as.factor)

#Separating tumor and CTC counts and getting them for each possible cytoband:value combination
tumor = samps %>% filter(sample == "TUMOR")
tumor = tumor[,c(3,5:length(tumor))]
tumor2 = tumor %>% gather(key = variable, value = count, cb.10p11.21:length(tumor))
tumor2$count = as.factor(tumor2$count)
tumor3 = as.data.frame(table(tumor2$variable, tumor2$count))
colnames(tumor3) = c("cytoband", "value", "count")
tumor3$freq = tumor3$count / 9
tumor3$type = "TUMOR"

#For CTCs, same thing
ctc = samps %>% filter(sample == "CTC")
ctc = ctc[,c(3,5:length(ctc))]
ctc2 = ctc %>% gather(key = variable, value = count, cb.10p11.21:length(ctc))
ctc2$count = as.factor(ctc2$count)
ctc3 = as.data.frame(table(ctc2$variable, ctc2$count))
colnames(ctc3) = c("cytoband", "value", "count")
ctc3$freq = ctc3$count / 9
ctc3$type = "CTC"

#Recombining the two for comparison:
combined = rbind(tumor3, ctc3) 
combined = combined[,-3]
combined = combined %>% spread(key = type, value = freq)
na.zero <- function (x) {
  x[is.na(x)] <- 0
  return(x)
}
combined = na.zero(combined)

#Getting rid of neutral because we don't care about it
combined = combined %>% filter(value == c("amp", "del"))

#Getting the absolute value of the difference between the two sample types
combined$diff = abs(combined$CTC - combined$TUMOR)
combined$ctc.diff = combined$CTC - combined$TUMOR

#Limiting to cytobands with over 20% difference
ctc.unique = combined %>% filter(ctc.diff > 0.2) #5 observations



# t-SNE space -------------------------------------------------------------

pca.cnv <- read.delim("pca.cnv.0.95cutoff.txt", sep = "\t") #11024 obs, 1049 vars (1043 genes + other vars)
###Trying t-SNE space
cnv.mat = pca.cnv[,c(7:1049)]
rownames(cnv.mat) = pca.cnv$sample.name
cnv.mat$samp = seq(0,1,length.out = 11024)
cnv.mat = as.matrix(cnv.mat)

###Double check for duplicates
head(cnv.mat[1:5,1:5])
table(duplicated(colnames(cnv.mat)))
unique.colnames = unique(colnames(cnv.mat)) #1044, cin + 1043 genes

x = as.data.frame(duplicated(cnv.mat))
row.names(x) = row.names(cnv.mat)
x$sample.name = row.names(x)
dups = left_join(x,cnv.mat)
table(x)


cnv.tsne = Rtsne(cnv.mat, theta = 0.4, verbose = T)

sne.vars = as.data.frame(cnv.tsne$Y)
