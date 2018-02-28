library("DESeq2") # for differntial gene expression analysis
library("arrayQualityMetrics") # to call outliers
library("Biobase") #for ExpressionSet
library("pheatmap") # for sample heatmap
library("VennDiagram") # for Venn diagram of #s of DEGs

setwd("~/Documents/davieslab/fgb/deseq2/")


# Load data ---------------------------------------------------------------
countdata = read.table("allcounts_FGB_OfavB1_jan2018.txt",header=TRUE,row.names = 1) 
head(countdata) 
length(countdata[,1])
# 72134 isogroups mapped

# make namers pretty
names(countdata) <- gsub("_comb.fastq.trim.sam.counts", "", names(countdata))
names(countdata)

# split coral and symbiont ---------------
head(countdata)
tail(countdata)

sym.rows <- grep("symbB1", row.names(countdata))
coral.count <- countdata[-sym.rows,]
head(coral.count)
tail(coral.count)

sym.count <- countdata[sym.rows,]
head(sym.count)
tail(sym.count)

# rows add up to total? didn't miss anything?
nrow(coral.count)+nrow(sym.count) == nrow(countdata)
#yep

# CORAL ANALYSIS -------------------------------------------
#----------Total counts?
totalCounts = colSums(coral.count)
head(totalCounts)
totalCounts
min(totalCounts) #53,148
mean(totalCounts) #271,631.6
max(totalCounts)  #464,432

# make conditions table --------
names(coral.count)

species <- ifelse( grepl("FR",names(coral.count)), "FR", "OF")

lesion <- ifelse( grepl("AH", names(coral.count)), "AH", ifelse( grepl("AL", names(coral.count)), "AL", "U"))

bank <- ifelse( grepl("W",names(coral.count)), "west", "east")

matches <- regmatches(names(coral.count), gregexpr("[[:digit:]]+", names(coral.count)))
colony <- as.numeric(unlist(matches))
colony <- paste(bank, colony, sep="_")
colony

coldata = as.data.frame(cbind(species, lesion, bank, colony))
row.names(coldata) = names(coral.count)
coldata

# Construct data object ---------------------------------------------------

dds <- DESeqDataSetFromMatrix(
  countData = coral.count,
  colData = coldata,
  design = ~ species + lesion + bank)

save(coldata, coral.count, dds, file="ddsCoral.Rdata")

# Call outliers -----------------------------------------------------------

rld = rlogTransformation(dds, blind=TRUE)
e = ExpressionSet(assay(rld), AnnotatedDataFrame(as.data.frame(colData(rld))))
arrayQualityMetrics(e, intgroup=c("species", "lesion"), force=T)

# Set base mean minimum ---------------------------------------------------
means = apply(coral.count,1,mean)
table(means>3)
# FALSE  TRUE 
# 13163  6992 

means3 = names(means[means>3])
head(means3)
length(means3)
#6992

coral.countFilt = coral.count[row.names(coral.count) %in% means3,]
head(coral.countFilt)

totalCountsFilt = colSums(coral.countFilt)
totalCountsFilt

min(totalCountsFilt) #52276
max(totalCountsFilt) #448716
mean(totalCountsFilt) #263389.3

# check sample order
test = cbind(names(coral.countFilt),as.vector(row.names(coldata)))
test

# Reconstruct data object (filtered) ------------------------------

ddsFilt <- DESeqDataSetFromMatrix(
  countData = coral.countFilt,
  colData = coldata,
  design = ~ species + lesion + bank)

# DESeq -------------------------------------------------------------------

#-------------DESeq pipeline in one step: makes large DESeqDataSet
deds <- DESeq(ddsFilt)
# -- replacing outliers and refitting for 40 genes
# -- DESeq argument 'minReplicatesForReplace' = 7 

# Load and save
save(coldata, coral.countFilt, ddsFilt, deds, file="ddsFilt.Rdata")
load("ddsFilt.Rdata")

#---Results

resSpecies = results(deds, independentFiltering = F, contrast=c("species","FR","OF"))
resSpecies

resBank = results(deds, independentFiltering = F, contrast=c("bank", "west", "east"))
resBank

resALvAH = results(deds, independentFiltering = F, contrast=c("lesion", "AL", "AH"))
resALvAH

resALvU = results(deds, independentFiltering = F, contrast=c("lesion", "AL", "U"))
resALvU

resAHvU = results(deds, independentFiltering = F, contrast=c("lesion", "AH", "U"))
resAHvU


# how many genes pass multiplicity-corrected 0.1 FDR cutoff?
table(resSpecies$padj < 0.1)
# FALSE  TRUE 
# 5774  1218 

table(resBank$padj < 0.1)
# FALSE  TRUE 
# 6937    55

table(resALvAH$padj < 0.1)
# FALSE  TRUE 
# 6835   157 

table(resALvU$padj < 0.1)
# FALSE  TRUE 
# 6624   368 

table(resAHvU$padj < 0.1)
# FALSE 
# 6992 


# Extract log-fold changes --------
table(abs(resSpecies$log2FoldChange)>1.5)
# FALSE  TRUE 
#   6915    77
write.csv( as.data.frame(resSpecies), file="resSpecies.csv" ) 

table(abs(resBank$log2FoldChange)>1.5)
# FALSE  TRUE 
# 6965    27 
write.csv( as.data.frame(resBank), file="resBank.csv" ) 

table(abs(resALvAH$log2FoldChange)>1.5)
# FALSE  TRUE 
# 6950    42 
write.csv( as.data.frame(resALvAH), file="resALvAH.csv" ) 

table(abs(resAHvU$log2FoldChange)>1.5)
# FALSE  TRUE 
# 6974    18 
write.csv( as.data.frame(resAHvU), file="resAHvU.csv" ) 

table(abs(resALvU$log2FoldChange)>1.5)
# FALSE  TRUE 
#  6923    69 
write.csv( as.data.frame(resALvU), file="resALvU.csv" ) 

# new rld ------

rld = rlogTransformation(ddsFilt, blind=TRUE)

# Save/Load Data ----------------------------------------------------------

# save(rld, coldata, coral.countFilt, ddsFilt, deds, resSpecies, resBank, resAHvU, resALvAH, resALvU, rld, file="results.Rdata")
load("results.Rdata")

# Explore with plots ------------------------------------------------------

#Sample distance heatmap
pheatmap(cor(assay(rld)),border_color=NA, main="SampleHeatmap")

# Diagnostics -------------------------------------------------------------

#Dispersions plot
plotDispEsts(deds, main="Dispersion Plot Response")

#MA plot
plotMA(resSpecies, ylim = c(-1, 1), main="MA Plot Species") 
plotMA(resBank, ylim = c(-1, 1), main="MA Plot Bank") 

plotMA(resALvAH, ylim = c(-1, 1), main="MA Plot AL v AH") 
plotMA(resALvU, ylim = c(-1, 1), main="MA Plot AL v U") 
plotMA(resAHvU, ylim = c(-1, 1), main="MA Plot AH v U") 

# Write results for making heatmaps ---------------------------------------

###--------------Get pvals
head(resSpecies)
valsSpecies = cbind(resSpecies$pvalue, resSpecies$padj)
head(valsSpecies)
colnames(valsSpecies)=c("pval.s", "padj.s")
length(valsSpecies[,1])
table(complete.cases(valsSpecies))

head(resBank)
valsBank = cbind(resBank$pvalue, resBank$padj)
head(valsBank)
colnames(valsBank)=c("pval.b", "padj.b")
length(valsBank[,1])
table(complete.cases(valsBank))

head(resALvAH)
valsALvAH = cbind(resALvAH$pvalue, resALvAH$padj)
head(valsALvAH)
colnames(valsALvAH)=c("pval.alah", "padj.alah")
length(valsALvAH[,1])
table(complete.cases(valsALvAH))

head(resALvU)
valsALvU = cbind(resALvU$pvalue, resALvU$padj)
head(valsALvU)
colnames(valsALvU)=c("pval.alu", "padj.alu")
length(valsALvU[,1])
table(complete.cases(valsALvU))

head(resAHvU)
valsAHvU = cbind(resAHvU$pvalue, resAHvU$padj)
head(valsAHvU)
colnames(valsAHvU)=c("pval.ahu", "padj.ahu")
length(valsAHvU[,1])
table(complete.cases(valsAHvU))

#Make rlogdata and pvals table
rldpvals=cbind(assay(rld),valsSpecies, valsBank, valsALvAH, valsALvU, valsAHvU)
head(rldpvals)
dim(rldpvals)
# 6992   86
table(complete.cases(rldpvals))
  
write.csv(rldpvals, "25jan2018_fgb_bm3_wald_RLDandPVALS.csv", quote=F)

# Venn diagram ---------
# rldpvals=read.csv("25jan2018_fgb_bm3_wald_RLDandPVALS.csv")
rldpvals=as.data.frame(rldpvals)

alah=row.names(rldpvals[rldpvals$padj.alah<0.1 & !is.na(rldpvals$padj.alah),])
alu=row.names(rldpvals[rldpvals$padj.alu<0.1 & !is.na(rldpvals$padj.alu),])
# ahu=row.names(rldpvals[rldpvals$padj.ahu<0.1 & !is.na(rldpvals$padj.ahu),])               # NO DEGs

degs10 = union(alah, alu)

length(degs10)


candidates=list("ALvAH"=alah, "ALvU"=alu)

prettyvenn=venn.diagram(
  x = candidates,
  filename=NULL,
  col = "transparent",
  fill = c("coral2", "forestgreen"),
  alpha = 0.5,
  label.col = c("darkred", "white", "darkgreen"),
  cex = 2.5,
  fontfamily = "sans",
  fontface = "bold",
  cat.default.pos = "text",
  cat.col = c("darkred", "darkgreen"),
  cat.cex = 2.5,
  cat.fontfamily = "sans",
  cat.dist = c(0.08, 0.08),
  cat.pos = 1
);
grid.draw(prettyvenn)

# what are the genes?
alah_unique <- alah[!alah %in% alu]
alah_unique
alah_unique_rld <- rldpvals[alah_unique,]
alah_unique_rld$X <- sub("c", "isogroup", alah_unique_rld$X)

head(gg)

alah_unique_annot <- merge(alah_unique_rld, gg, by = 1)
alah_unique_annot$V2

# Write results for GO/KOG analysis -------------------------------------------

# by -log p-value
logs=data.frame(cbind("gene"=row.names(resSpecies),"logP"=round(-log(resSpecies$pvalue+1e-10,10),1)))
logs$logP=as.numeric(as.character(logs$logP))
sign=rep(1,nrow(logs))
sign[resSpecies$log2FoldChange<0]=-1  ##change to correct model
table(sign)
# -1    1 
# 3096 3896 
logs$logP=logs$logP*sign
logs$gene=gsub("c", "isogroup", logs$gene)
write.table(logs,quote=F,row.names=F,file="GO_Species_logP.csv",sep=",")

logs=data.frame(cbind("gene"=row.names(resBank),"logP"=round(-log(resBank$pvalue+1e-10,10),1)))
logs$logP=as.numeric(as.character(logs$logP))
sign=rep(1,nrow(logs))
sign[resBank$log2FoldChange<0]=-1  ##change to correct model
table(sign)
# -1    1 
# 3381 3611 
logs$logP=logs$logP*sign
logs$gene=gsub("c", "isogroup", logs$gene)
write.table(logs,quote=F,row.names=F,file="GO_Bank_logP.csv",sep=",")

logs=data.frame(cbind("gene"=row.names(resAHvU),"logP"=round(-log(resAHvU$pvalue+1e-10,10),1)))
logs$logP=as.numeric(as.character(logs$logP))
sign=rep(1,nrow(logs))
sign[resAHvU$log2FoldChange<0]=-1  ##change to correct model
table(sign)
# -1    1 
# 3315 3677 
logs$logP=logs$logP*sign
logs$gene=gsub("c", "isogroup", logs$gene)
write.table(logs,quote=F,row.names=F,file="GO_AHvU_logP.csv",sep=",")

logs=data.frame(cbind("gene"=row.names(resALvAH),"logP"=round(-log(resALvAH$pvalue+1e-10,10),1)))
logs$logP=as.numeric(as.character(logs$logP))
sign=rep(1,nrow(logs))
sign[resALvAH$log2FoldChange<0]=-1  ##change to correct model
table(sign)
# -1    1 
# 3895 3097 
logs$logP=logs$logP*sign
logs$gene=gsub("c", "isogroup", logs$gene)
write.table(logs,quote=F,row.names=F,file="GO_ALvAH_logP.csv",sep=",")

logs=data.frame(cbind("gene"=row.names(resALvU),"logP"=round(-log(resALvU$pvalue+1e-10,10),1)))
logs$logP=as.numeric(as.character(logs$logP))
sign=rep(1,nrow(logs))
sign[resALvU$log2FoldChange<0]=-1  ##change to correct model
table(sign)
# -1    1 
# 3683 3309 
logs$logP=logs$logP*sign
logs$gene=gsub("c", "isogroup", logs$gene)
write.table(logs,quote=F,row.names=F,file="GO_ALvU_logP.csv",sep=",")

####################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################
# heatmaps ------------
library(pheatmap)

rd <- read.csv("25jan2018_fgb_bm3_wald_RLDandPVALS.csv")
row.names(rd) <- rd$X
rd$X <- NULL
head(rd)

# check to make sure this is correct, but I should be able to replace "c" with "isogroup" to match gene file....

row.names(rd) <- gsub("c", "isogroup", row.names(rd))
head(rd)

names(rd)
exp <- rd[c(1:76)]
head(exp)

#     Load annotations
gg = read.delim("../tagseq/orb_fav_iso2gene.tab", header = F)
head(gg)

#     Make p-value cut-offs
sig.s = row.names(rd[rd$padj.s<0.000001 & !is.na(rd$padj.s),])
length(sig.s)

sig.b = row.names(rd[rd$padj.b<0.05 & !is.na(rd$padj.b),])
length(sig.b)

sig.alah = row.names(rd[rd$padj.alah<0.01 & !is.na(rd$padj.alah),])
length(sig.alah)

sig.alu = row.names(rd[rd$padj.alu<0.001 & !is.na(rd$padj.alu),])
length(sig.alu)

sig.ahu = row.names(rd[rd$padj.ahu<0.05 & !is.na(rd$padj.ahu),])
length(sig.ahu)

# sig expression
exp.s <- exp[row.names(exp) %in% sig.s,]
nrow(exp.s)

exp.b <- exp[row.names(exp) %in% sig.b,]
nrow(exp.b)

exp.alah <- exp[row.names(exp) %in% sig.alah,]
nrow(exp.alah)

exp.alu <- exp[row.names(exp) %in% sig.alu,]
nrow(exp.alu)

exp.ahu <- exp[row.names(exp) %in% sig.ahu,]
nrow(exp.ahu)

# species heatmap ---------
# naming the rows by gene names
gnames=c();expg=c()
for(i in row.names(exp.s)){
  s=subset(gg,V1==i)
  if (length(s[,1])>0){
    gnames=append(gnames,paste(s$V2[1],i,sep="."))
    expg=rbind(expg,exp.s[i,])
  } 
}
row.names(expg)=gnames
expl=expg
means=apply(expl,1,mean) # means of rows
explc=expl-means # subtracting them

# make colors 
heat.colors = colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")),bias=1.1)(100)

# cluster plot
quartz()
pheatmap(as.matrix(explc),color = heat.colors,cex=0.9,border_color=NA,clustering_distance_rows="correlation", cluster_cols=T)

# bank heatmap ---------
# naming the rows by gene names
gnames=c();expg=c()
for(i in row.names(exp.b)){
  s=subset(gg,V1==i)
  if (length(s[,1])>0){
    gnames=append(gnames,paste(s$V2[1],i,sep="."))
    expg=rbind(expg,exp.b[i,])
  } 
}
row.names(expg)=gnames
expl=expg
means=apply(expl,1,mean) # means of rows
explc=expl-means # subtracting them

# make colors 
heat.colors = colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")),bias=1.1)(100)

# cluster plot
quartz()
pheatmap(as.matrix(explc),color = heat.colors,cex=0.9,border_color=NA,clustering_distance_rows="correlation", cluster_cols=T)

# ALvAH heatmap ---------
# naming the rows by gene names
gnames=c();expg=c()
for(i in row.names(exp.alah)){
  s=subset(gg,V1==i)
  if (length(s[,1])>0){
    gnames=append(gnames,paste(s$V2[1],i,sep="."))
    expg=rbind(expg,exp.alah[i,])
  } 
}
row.names(expg)=gnames
expl=expg
means=apply(expl,1,mean) # means of rows
explc=expl-means # subtracting them

# make colors 
heat.colors = colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")),bias=1.9)(100)

# cluster plot
quartz()
pheatmap(as.matrix(explc),color = heat.colors,cex=0.9,border_color=NA,clustering_distance_rows="correlation", cluster_cols=T)

# ALvU heatmap ---------
# naming the rows by gene names
gnames=c();expg=c()
for(i in row.names(exp.alu)){
  s=subset(gg,V1==i)
  if (length(s[,1])>0){
    gnames=append(gnames,paste(s$V2[1],i,sep="."))
    expg=rbind(expg,exp.alu[i,])
  } 
}
row.names(expg)=gnames
expl=expg
means=apply(expl,1,mean) # means of rows
explc=expl-means # subtracting them

# make colors 
heat.colors = colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")),bias=1.9)(100)

# cluster plot
quartz()
pheatmap(as.matrix(explc),color = heat.colors,cex=0.9,border_color=NA,clustering_distance_rows="correlation", cluster_cols=T)

#Constrained analysis of proximities -------

library(vegan)
library(vegan3d)

load("results.Rdata")

# specify groups
head(coldata)
speciesv <- as.vector(coldata$species)
bankv <- as.vector(coldata$bank)
bankshape <- ifelse("east", "15", "16")
lesionv <- as.vector(coldata$lesion)
lesioncolv <- ifelse(lesionv=="U", "blue", ifelse(lesionv=="AH", "orange", "red"))

metav <- data.frame(cbind(speciesv, bankv, lesionv))
metav

# format expression data
rldd <- assay(rld)
head(rldd)

# constrained by a model:
cmd <- capscale(t(rldd)~speciesv+bankv+lesionv,metav,distance="manhattan")
cmd

# I'm not sure this is the correct way to calculate this...
eig <- eigenvals(cmd)[c(1:4)]
eig
perc.var <- eig/sum(eig)

cmd.scores <- scores(cmd)
cmd.scores

axes2plot <- c(1,2)  # try 2,3 too

plot(cmd,choices = axes2plot, display = "sites", xlab = paste("CAP1 ", round(perc.var[1],2)*100, "% variance explained",sep=""), ylab = paste("CAP2 ", round(perc.var[2],2)*100, "% variance explained",sep=""))
# points(cmd.scores$sites[bankv=="east",], pch=15)
# points(cmd.scores$sites[bankv=="west",], pch=16)
ordispider(cmd,choices = axes2plot, groups = metav$speciesv, col="grey50")
ordiellipse(cmd, choices = axes2plot, groups = metav$lesionv, draw="polygon", col = c("orange", "red", "blue"), label=T)

anova(cmd)
step(cmd)
# about ____% of variation is due to constraints (i.e. model factors)
adonis(t(rldd)~speciesv+bankv+lesionv,metav,distance="manhattan")

axes3plot <- c(1,2,3)
ordirgl(cmd, choices = axes3plot, pch = 16, col = c("orange", "red", "blue"))

# unconstrained (this is our usual PCoA):
cmd0 <- capscale(t(rldd)~1,distance="manhattan")

quartz()
axes2plot = c(1,2)  # try 2,3 too
plot(cmd0, choices = axes2plot, display="species") # choices - axes to display
ordispider(cmd0, choices = axes2plot, groups = metav$speciesv,col="grey50")
ordiellipse(cmd0, choices = axes2plot, groups = metav$lesionv, draw = "polygon", col = c("orange", "red", "blue"), label=T)


