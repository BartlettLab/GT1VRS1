##RUVseq-assited DE analysis
library(RUVSeq)
library(DESeq2)

##Tools from EDASeq are included in RUVSeq to explore the gene expression data
cts <- as.matrix(read.table("htseq_counts.allLibs.txt",header=T,row.names="gene"))
colnames(cts) <- gsub(".trimAligned.sortedByCoord.out.bam","",colnames(cts))
#remove special rows from HTseq files
cts <- cts[-grep("__",rownames(cts)),]
colSums(cts)

#set up colData object
coldata <- read.table("RNA-seq_metadata.txt",sep="\t",header=T,row.names=1)
coldata$gt1 <- factor(coldata$gt1)
coldata$vrs1 <- factor(coldata$vrs1)
coldata$mut <- factor(coldata$mut)
coldata$mut <- factor(coldata$mut, levels=c("wt","vrs1-het","gt1-het","gt1-het_vrs1-het","vrs1","gt1-het_vrs1","gt1","gt1_vrs1-het","gt1_vrs1"))
coldata$ear.size.f <- "B"
coldata$ear.size.f[coldata$ear.size > 1.4] <- "C"
coldata$ear.size.f[coldata$ear.size < 0.7] <- "A"
#remove un-genotyped individuals
noGeno <- is.na(coldata$mut)
cts <- cts[,!noGeno]
coldata <- coldata[!noGeno,]
coldata$index <- paste0(coldata$mut,coldata$ear.size.f)

#filter data for gene expression minimum
dim(cts)
#[1] 39756    53
filter<-apply(cts,1,function(x)length(x[x>5])>=2)
ctsFilter<-cts[filter,]
dim(ctsFilter)
#[1] 29262    53

#make EDAseq SeqExpressionSet
set<-newSeqExpressionSet(as.matrix(ctsFilter), phenoData=coldata)

library(RColorBrewer) 
colors<-brewer.pal(9,"Set1") 
plotRLE(set,outline=FALSE,ylim=c(-4,4),col=colors[set$mut]) 
plotPCA(set,col=colors[set$mut],cex=1.2)

#RUVSeq suggests an upper quartile normalization
set <-betweenLaneNormalization(set,which="upper")
plotRLE(set,outline=FALSE,ylim=c(-4,4),col=colors[set$mut]) 
plotPCA(set,col=colors[set$mut],cex=1.2)

#Convert data to DESeq2 object for preliminary DE to get "negative" control genes
dds<-DESeqDataSetFromMatrix(countData=counts(set), 
				    colData=pData(set),
				    design= ~ index)
dds <- DESeq(dds)

#get the least differentially expressed genes
summary(results(dds,contrast=c("index","gt1_vrs1B","wtB")),alpha=0.05)
res <- results(dds,contrast=c("index","gt1_vrs1B","wtB"))
resOrdered <- res[order(res$pvalue,decreasing=TRUE),]
empirical <- rownames(resOrdered)[1:5000] #the top 5000 least DE genes in the gt1 vrs1 vs wt comparison

setRUV <- RUVg(set,empirical,k=1)
head(pData(setRUV))
#      ear.size ear.size.f vrs1 gt1      mut     index         W_1
#jg107      2.4          C    2   0     vrs1     vrs1C -0.55938897
#jg121      1.0          B    0   0       wt       wtB  0.04337470
#jg145      1.0          B    0   1  gt1-het  gt1-hetB  0.04997112
#jg153      0.5          A    2   2 gt1_vrs1 gt1_vrs1A  0.04801004
#jg163      0.6          A    0   1  gt1-het  gt1-hetA  0.02518301
#jg164      1.2          B    0   0       wt       wtB  0.04700921
plotRLE(setRUV,outline=FALSE,ylim=c(-4,4),col=colors[setRUV$mut]) 
plotPCA(setRUV,col=colors[factor(setRUV$ear.size.f)],cex=1.2)
#Note: the PCA method in DESeq2 cannot easily represent the data from the RUVseq
#or from the EDAseq objects, so it is not worth trying to squeeze that data in

#Redo DESeq2 with weights added in the formula
dds<-DESeqDataSetFromMatrix(countData=counts(setRUV), 
				    colData=pData(setRUV),
				    design= ~ W_1 + index)
levels(dds$mut) <- gsub("-",".",levels(dds$mut))
levels(dds$index) <- gsub("-",".",levels(dds$index))

dds <- DESeq(dds)

#get all pairwise comparisons
res <- results(dds,contrast=c("index","gt1_vrs1B","wtB"))
resOrdered <- res[order(res$padj),]
summary(res,alpha=0.05)
write.csv(as.data.frame(resOrdered),file = "gt1vrs1B_vs_wtB.csv")

res.gt1 <- results(dds,contrast=c("index","gt1B","wtB"))
resOrdered <- res[order(res$padj),]
summary(res,alpha=0.05)
write.csv(as.data.frame(resOrdered),file = "gt1B_vs_wtB.csv")

res.vrs1 <- results(dds,contrast=c("index","vrs1B","wtB"))
resOrdered <- res[order(res$padj),]
summary(res,alpha=0.05)
write.csv(as.data.frame(resOrdered),file = "vrs1B_vs_wtB.csv")

res.gt1vrs1 <- results(dds,contrast=c("index","gt1_vrs1B","gt1B"))
resOrdered <- res[order(res$padj),]
summary(res,alpha=0.05)
write.csv(as.data.frame(resOrdered),file = "gt1vrs1B_vs_gt1B.csv")

res <- results(dds,contrast=c("index","gt1_vrs1B","vrs1B"))
resOrdered <- res[order(res$padj),]
summary(res,alpha=0.05)
write.csv(as.data.frame(resOrdered),file = "gt1vrs1B_vs_vrs1B.csv")

res <- results(dds,contrast=c("index","gt1B","vrs1B"))
resOrdered <- res[order(res$padj),]
summary(res,alpha=0.05)
write.csv(as.data.frame(resOrdered),file = "gt1B_vs_vrs1B.csv")

res <- results(dds,contrast=c("index","gt1B","gt1.hetB"))
resOrdered <- res[order(res$padj),]
summary(res,alpha=0.05)
write.csv(as.data.frame(resOrdered),file = "gt1B_vs_gt1hetB.csv")

res <- results(dds,contrast=c("index","gt1.hetB","wtB"))
resOrdered <- res[order(res$padj),]
summary(res,alpha=0.05)
write.csv(as.data.frame(resOrdered),file = "gt1hetB_vs_wtB.csv")

res <- results(dds,contrast=c("index","gt1.het_vrs1B","vrs1B"))
resOrdered <- res[order(res$padj),]
summary(res,alpha=0.05)
write.csv(as.data.frame(resOrdered),file = "gt1het_vrs1B_vs_vrs1B.csv")

res <- results(dds,contrast=c("index","gt1.het_vrs1B","wtB"))
resOrdered <- res[order(res$padj),]
summary(res,alpha=0.05)
write.csv(as.data.frame(resOrdered),file = "gt1het_vrs1B_vs_wt1B.csv")

#add in gene data from genome
geneData <- read.delim("Zm00001eb.1.fulldata.txt",header=T,sep="\t")#make sure this file is complete
sigDiff <- resOrdered$padj < 0.05 
sigDiff[is.na(sigDiff)] = FALSE
res.sigDiff <- resOrdered[sigDiff,]

sigDiff.geneData <- geneData[geneData$gene_model %in% rownames(res.sigDiff),]
write.table(sigDiff.geneData,file="gt1B_vs_gt1hetB.sigDiff.txt",quote=FALSE,sep="\t",row.names=FALSE)
#do this for all the comparisons in results section

#make list of all DE genes
sigGeneFiles <- list.files(path="./Results/RUVSeq_results",pattern="wtB.sigDiff")
sigGeneFiles <- paste0("./Results/RUVSeq_results/",sigGeneFiles)
sigGeneList <- list()
sigGenes <- vector()
for(i in 1:length(sigGeneFiles)){
	temp.df <- read.table(file = sigGeneFiles[i],sep="\t",header = TRUE)
	sigGeneList[[i]] <- temp.df$gene_model
	sigGenes <- append(sigGenes,temp.df$gene_model)
}
uniqSigGenes <- unique(sort(sigGenes))

#all shared DE genes
sigDiff.all.geneData <- geneData[geneData$gene_model %in% uniqSigGenes[table(sigGenes) > 4],]
dim(sigDiff.all.geneData)#52 rows
write.table(sigDiff.all.geneData,file="DEG.all.sigDiff.txt",quote=FALSE,sep="\t",row.names=FALSE)

#make Venn diagram of gt1, vrl1, and gt1 vrl1 in 
library(ggvenn)
ggvenn(sigGeneList,columns=c(1,4,5))

shared.sigDE <- intersect(intersect(sigGeneList[[4]],sigGeneList[[1]]),sigGeneList[[5]])
gt1vrs1_uniq.sigDE <- setdiff(setdiff(sigGeneList[[4]],sigGeneList[[1]]),sigGeneList[[5]])
gt1vrs1all <- union(shared.sigDE,gt1vrs1_uniq.sigDE)

sigDiff.gt1vrs1.geneData <- geneData[geneData$gene_model %in% gt1vrs1_uniq.sigDE,]
dim(sigDiff.gt1vrs1.geneData)#52 rows
write.table(sigDiff.gt1vrs1.geneData,file="DEG.gt1vrs1_uniq.sigDiff.txt",quote=FALSE,sep="\t",row.names=FALSE)

res.gt1vrs1.uniq <- res.sigDiff[rownames(res.sigDiff) %in% gt1vrs1_uniq.sigDE,]
sigDE.gt1vrs1.up <- res.gt1vrs1.uniq[res.gt1vrs1.uniq$log2FoldChange > 0,]
sigDE.gt1vrs1.down <- res.gt1vrs1.uniq[res.gt1vrs1.uniq$log2FoldChange < 0,]

#GO term analysis - revisit this after building shared DE list

#GO terms are in the geneData data frame (listed as obo_terms)
library(topGO)
head(geneData)

obo_terms <- geneData[,c("gene_model","obo_terms")]
write.table(obo_terms,file="gene_GO_terms.Zm00001eb.1.map",quote=FALSE,sep="\t",row.names=FALSE)
#edit file to below format for reading in with readMappings()
#gene_ID<TAB>GO_ID1, GO_ID2, GO_ID3, ....

#read file back in
geneID2GO <- readMappings(file = "gene_GO_terms.Zm00001eb.1.map")
geneNames <- names(geneID2GO)
geneList <- factor(as.integer(geneNames %in% rownames(sigDE.gt1vrs1.up)))
names(geneList) <- geneNames
GOData <- new("topGOdata", ontology = "BP", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
resultFisher <- runTest(GOData, algorithm = "classic", statistic = "fisher")
allRes <- GenTable(GOData, classicFisher = resultFisher, ranksOf = "classicFisher", topNodes = length(resultFisher@score))
allRes$qval.bh<-p.adjust(as.numeric(allRes[,"classicFisher"]),method="BH")
write.table(allRes, file = "sigDE.gt1vrs1.BP_GO.up.txt")

#make heatmap of logFC 
library(pheatmap)
#initialize with gt1vrs1 DE genes
logFC.df <- data.frame(Gene = rownames(res.gt1vrs1.uniq), gt1vrs1.log2FC = res.gt1vrs1.uniq$log2FoldChange)
#add columns to logFC.df
for(i in 1:length(logFC.df$Gene)){
	logFC.df$gt1.log2FC[i] <- res$log2FoldChange[rownames(res) %in% logFC.df$Gene[i]]
}

for(i in 1:length(logFC.df$Gene)){
	logFC.df$vrs1.log2FC[i] <- res$log2FoldChange[rownames(res) %in% logFC.df$Gene[i]]
}

for(i in 1:length(logFC.df$Gene)){
	logFC.df$gt1hetvrs1.log2FC[i] <- res$log2FoldChange[rownames(res) %in% logFC.df$Gene[i]]
}

rownames(logFC.df) <- logFC.df$Gene
logFC.df <- logFC.df[,-1]
pheatmap(logFC.df)
for(i in 1:length(rownames(logFC.df))){
	logFC.df$Locus[i] <- geneData$locus_symbol[geneData$gene_model %in% rownames(logFC.df)[i]]
}

logFC.locus.df <- logFC.df
replace.list <- grep(pat = "[A-Za-z]",logFC.df$Locus,invert=TRUE)
for(i in replace.list){
	logFC.locus.df$Locus[i] <- rownames(logFC.locus.df)[i]
}
rownames(logFC.locus.df) <- logFC.locus.df$Locus
pheatmap(logFC.locus.df[,c(2,3,4,1)],cluster_cols=F)

pdf("heatmap_of_gt1vrs1_DEGs.log2FC.pdf",height = 11, width = 8.5)
pheatmap(logFC.locus.df[,c(2,3,4,1)], cluster_cols = FALSE)
dev.off()

#alternate graph with low expression differences excluded
table(abs(logFC.locus.df$gt1vrs1.log2FC) >= 1)
pheatmap(logFC.locus.df[abs(logFC.locus.df$gt1vrs1.log2FC) >= 1.0,c(2,3,4,1)], cluster_cols = FALSE)
