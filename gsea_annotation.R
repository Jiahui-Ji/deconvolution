#######
#######
#######Cell-by-cell annotation by GSEA
#######
#######



#load required package
library("fgsea")
library("ggplot2")




#example data check
data(examplePathways) #class list with$
data(exampleRanks) #numeric






#make pathway list for gsea
#50 gene signatures for pn, mes, cl and mig; 39 gene signatures for met

#load gene signatures
pn
mes
cl

met
mig

met=met[1:39]

#make list for gene signatures for 5 GBM subtypes
marker=list()
marker$pn=pn
marker$mes=mes
marker$cl=cl
marker$met=met
marker$mig=mig






#input file
#single cell count data
gbm_normal #23368*3589

#marker list
load('gbm_subtype_marker.rds')




#test how many gene signarues in sc dataset
table(rownames(gbm_normal) %in% marker$pn)
table(rownames(gbm_normal) %in% marker$mes)
table(rownames(gbm_normal) %in% marker$cl)
table(rownames(gbm_normal) %in% marker$met)
table(rownames(gbm_normal) %in% marker$mig)




#loop to run GSEA for each cell
celltype=c()
for (i in 1:length(colnames(gbm_normal)))
{

	#gene list ordered by gene expression level
	gene_list=gbm_normal[,i]
	gene_list=sort(gene_list, decreasing=T)


	#run gsea 
	fgseaRes=fgsea(pathway=marker,stats=gene_list,minSize=15,
		maxSize=500,nperm=10000)


	tmp=as.data.frame(fgseaRes)
	tmp=tmp[tmp$ES > 0,]
	tmp=tmp[order(tmp$pval),]

	#select for p-value lower than 0.05 and NES over 0
	p_score=tmp$pval[1]
	nes_score=tmp$NES[1]

	if (p_score<0.05 & nes_score>0)
		{celltype=append(celltype,tmp$pathway[1])}
	else
		{celltype=append(celltype,'unclassified')}
}


#######done!













#######
#######
#######Show the annotated results in UMAP
#######
#######

#load package
library(colorout)
library(DropletUtils)
library(Seurat)
library(DoubletFinder)
library(dplyr)
library(ggplot2)
library(conos)

#seurat object
ControlSample1.data.B$gbm_subtype=celltype


pdf('gbm_subtype_v3.pdf')
DimPlot(ControlSample1.data.B, reduction = "umap",group.by="gbm_subtype",cols=c("DarkOrange","CornflowerBlue","LightCoral","SeaGreen","DarkOrchid","Gainsboro"))
dev.off()
#######done!


#0000FF00
pdf('gbm_subtype_v4.pdf')
DimPlot(ControlSample1.data.B, reduction = "umap",group.by="gbm_subtype",cols=c("#0000FF00","#0000FF00","LightCoral","SeaGreen","#0000FF00","#0000FF00"))
dev.off()





