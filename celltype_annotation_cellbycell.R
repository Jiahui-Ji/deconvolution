######load packages
library(colorout)
library(DropletUtils)
library(Seurat)
library(DoubletFinder)
library(dplyr)
library(ggplot2)
library(conos)






#######
#######
####### cell-by-cell celltype annotatioin by Fisher's exact test 
#######
#######






#background gene, genes that express in over 10% cells

#in GBM_data_and_metadata
#GBM_normalized_gene_counts.csv


gbm_normal=read.table('GBM_normalized_gene_counts.csv',header=T)
gbm_normal=as.matrix(gbm_normal)  #dim, 23368x3589 


all.genes=rownames(gbm_normal)



#background gene, genes that express in at least 10% of cells
background.gene=c()
for (i in 1:length(rownames(gbm_normal)))
{

	n=0
	for (j in 1:length(colnames(gbm_normal)))
	{

		if (gbm_normal[i,j]>0)
			{n=n+1}

	}


	propotion=n/length(colnames(gbm_normal))
	if (propotion >= 0.1)
	{background.gene=append(background.gene,rownames(gbm_normal)[i])}	


}



#the number of background gene is 8401











load('mes.rds') #mes, 50 mes marker gene
load('pn.rds') #pn, 50 pn marker gene 
load('cl.rds') #cl, 50 cl marker gene
load('met.rds')
load('mig.rds')

met.select=as.factor(met.select)
mig.select=as.factor(mig.select)

#do Fisher's exact test
All.genes=background.gene
Clinical.marker=mig.select  #change the marker gene list here

TMP=matrix(ncol=length(colnames(gbm_normal)),nrow=1)
for (i in 1:length(colnames(gbm_normal)))
{
 
	number=c()
		for (j in 1:length(rownames(gbm_normal)))
		{
			if (gbm_normal[j,colnames(gbm_normal)[i]]>0)
			
			{number=append(number,j)}

		}

	TMP.Genes=rownames(gbm_normal)[number] 
	#print(length(TMP.Genes))
	#consider all genes expressed expressed over 0 in one cell

	Overlap.Genes=intersect(Clinical.marker,All.genes)

	TMP.MAT=matrix(ncol=2,nrow=2)
    TMP.MAT[1,1]=length(intersect(TMP.Genes,Overlap.Genes))
    TMP.MAT[1,2]=length(setdiff(Overlap.Genes,TMP.Genes))
    TMP.MAT[2,1]=length(setdiff(TMP.Genes,Overlap.Genes))
    TMP.MAT[2,2]=length(All.genes)-TMP.MAT[1,1]-TMP.MAT[1,2]-TMP.MAT[2,1]
    
    if(TMP.MAT[2,2]>0)
    	{TMP[1,i]=fisher.test(TMP.MAT,alternative="greater")$p.value}
	else
		{TMP[1,i]=NA}

}


mes.annotation=TMP
pn.annotation=TMP
cl.annotation=TMP
met.annotation=TMP
mig.annotation=TMP

mes.cells=colnames(gbm_normal)[which(mes.annotation<0.05)]
pn.cells=colnames(gbm_normal)[which(pn.annotation<0.05)]
cl.cells=colnames(gbm_normal)[which(cl.annotation<0.05)]
met.cells=colnames(gbm_normal)[which(met.annotation<0.05)]
mig.cells=colnames(gbm_normal)[which(mig.annotation<0.05)]




save(met.cells,file='met_cells.rds')
save(mig.cells,file='mig_cells.rds')
#######done!



























#######
#######
#######MES, CL and PN type UMAP 
#######
#######



#annotate cell types
cell.name=colnames(gbm_normal)

subtype=rep("unknown",length(cell.name))
for (i in 1:length(cell.name))
{
	if (cell.name[i] %in% mes.cells)
	{subtype[i]='MES type'}
	if (cell.name[i] %in% pn.cells)
	{subtype[i]='PN type'}
	if (cell.name[i] %in% cl.cells)
	{subtype[i]='CL type'}
	
}




#annotate cell types after adding metabolic and migratory subtypes
subtype=rep("unknown",length(cell.name))
for (i in 1:length(cell.name))
{
	if (cell.name[i] %in% mes.cells)
	{subtype[i]='MES type'}
	if (cell.name[i] %in% pn.cells)
	{subtype[i]='PN type'}
	if (cell.name[i] %in% cl.cells)
	{subtype[i]='CL type'}
	if (cell.name[i] %in% met.cells)
	{subtype[i]='MET type'}
	if (cell.name[i] %in% mig.cells)
	{subtype[i]='MIG type'}
	
}




#create seurat object and plot UMAP
ControlSample1.data= CreateSeuratObject(counts =gbm_normal, min.cells = 5,project = "pd")
ControlSample1.data

ControlSample1.data.B= NormalizeData(ControlSample1.data, normalization.method = "LogNormalize", scale.factor = 10000)


ControlSample1.data.B = FindVariableFeatures(ControlSample1.data.B, selection.method = "mean.var.plot", nfeatures = 2000)
All.genes = rownames(ControlSample1.data.B)

ControlSample1.data.B = ScaleData(ControlSample1.data.B, features = All.genes)
ControlSample1.data.B = RunPCA(ControlSample1.data.B,npcs = 30, verbose = FALSE)
ControlSample1.data.B = FindNeighbors(ControlSample1.data.B, dims = 1:10)
ControlSample1.data.B= FindClusters(ControlSample1.data.B, resolution = 1)
ControlSample1.data.B = RunUMAP(ControlSample1.data.B, dims = 1:10)




ControlSample1.data.B$celltype=subtype



#cell-type annotaion in paper
tsne=read.table(file='GBM_TSNE.csv',header=T)
meta=read.table(file='GBM_metadata.csv',header=T)


celltype.color=meta$Selection
location=meta$Location
clusters=meta$Cluster_2d



color=c('#C5B0D5','#D62728','#FF9896',
		'#1F77B4','#FFBB78','#8C564B',
		'#2CA02C','#AEC7E8','#9467BD',
		'#FF7F0E','#98DF8A','#C49C94'
	    )

meta$Splice_sites_Annotate



#change cluster label to name in paper
cluster.name=c()

for (i in clusters)
{
	if (i == 1)
		{cluster.name=append(cluster.name,"Neoplastic cell 1")}

	if (i == 2)
		{cluster.name=append(cluster.name,"Oligodendrocytes")}

	if (i == 3)
		{cluster.name=append(cluster.name,"Vascular cells 1")}

	if (i == 4)
		{cluster.name=append(cluster.name,"Neoplastic cell 2")}

	if (i == 5)
		{cluster.name=append(cluster.name,"Neurons")}

	if (i == 6)
		{cluster.name=append(cluster.name,"Vascular cells 2")}

	if (i == 7)
		{cluster.name=append(cluster.name,"Myeloid cells 1")}

	if (i == 8)
		{cluster.name=append(cluster.name,"Myeloid cells 2")}

	if (i == 9)
		{cluster.name=append(cluster.name,"OPCs")}

	if (i == 10)
		{cluster.name=append(cluster.name,"Astrocytes")}

	if (i == 11)
		{cluster.name=append(cluster.name,"Neoplastic cell 3")}

	if (i == 12)
		{cluster.name=append(cluster.name,"Vascular cells 3")}


}





ControlSample1.data.B$antibody=celltype.color
ControlSample1.data.B$location=location
ControlSample1.data.B$cluster=cluster.name













pdf('gbm_subtype.pdf')
DimPlot(ControlSample1.data.B, reduction = "umap",group.by="celltype",cols=c("DarkOrange","CornflowerBlue","LightCoral","Gainsboro"))
dev.off()

pdf('gbm_paper_annotation.pdf')
DimPlot(ControlSample1.data.B, reduction = "umap",group.by="antibody")
dev.off()

pdf('gbm_location.pdf')
DimPlot(ControlSample1.data.B, reduction = "umap",group.by="location")
dev.off()

pdf('gbm_paper_cluster.pdf')
DimPlot(ControlSample1.data.B, reduction = "umap",group.by="cluster",cols=color)
dev.off()

pdf('gbm_subtype_v2.pdf')
DimPlot(ControlSample1.data.B, reduction = "umap",group.by="celltype",cols=c("DarkOrange","CornflowerBlue","LightCoral","SeaGreen","DarkOrchid","Gainsboro"))
dev.off()


























#######
#######
####### Get the gene signatures for MIG and MEG
#######
#######



write.table(count.data,file="gbm_signature.csv",sep=',',quote=F)



gbm_count=read.table(file='gbm_signature.txt',sep='\t',header=T)




#remove the NA value
data=na.omit(gbm_count)  #25218*153




row=data[,1]
col=colnames(data)
col=col[1:152]


gbm_data=data[1:25218,2:153]
colnames(gbm_data)=col

gbm_data=as.matrix(gbm_data)
rownames(gbm_data)=row  #gbm_data, 25218*152
gbm_data=unique(gbm_data) #24348*152 final gbm count data



#load cluster1.name,cluster2.name,cluster3.name,cluster4.name,cluster5.name







#metabolic subtype amd migrotory label
col.new=gsub('[.]','-',col)
colnames(gbm_data)=col.new

met.label=rep(0,152)
for (i in 1:length(col.new))
{
	if (col.new[i] %in% cluster3.name)
		{met.label[i]=1}
	
}

mig.label=rep(0,152)
for (i in 1:length(col.new))
{
	if (col.new[i] %in% cluster5.name)
		{mig.label[i]=1}
	
}


load('cluster1_name.rds')
load('cluster2_name.rds')
load('cluster3_name.rds')
load('cluster4_name.rds')
load('cluster5_name.rds')


pn_cell=gbm_data[,cluster1.name]
mes_cell=gbm_data[,cluster2.name]
met_cell=gbm_data[,cluster3.name]
cl_cell=gbm_data[,cluster4.name]
mig_cell=gbm_data[,cluster5.name]



write.table(pn_cell,file="pn_cell.txt",quote=F,sep='\t')
write.table(mes_cell,file="mes_cell.txt",quote=F,sep='\t')
write.table(met_cell,file="met_cell.txt",quote=F,sep='\t')
write.table(cl_cell,file="cl_cell.txt",quote=F,sep='\t')
write.table(mig_cell,file="mig_cell.txt",quote=F,sep='\t')




#prepare for anova test
final.result=matrix(ncol=2,nrow=length(colnames(data.anova))-2)
rownames(final.result)=colnames(data.anova)[1:24348]
colnames(final.result)=c('MET','MIG')

data.anova=t(gbm_data)

data.anova=cbind(data.anova,met.label,mig.label)
data.anova=as.data.frame(data.anova)


#loop for anova test
for (i in 1:24348)
{
	q=aov(data.anova[,i]~data.anova$met.label,data.anova)
	final.result[i,1]=summary(q)[[1]][["Pr(>F)"]][[1]]

	w=aov(data.anova[,i]~data.anova$mig.label,data.anova)
	final.result[i,2]=summary(w)[[1]][["Pr(>F)"]][[1]]
}

save(final.result,file="final_meg_pvalue.rds")


#selected MET genes
met.genes=rownames(final.result[which(final.result[,1]<1e-03),]) #2818

#selected MIG genes
mig.genes=rownames(final.result[which(final.result[,2]<1e-03),]) #1875









#MET table
met.table=gbm_data[met.genes,]

met.samples=which(met.label==1)
other.not.met=which(met.label==0)


met.expression=matrix(nrow=length(rownames(met.table)),ncol=2)
colnames(met.expression)=c('MET expreesion','rest expression')
rownames(met.expression)=rownames(met.table)


for (i in 1:length(rownames(met.table)))
{
	met.expression[i,1]=mean(met.table[i,met.samples])
	met.expression[i,2]=mean(met.table[i,other.not.met])
}


met.select=rownames(met.expression[which(met.expression[,1]>met.expression[,2]),])
#1493 genes




#MIG table
mig.table=gbm_data[mig.genes,]

mig.samples=which(mig.label==1)
other.not.mig=which(mig.label==0)


mig.expression=matrix(nrow=length(rownames(mig.table)),ncol=2)
colnames(mig.expression)=c('MIG expreesion','rest expression')
rownames(mig.expression)=rownames(mig.table)


for (i in 1:length(rownames(mig.table)))
{
	mig.expression[i,1]=mean(mig.table[i,mig.samples])
	mig.expression[i,2]=mean(mig.table[i,other.not.mig])
}


mig.select=rownames(mig.expression[which(mig.expression[,1]>mig.expression[,2]),])
#1870 genes













#gene expresion data and do t-test to compare mean gene expression
all.express=matrix(nrow=length(rownames(gbm_data)),ncol=6)
colnames(all.express)=c('MET','rest not MET','MIG','rest not MIG','MET p-value','MIG p-value')
rownames(all.express)=rownames(gbm_data)

for (i in 1:length(rownames(gbm_data)))
{
	all.express[i,1]=mean(gbm_data[i,met.samples])
	all.express[i,2]=mean(gbm_data[i,other.not.met])
	all.express[i,3]=mean(gbm_data[i,mig.samples])
	all.express[i,4]=mean(gbm_data[i,other.not.mig])
}


all.express=all.express[which(all.express[,1]!=all.express[,2]),]
all.express=all.express[which(all.express[,3]!=all.express[,4]),]

for (i in 1:length(rownames(all.express)))
{
	all.express[i,5]=t.test(all.express[i,1:2])$p.value
	all.express[i,6]=t.test(all.express[i,3:4])$p.value

}


met.genes=rownames(all.express[which(all.express[,5]<1e-03 & all.express[,1]>all.express[,2]),])  #62 
mig.genes=rownames(all.express[which(all.express[,6]<1e-03 & all.express[,3]>all.express[,4]),])  #82



#the 50 gene signatures for pn, mes and cl 
#no overlap between them


met.matrix=all.express[met.genes,]
met.matrix=met.matrix[order(met.matrix[,5]),]

met.select=rownames(met.matrix)[1:50]
save(met.select,file='met.rds')



mig.matrix=all.express[mig.genes,]
mig.matrix=mig.matrix[order(mig.matrix[,6]),]

mig.select=rownames(mig.matrix)[1:50]
save(mig.select,file='mig.rds')

#######get the 50 gene signatures for the novel MET and MIG subtype
#######done!








