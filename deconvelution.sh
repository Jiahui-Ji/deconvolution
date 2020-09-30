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
#######Generate normlaized counts and associated cell type labels
#######
#######




#cluster.name -- cluster label
#gbm_normal -- count data
#ControlSample1.data.B -- seurat object

#generate cell label file
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



cell.name=colnames(ControlSample1.data.B)
cell.label.matrix=matrix(nrow=length(cell.name),ncol=2)
colnames(cell.label.matrix)=c('Cell name','Celltype') #Attention! have to name the second column as 'Celltype'
cell.label.matrix[,1]=cell.name
cell.label.matrix[,2]=cluster.name

write.table(cell.label.matrix,file='gbmdata_celltypes.txt',quote=F,sep='\t') #cell label information






'''
mat=read.table(file='GBM_normalized_gene_counts.csv')

write.table(t(mat),'adata_norm_counts_all.txt',quote=F,sep='\t')
'''




#load python
import os 
import numpy as np
import pandas as pd
import time as tm

import scanpy as sc



#GBM_raw_gene_counts.csv
#GBM_normalized_gene_counts.csv

cell_label=pd.read_csv('gbm_celltype_anotation.csv',sep=',')
gbm_sc=pd.read_csv('GBM_raw_gene_counts.csv',sep='\s')

gbm_sc=np.transpose(gbm_sc)

#normalization
adata=gbm_sc
adata_normal=sc.pp.normalize_per_cell(adata,copy=True)

adata_normal.index=cell_label['cell name']



#save file with the required file name
adata_normal.to_csv('gbmdata_norm_counts_all.txt', sep='\t')






#######done!



























#######
#######
####### start runing scaden, bulk simulation
#######
#######


#bulk simulation
python bulk_simulation.py --cells 100 --samples 1000 --data /data/jj1419/gbm/GBM_data_and_metadata/preprocessing/


#create h5ad file, there is one bug, change as_matrix() to values in the original py file
python create_h5ad_file.py --data /data/jj1419/gbm/GBM_data_and_metadata/preprocessing/ --out gbm.h5ad



########done!



























#######
#######
####### train and test
#######
#######






#preprocess your traning data, training data will be log2-transformmed
#and scaled to the range [0,1]
scaden process gbm.h5ad /Users/jijiahui/desktop/preprocessing/pn_cell.txt


#training, train a model for 20,000 steps
scaden train processed.h5ad 


#prediction
scaden predict pn_cell.txt
scaden predict mes_cell.txt
scaden predict met_cell.txt
scaden predict cl_cell.txt
scaden predict mig_cell.txt



#Add the following codes into the functions.py
'''
new=[0 for x in range(0,len(available_genes))]
for i in range(0,len(available_genes)):
	x='"'
	y=available_genes[i]
	z='"'
	q=x+y+z
	new[i]=q


available_genes=new

new_sig_genes = list(set(available_genes).intersection(sig_genes_complete))
'''


#Add the following codes into the scaden.py
'''
data_index = list(data.index)
new=[0 for x in range(0,len(data_index))]
for i in range(0,len(data_index)):
	x='"'
	y=data_index[i]
	z='"'
	q=x+y+z
	new[i]=q


data.index=new


data = data.loc[sig_genes]
'''



#######done!










