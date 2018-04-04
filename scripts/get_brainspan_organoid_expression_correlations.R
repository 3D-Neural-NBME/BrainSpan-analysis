library(reshape2)
library(ggplot2)
library(plyr)
library(dplyr)
options(stringsAsFactors = F)

# This script modifies and merges the BrainSpan expression data with that generated for the manuscript and subsets down to genes with high expression and variance, as described in the manuscript. The resulting file is available as final.merged.brainspan.organoid.correlation.matrix.txt

# Data is available at : https://www.dropbox.com/sh/t0z9n8bk4moal1c/AAAVPWBookwqObd0ojXRl6TGa?dl=0
# BrainSpan data is available at : http://www.brainspan.org/static/download.html (RNA-Seq Gencode v10 summarized to genes)

sem = function(x){sd(x)/sqrt(length(x))}

# Load in reprocessed data 
artificial_tissue = read.delim("~/Dropbox (Personal)/3D_Neural_NBME:BrainSpan_analysis/gencode.v10.expression.matrix.txt")
rownames(artificial_tissue) = artificial_tissue$Row.names
artificial_tissue$Row.names = NULL
artificial_tissue.sample_names = read.delim("~/Dropbox (Personal)/3D_Neural_NBME:BrainSpan_analysis/samples_described_detail_sample_info.txt")

# Load in brainspan data 
brainspan = read.delim("genes_matrix_csv/expression_matrix.csv",header=F,row.names=1,sep=",")
colnames = read.delim("genes_matrix_csv/columns_metadata.csv",sep=",")
rownames = read.delim("genes_matrix_csv/rows_metadata.csv",sep=",")

# Format brainspan data 
names(brainspan) = paste(colnames$donor_id,gsub(" ",".",colnames$age),colnames$structure_acronym,sep="-")
rownames(brainspan) = rownames$ensembl_gene_id

# Subselect region and age data from Brainspan data 
tissues_of_interest = c("VFC","DFC","HIP","STR","AMY","MD","CBC","MFC","OFC","ITC","STC","A1C","V1C","M1C","IPC","S1C" )
timepoints_of_interest  =  c("8 pcw","9 pcw","12 pcw","13 pcw","16 pcw","17 pcw", "19 pcw","21 pcw","24 pcw","37 pcw","4 mos","1 yrs")         
columns_to_select = names(brainspan)[grep(paste(tissues_of_interest,collapse='|'),names(brainspan))]
columns_to_select_v2 = columns_to_select[grep(paste(gsub(" ",".",timepoints_of_interest),collapse='|'),columns_to_select)]
brainspan = subset(brainspan, select=columns_to_select_v2)

# Merge the datasets 
both = merge(artificial_tissue,brainspan,by="row.names")
rownames(both) = both$Row.names
both$Row.names = NULL

# Select high expressed and high variance genes from  brainspan
brainspan_genes_expression = data.frame(MeanExp = apply(brainspan,1,mean),VarExp = apply(brainspan,1,var),stringsAsFactors = F)
genes_to_select = rownames(subset(brainspan_genes_expression, MeanExp> 3 & VarExp > 1))
both = both[rownames(both) %in% genes_to_select,]

rm(brainspan,artificial_tissue,colnames,rownames)

# Get correlations
both_numeric = apply(both, 2, as.numeric)
both_numeric = both_numeric+1
both_numeric = log2(both_numeric)

both_cor = cor(both_numeric)
cor_df = data.frame(both_cor,stringsAsFactors = F)

# Select out correlations to only be between Halil's data and BrainSpan data 
cor_df$tissue = sapply(strsplit(rownames(cor_df),split="\\-"), function(x) x[3])

brain_regions = cor_df[cor_df$tissue %in% tissues_of_interest,]
brain_regions$tissue = NULL
brain_region_artificial = brain_regions[,-grep(paste(tissues_of_interest,collapse="|"),names(brain_regions))]

rm(both,both_cor,cor_df,both_numeric,brain_regions)

brain_region_artificial$brainspan = rownames(brain_region_artificial)

correlations = melt(brain_region_artificial)

correlations$variable = gsub("\\.","_",correlations$variable)
correlations$variable = gsub("^X","",correlations$variable)
correlations$variable = gsub("_R1_001$","",correlations$variable)
correlations = correlations[correlations$variable!="20C061516_S63",]
a = merge(correlations,artificial_tissue.sample_names,by.x="variable",by.y="X",all.x=TRUE)
a$tissue = sapply(strsplit(a$brainspan,split="-"),function(x) x[3])
a$timepoint = sapply(strsplit(a$brainspan,split="-"),function(x) x[2])

final = data.frame(a %>% group_by(condition,tissue,timepoint) %>% summarise(Mean=mean(value),Median=median(value),SEM=sem(value)) ,stringsAsFactors = F)

# Add in desired order for plots
orders = data.frame(tissue=tissues_of_interest,tissue_order = rev(letters[1:16]),stringsAsFactors = F)
final_annotated_tissue = merge(final,orders,by="tissue",all.x=T,all.y=T)
orders_timepoint = data.frame(timepoint = gsub(" ",".",timepoints_of_interest),timepoint_order=letters[1:12],stringsAsFactors = F)

final_annotated = merge(final_annotated_tissue,orders_timepoint,by="timepoint",all.x=T)
final_annotated = subset(final_annotated,timepoint!="21.yrs")
final_annotated = subset(final_annotated,timepoint!="11.yrs")

rm(artificial_tissue.sample_names,brain_region_artificial,final, final_annotated_tissue,orders,columns_to_select_v2,columns_to_select,genes_to_select,timepoints_of_interest,tissues_of_interest,brainspan_genes_expression,correlations)

write.table(final_annotated, "~/Dropbox (Personal)/3D_Neural_NBME:BrainSpan_analysis/final.merged.brainspan.organoid.correlation.matrix.txt", sep = "\t", quote=F, row.names=F, col.names=T)

write.table(orders_timepoint, "~/Dropbox (Personal)/3D_Neural_NBME:BrainSpan_analysis/orders_timepoint.txt", sep = "\t", quote=F, row.names=F, col.names=T)

