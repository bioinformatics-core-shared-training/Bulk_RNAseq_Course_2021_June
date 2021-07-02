library(DESeq2)
library(tidyverse)


# load data

txi <- readRDS( 'RObjects/txi.rds' )
class(txi)
View(txi$counts)
dim(txi$counts)
View(txi$length)

# Sample info
sampleinfo <- read_tsv( 'data/samplesheet_corrected.tsv')

# check sample order

all(sampleinfo$SampleName == colnames(txi$counts))

#  create DESeqDataSet object
# y = 2x
# y ~ 2x
simple.model <- as.formula(~Status)
model.matrix(simple.model, data = sampleinfo)

as.factor(sampleinfo$Status)

sampleinfo <- mutate( sampleinfo, Status = fct_relevel( Status, 'Uninfected', 'Infected'))

sampleinfo$Status

model.matrix( simple.model, data=sampleinfo)

# Build DDS object
ddsObj.raw <- DESeqDataSetFromTximport( txi=txi, colData = sampleinfo, design = simple.model)

# filter low count genes
keep <- rowSums( counts(ddsObj.raw)) > 5

ddsObj.filt <- ddsObj.raw[ keep, ]

dim(counts(ddsObj.raw))

dim(counts(ddsObj.filt))

# Estimate size factors
ddsObj <- estimateSizeFactors( ddsObj.filt)

normalizationFactors(ddsObj.filt)
dim(normalizationFactors(ddsObj))

apply(normalizationFactors(ddsObj), 2, median)

# MA plot

logcounts <- log2( counts(ddsObj, normalized =FALSE ) + 1)
head(logcounts)

limma::plotMA( logcounts, array=5, ylim=c(-5,5))
abline( h = 0, col='red')

# MA plot- Normalized counts

logcounts <- log2( counts(ddsObj, normalized =TRUE ) + 1)

limma::plotMA( logcounts, array=5, ylim=c(-5,5))
abline( h = 0, col='blue')

# estimate dispersions
ddsObj <- estimateDispersions(ddsObj)
plotDispEsts(ddsObj)

# fitting model and testing
ddsObj <- nbinomWaldTest( ddsObj)

# DESeq = estimateSizeFactors + estimateDispersions + nbinomWaldTest

ddsObj <- DESeq( ddsObj.filt)

# Results tables
results.simple <- results(ddsObj, alpha = 0.05)
results.simple

sum(results.simple$padj < 0.05, na.rm=TRUE)

View( as.data.frame(results.simple))

# Upregulated
sum(results.simple$padj < 0.05 & results.simple$log2FoldChange > 0, na.rm=TRUE )

# Down regulated
sum(results.simple$padj < 0.05 & results.simple$log2FoldChange < 0, na.rm=TRUE )

###########################################################################################
# Exercise 1
# var = SS / df
# df = n  
additive.model <- as.formula( ~ TimePoint + Status   )
model.matrix(model.additive, data= sampleinfo)

ddsObj.raw <- DESeqDataSetFromTximport( txi=txi, colData = sampleinfo, design = additive.model)

keep <- rowSums(counts(ddsObj.raw)) > 5
ddsObj.filt <- ddsObj.raw[keep,]

ddsObj <- DESeq( ddsObj.filt)

# extract results
results.additive <- results(ddsObj, alpha=0.05)
results.additive

sum(results.additive$padj < 0.05, na.rm=TRUE)


model.matrix(model.additive, data= sampleinfo)

resultsNames(ddsObj)

results.InfectedvUninfected <- results.additive
rm(results.additive)


# Exercise 2
resultsNames(ddsObj)

results.d33vd11 <- results( ddsObj, name='TimePoint_d33_vs_d11', alpha = 0.05)
results.d33vd11
sum(results.d33vd11$padj < 0.05, na.rm=TRUE)


# PCA
vstcounts <- vst( ddsObj.raw, blind=TRUE)
dim(vstcounts)
plotPCA( vstcounts, intgroup = c( 'Status', 'TimePoint'))


# Comparing two design models
# Additive vs Simple model
# ~ TimePoint + Status
# ~ Status

design(ddsObj)
ddsObj.LRT <- DESeq( ddsObj, test='LRT', reduced = simple.model)

results.Additive_v_Simple <- results(ddsObj.LRT)
results.Additive_v_Simple
sum(results.Additive_v_Simple$padj < 0.05, na.rm=TRUE)

# Exercise 3
interaction.model <- as.formula( ~ TimePoint + Status + TimePoint:Status)


interaction.model <- as.formula( ~ TimePoint * Status)
model.matrix(interaction.model, data=sampleinfo)

ddsObj.raw <- DESeqDataSetFromTximport(txi = txi,
                                       colData = sampleinfo,
                                       design = interaction.model)

keep <- rowSums(counts(ddsObj.raw)) > 5
ddsObj.filt <- ddsObj.raw[keep,]

ddsObj.interaction <- DESeq( ddsObj.filt )

ddsObj.LRT <- DESeq( ddsObj.interaction, test='LRT', reduced = additive.model)

results.Interaction_v_Additive <- results(ddsObj.LRT)

sum(results.Interaction_v_Additive$padj < 0.05, na.rm=TRUE)

# Extract results from interaction model
resultsNames(ddsObj.interaction)

results.interaction.11 <- results(ddsObj.interaction, alpha = 0.05,
                                  name='Status_Infected_vs_Uninfected')

results.interaction33 <- results(ddsObj.interaction,
                                 alpha = 0.05,
                                 contrast = list( c( 'Status_Infected_vs_Uninfected', 'TimePointd33.StatusInfected') )
                                 
)

sum(results.interaction.11$padj < 0.05, na.rm=T)

sum(results.interaction33$padj < 0.05, na.rm=T)


# save the results
saveRDS(ddsObj.interaction, "results/DESeqDataSet.interaction.rds")
saveRDS(results.interaction.11, "results/DESeqResults.interaction_d11.rds")
saveRDS(results.interaction33, "results/DESeqResults.interaction_d33.rds")
