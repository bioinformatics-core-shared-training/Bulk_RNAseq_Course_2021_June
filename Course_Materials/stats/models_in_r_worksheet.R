# ==== Worksheet for Linear Models Section ==== #

# load packages
library(DESeq2)
library(tidyverse)


# Model Specification - Formula Syntax ------------------------------------

# create example data
example_samples <- tibble(
  sample = paste0("sample", 1:6),
  treatment = rep(c("A", "B"), each = 3)
)
example_samples

# define a model formula
temperature_model <- formula(~ treatment)

# model matrix - creates indicator variables for us
model.matrix(temperature_model, data = example_samples)



# Exercise 2 --------------------------------------------------------------


# Use read_tsv() to read table in "data/samplesheet_corrected.tsv"
# Store it in an object called sample_info


# Create a model formula to investigate the effect of “Status” on gene expression.


# Look at the model matrix and identify which is the reference group in the model.




# Exercise 3 --------------------------------------------------------------


# Using sample_info, create a new design formula specifying an additive model between “Status” and “TimePoint”.


# How many coefficients do you have with this model?


# What is your reference group?




# Models in DESeq2 --------------------------------------------------------

# simulate example DESeqDataSeq object
dds <- makeExampleDESeqDataSet(n = 100, m = 9)
colData(dds) <- DataFrame(treatment = factor(rep(c("A", "B", "C"), each = 3)),
                          row.names = paste0("sample", 1:9))

# check the sample information
colData(dds)

# add design to DESeqDataSet object
design(dds) <- formula(~ treatment)

# check model matrix
model.matrix(design(dds), data = colData(dds))

# fit the DESeq statistical model
dds <- DESeq(dds)

# obtain the results of differential expression
results(dds)

# check what coefficients are available to us
resultsNames(dds)




# Exercise 5 --------------------------------------------------------------

# simulate example DESeqDataSeq object
dds <- makeExampleDESeqDataSet(n = 100, m = 12)
colData(dds) <- DataFrame(treatment = factor(rep(c("A", "B"), each = 6)),
                          genotype = factor(rep(c("WT", "MUT"), 6),
                                            levels = c("WT", "MUT")),
                          row.names = paste0("sample", 1:12))

# check the sample information
colData(dds)

# fix the code below to model the effects of genotype, treatment and their interaction.
design(dds) <- formula(~ FIXME)

# fit the model
dds <- DESeq(dds)

# use resultsNames() to see the coefficients that DESeq created
# relate them to the "beta" parameters shown in the exercise page
