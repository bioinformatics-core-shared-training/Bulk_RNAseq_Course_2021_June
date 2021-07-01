library(tximport)
library(DESeq2)
library(tidyverse)

# load sample metadata

sampleinfo <- read_tsv("data/samplesheet.tsv", col_types = "cccc")

# load the count data

files <- str_c("salmon/", sampleinfo$SampleName, "/quant.sf")
files <- set_names(files, sampleinfo$SampleName)

tx2gene <- read_tsv("references/tx2gene.tsv")
head(tx2gene)

txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

saveRDS(txi, file = "salmon_output/txi.rds")

# Exercise 1

tpm <- tximport(files, type = "salmon", 
                tx2gene = tx2gene,
                countsFromAbundance = "lengthScaledTPM")

# Create a raw counts matrix

rawCounts <- round(txi$counts, 0)

dim(rawCounts)

keep <- rowSums(rawCounts) > 5
table(keep)

filtCounts <- rawCounts[keep, ]
dim(filtCounts)

# Data transformation and visualisations

summary(filtCounts)

## Boxplot

boxplot(filtCounts)

## Variance and counts

plot(rowMeans(filtCounts), rowSds(filtCounts),
     xlim = c(0, 10000),
     ylim = c(0, 5000))

# log2 transformation

logcounts <- log2(filtCounts + 1)

statusCols <- str_replace_all(sampleinfo$Status, c(Infected="red", Uninfected="orange"))

## Distribution

boxplot(logcounts,
        col=statusCols,
        las=2)

## mean v variance

plot(rowMeans(logcounts), rowSds(logcounts))

# DESeq2: Variance Stabilising Transformation

vstcounts <- vst(filtCounts)

## boxplot

boxplot(vstcounts,
        col = statusCols,
        las = 2)

## mean v variance

plot(rowMeans(vstcounts), rowSds(vstcounts))

# Exercise 2

rlogcounts <- rlog(filtCounts)
boxplot(rlogcounts,
        col = statusCols,
        las = 2)

# Principal Component Analysis

library(ggfortify)

rlogcounts <- rlog(filtCounts)


pcDat <- prcomp(t(rlogcounts))

autoplot(pcDat)

# add colour and shape

autoplot(pcDat,
         data = sampleinfo,
         colour = "Status",
         shape = "TimePoint",
         size = 5)

# Exercise 3

autoplot(pcDat,
         data = sampleinfo,
         colour = "Status",
         shape = "TimePoint",
         x = 2,
         y = 3,
         size = 5)

# Identify the swapped samples

library(ggrepel)

autoplot(pcDat,
         data = sampleinfo,
         colour = "Status",
         shape = "TimePoint",
         size = 5) +
  geom_text_repel(aes(x = PC1, y = PC2, label = SampleName))

# fix the sample labelling

sampleinfo <- mutate(sampleinfo, Status = case_when(
                                      SampleName=="SRR7657882" ~ "Uninfected",
                                      SampleName=="SRR7657873" ~ "Infected",
                                      TRUE ~ Status))


autoplot(pcDat,
         data = sampleinfo,
         colour = "Status",
         shape = "TimePoint",
         size = 5) +
  geom_text_repel(aes(x = PC1, y = PC2, label = SampleName))


# save the corrected sample sheet

write_tsv(sampleinfo, "results/SampleInfo_Corrected.txt")




