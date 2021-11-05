setwd("E:/RNA Seq Data Analysis/siva_output_files")
library(readr)
library(tidyr)

library("tximport", lib.loc="~/R/win-library/3.5")
dirs <- list.files("siva_SALMON_OUT/WT")
dirs
quant_files <- list.files("siva_SALMON_OUT/WT",pattern="quant.sf",recursive = TRUE,full.names = TRUE)
names(quant_files) <- dirs
head(quant_files)
txi.tx <- tximport(quant_files, type = "salmon", txOut = TRUE)
quants <- read_tsv(quant_files[1])

library(GenomicFeatures)
gtf_file <- "transcriptome.gtf"
file.exists(gtf_file)
txdb <- makeTxDbFromGFF(gtf_file)
keytypes(txdb)
columns(txdb)
k <- keys(txdb, keytype="TXNAME")

library(tximport)
tx_map <- select(txdb, keys = k, columns="GENEID", keytype = "TXNAME")
head(tx_map)
tx2gene <- tx_map
write.csv(tx2gene,file="tx2gene-WT.csv",row.names = FALSE,quote=FALSE)

txi_lengthscaledTPMvalues = tximport::tximport(quant_files,type = "salmon",txIn = TRUE, countsFromAbundance = "lengthScaledTPM",tx2gene = tx2gene,dropInfReps = TRUE)
write.csv(txi_lengthscaledTPMvalues,file="txi_lengthscaledTPM_WT.csv",row.names = TRUE,quote=TRUE)
