# Set working directory
setwd("E:/RNA Seq Data Analysis/siva_output_files")

# Load libraries
library(readr)
library(tidyr)
library(tximport)
library(GenomicFeatures)

# Get the list of directories
dirs <- list.files("siva_SALMON_OUT/WT")
print(dirs)

# Get the list of quant files
quant_files <- list.files("siva_SALMON_OUT/WT", pattern = "quant.sf", recursive = TRUE, full.names = TRUE)
names(quant_files) <- dirs
print(head(quant_files))

# Perform tximport with type = "salmon" and txOut = TRUE
txi.tx <- tximport(quant_files, type = "salmon", txOut = TRUE)

# Read the first quant file
quants <- read_tsv(quant_files[1])

# Create txdb from gtf file
gtf_file <- "transcriptome.gtf"
if (file.exists(gtf_file)) {
  txdb <- makeTxDbFromGFF(gtf_file)
  print(keytypes(txdb))
  print(columns(txdb))
  
  # Get keys from txdb
  k <- keys(txdb, keytype = "TXNAME")
  
  # Generate tx2gene mapping
  tx_map <- select(txdb, keys = k, columns = "GENEID", keytype = "TXNAME")
  head(tx_map)
  
  # Save tx2gene mapping to a csv file
  write.csv(tx_map, file = "tx2gene-WT.csv", row.names = FALSE, quote = FALSE)
  
  # Perform tximport with lengthScaledTPM
  txi_lengthscaledTPMvalues <- tximport(quant_files, type = "salmon", txIn = TRUE,
                                        countsFromAbundance = "lengthScaledTPM",
                                        tx2gene = tx_map, dropInfReps = TRUE)
  
  # Save txi_lengthscaledTPMvalues to a csv file
  write.csv(txi_lengthscaledTPMvalues, file = "txi_lengthscaledTPM_WT.csv",
            row.names = TRUE, quote = TRUE)
} else {
  print("GTF file not found.")
}
