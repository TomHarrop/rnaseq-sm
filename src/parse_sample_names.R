#!/usr/bin/env Rscript

library(data.table)

set.seed(42)

# get a data.table of file names
read_dir = "data/reads"
read_file_table <- data.table(
    path = list.files(read_dir,
           pattern = "fastq.gz",
           full.names = TRUE,
           recursive = TRUE))

# split the path into components
read_file_table[, fullname := basename(path)]
read_file_table[, filename := sub(".fastq.gz", "", basename(path))]
read_file_table[, c("sample", "barcode", "lane", "read", "filenumber") :=
                    tstrsplit(filename, "_")]

# generate dummy sample names
sample_key <- read_file_table[, .(
    samplename = paste0(letters[sample(c(1:26), 4, replace = TRUE)],
                        collapse = "")
), by = sample]

# save sample names
fwrite(sample_key, "data/sample_key.csv", sep = ",")


unique(read_file_table[, .(barcode), by = sample])
