
library(dplyr)
library(tidyr)
library(tools)

input_files <- unlist(snakemake@input) %>%
    setNames(file_path_sans_ext(basename(.)))
output_file <- unlist(snakemake@output)

evalue_threshold <- snakemake@params["evalue"]
pident_threshold <- snakemake@params["pident"]

file_sizes <- file.info(input_files) %>%
    pull(size) %>%
    setNames(names(input_files))

col_names <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
names(file_sizes[file_sizes > 0]) %>%
    `[`(input_files, .) %>%
    lapply(read.table, col.names = col_names) %>%
    bind_rows(.id = "genome") %>%
    separate(qseqid, into = "gene", extra = "drop", sep = "_") %>%
    arrange(-bitscore) %>%
    distinct(sseqid, .keep_all = T) %>%
    filter(evalue < evalue_threshold, pident > pident_threshold) %>%
    select(genome, gene, protein = sseqid) %>%
    write.table(file = output_file, na = "", quote = F, sep = "\t", col.names = T, row.names = F)
