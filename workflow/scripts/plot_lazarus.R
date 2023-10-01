library(tidyr)
library(readr)
library(dplyr)
library(ggseqlogo)
library(ggplot2)

dir_name <- unlist(snakemake@input)
output_file <- unlist(snakemake@output)
positions <- unlist(snakemake@params["positions"])

dat_files <- Sys.glob(file.path(dir_name, "node*.dat"))
nodes <- basename(dat_files) %>%
    gsub("node|.dat", "", .)

freqs <- lapply(dat_files, read_delim, delim = "  ", col_names = c("pos", "res"), show_col_types = F) %>%
    setNames(nodes) %>%
    bind_rows(.id = "node") %>%
    mutate(node = as.numeric(node)) %>%
    separate_rows(res) %>%
    mutate(frequency = lead(res)) %>%
    filter(grepl("[A-Z]", res)) %>%
    mutate(frequency = as.numeric(frequency)) %>%
    filter(pos %in% positions) %>%
    spread(res, frequency, fill = 0) %>%
    split(.$node) %>%
    lapply(select, -pos, -node) %>%
    lapply(t)

p <- ggseqlogo(freqs, ncol = 4, method = 'prob')

ggsave(output_file, p)
