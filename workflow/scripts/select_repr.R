
library(dplyr)
library(tidyr)
library(treeio)

with(snakemake@input, {
    ingroup_file <<- ingroup
    clstr_file <<- clstr
})

read.cdhit.clstr <- function(fname) {
    data.fields <- c("E.Value", "Aln", "Identity")
    read.table(fname, sep = "\t", comment.char = "", quote = "", fill = T, stringsAsFactors = F, col.names = c("Col1", "Col2")) %>%
        separate(Col1, into = c("Seq.Num", "Cluster"), sep = " ", fill = "right") %>%
        fill(Cluster) %>%
        filter(!grepl(">", Seq.Num)) %>%
        separate(Col2, into = c("Seq.Len", "Col2"), sep = "aa, >", convert = T) %>%
        extract(Col2, into = c("Seq.Name", "Is.Representative", "Col2"), regex = "(.*?)[.]{3} ([*]|at) ?(.*)") %>%
        mutate(Is.Representative = Is.Representative == "*", Col2 = ifelse(Is.Representative, "100%", Col2)) %>%
        group_by(Cluster) %>%
        mutate(Representative = Seq.Name[which(Is.Representative)]) %>%
        separate_rows(Col2, sep = ",") %>%
        separate(Col2, into = data.fields, sep = "/", fill = "left", convert = T) %>%
        mutate(Identity = sub("%", "", Identity) %>% as.numeric) %>%
        group_by(Seq.Name) %>%
        mutate(level.rank = paste0(".", 1:n() - 1), level.rank = ifelse(level.rank == ".0", "", level.rank)) %>%
        pivot_wider(names_from = level.rank, values_from = data.fields, names_sep = "") %>%
        ungroup
}

ingroup <- read.fasta(ingroup_file)

reprs <- read.cdhit.clstr(clstr_file) %>%
    mutate(ingroup_is_identical = Seq.Name %in% names(ingroup) & Identity == 100) %>%    
    group_by(Representative) %>%
    arrange(Identity) %>%
    summarize(Seq.Name = ifelse(any(ingroup_is_identical), first(Seq.Name), Representative)) %>%
    pull(Seq.Name) %>%
    c(names(ingroup)) %>%
    unique
cat(reprs, sep = "\n", file = unlist(snakemake@output))
