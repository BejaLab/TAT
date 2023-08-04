
snakemake@source("tree_functions.R")
library(dplyr)
library(tidyr)

if (interactive()) {
    tree_file <- "analysis/TAT/root_digger/TwR.mafft.cds.trim.raxml.bestTree.rooted.tree"
    gtdb_file <- "analysis/metadata/gtdb_filtered.tsv"
    lanclos_file <- "analysis/metadata/lanclos.tsv"
    oceandna_file <- "analysis/metadata/oceandna_filtered.tsv"
    output_file <- "tmp.pdf"
    genera <- c(HIMB114 = "IIIa.1", `GCA-002704185` = "IIIa.2", IMCC9063 = "IIIa.3")
    aln_file <- "analysis/TAT/rhodopsins.mafft"
    pos_file <- "input/pos.txt"
    ingroup_file <- "input/ingroup.tsv"
}
pos <- readLines(pos_file)
ref <- pos[1]
dte <- as.numeric(pos[-1])

aln <- read.fasta(aln_file, type = "AA") %>%
    as.character

motifs <- data.frame(res = aln[[ref]]) %>%
    mutate(aln_pos = 1:n()) %>%
    filter(res != "-") %>%
    filter(1:n() %in% dte) %>%
    pull(aln_pos) %>%
    lapply(aln, `[`, .) %>%
    lapply(paste, collapse = "") %>%
    {data.frame(label = names(.), motif = unlist(.))} %>%
    mutate(label = gsub(" .+", "", label), motif = toupper(motif))

oceandna <- read.table(oceandna_file, sep = "\t", col.names = c("genome", "taxonomy"))
lanclos <- read.table(lanclos_file, sep = "\t", header = T, fill = T, comment.char = "") %>%
    extract(GTDB_accession, into = "acc", regex = "_(\\d+)")
taxonomy <- read.table(gtdb_file, sep = "\t", col.names = c("genome", "taxonomy", "biosamples")) %>%
    extract(genome, into = "acc", regex = "_(\\d+)", remove = F) %>%
    full_join(lanclos, by = "acc") %>%
    bind_rows(oceandna) %>%
    extract(taxonomy, into = c("genus", "species"), regex = "g__([^;]+);s__([^;]*)") %>%
    mutate(genome = ifelse(is.na(genome), Name, genome)) %>%
    select(Name, genome, genus, species, clade = SubcladeA) %>%
    mutate(clade = recode(genus, !!!as.list(genera), .default = NA_character_, .missing = clade)) %>%
    group_by(clade) %>%
    mutate(clade_label = paste(clade, "=", genus))
ingroup <- read.table(ingroup_file, sep = "\t", header = T, na.strings = "")

metadata <- data.frame(label = names(aln)) %>%
    separate(label, into = c("genome", "gene"), sep = "@", fill = "left", remove = F) %>%
    left_join(taxonomy, by = "genome") %>%
    left_join(ingroup, by = c(label = "Accession")) %>%
    left_join(motifs, by = "label") %>%
    mutate(tree_label = gsub("[@]", "_", label)) %>%
    distinct(tree_label, .keep_all = T)

tree <- read.tree(tree_file) %>%
    phangorn::midpoint(node.labels = "support") %>%
    as_tibble %>%
    left_join(metadata, by = c(label = "tree_label")) %>%
    mutate(support = ifelse(node %in% parent & node != parent, label, NA)) %>%
    # separate(support, into = c("SH_aLRT", "UFboot"), sep = "/", convert = T) %>%
    add_mrca(clade_label)

p <- ggtree(to_treedata(tree), layout = "rectangular") +
    # geom_point2(aes(subset = !is.na(UFboot) & UFboot >= 90, x = branch), color = "gray", size = 1) +
    geom_tippoint(aes(subset = !is.na(Expressed)), color = "red") +
    geom_tiplab(aes(subset = !is.na(Alias), label = Alias), offset = 0.05) +
    geom_tiplab(aes(label = motif), align = T, offset = 0.6) +
    geom_cladelab(mapping = aes(subset = !is.na(clade_label_mrca), node = node, label = clade_label_mrca), align = T, offset = 0.8) +
    #geom_cladelab(mapping = aes(subset = !is.na(species_mrca), node = node, label = species_mrca), offset = 0.05) +
    geom_treescale(width = 0.2)

ggsave(output_file, p, width = 3, height = 5)
