
library(dplyr)
library(tidyr)

if (interactive()) {
    tree_file <- "analysis/beast2/rootAnnotator_annotatedMCCTree.nexus_fixed"
    gtdb_file <- "analysis/metadata/gtdb_filtered.tsv"
    lanclos_file <- "analysis/metadata/lanclos.tsv"
    oceandna_file <- "analysis/metadata/oceandna_filtered.tsv"
    output_file <- "tmp.pdf"
    genera <- c(HIMB114 = "IIIa.1", `GCA-002704185` = "IIIa.2", IMCC9063 = "IIIa.3")
    aln_file <- "analysis/TAT/rhodopsins.mafft"
    pos_file <- "input/pos.txt"
    ingroup_file <- "input/ingroup.tsv"
    habitats_file <- "input/habitats.tsv"
    source("workflow/scripts/tree_functions.R")
} else {
    with(snakemake@input, {
        tree_file <<- tree
        gtdb_file <<- gtdb
        lanclos_file <<- lanclos
        oceandna_file <<- oceandna
        aln_file <<- aln
        pos_file <<- pos
        ingroup_file <<- ingroup
        habitats_file <<- habitats
    })
    with(snakemake@params, {
        genera <<- genera %>%
            {setNames(names(.), gsub("g__", "", .))}
    })
    output_file <- unlist(snakemake@output)
    snakemake@source("tree_functions.R")
}
habitats <- read.table(habitats_file, sep = "\t", comment.char = "", header = T)

pos <- readLines(pos_file)
ref <- pos[1]
dte <- as.numeric(pos[-1])

aln <- read.fasta(aln_file, type = "AA") %>%
    as.character %>%
    `names<-`(sub(" .+", "", names(.)))

motifs <- data.frame(res = aln[[ref]]) %>%
    mutate(aln_pos = 1:n()) %>%
    filter(res != "-") %>%
    filter(1:n() %in% dte) %>%
    pull(aln_pos) %>%
    lapply(aln, `[`, .) %>%
    lapply(toupper) %>%
    lapply(t) %>%
    lapply(data.frame) %>%
    lapply(setNames, c("D85", "T89", "D96")) %>%
    bind_rows(.id = "label")

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
    left_join(habitats, by = "label") %>%
    distinct(label, .keep_all = T)

tree <- read.beast(tree_file) %>%
    as_tibble %>%
    left_join(metadata, by = "label") %>%
    add_mrca(clade_label)

wrap_float <- function(x) {
    recode(format(round(x, 2), nsmall = 2), "  NA" = "")
}

p <- ggtree(to_treedata(tree), layout = "rectangular") +
    geom_text(aes(x = branch, label = wrap_float(posterior)), color = "black", size = 2, vjust = -0.5) +
    geom_text(aes(x = branch, label = wrap_float(Root_Probability)), color = "red", size = 2, vjust = 1.5) +
    geom_tippoint(aes(subset = !is.na(Expressed), x = x + 0.025), color = "red") +
    geom_tippoint(aes(subset = !is.na(type), color = type, x = x + 0.05)) +
    geom_tiplab(aes(subset = !is.na(Alias), label = Alias), offset = 0.1) +
    geom_tiplab(aes(label = D85), align = T, offset = 0.4) +
    geom_tiplab(aes(label = T89), align = T, offset = 0.45) +
    geom_tiplab(aes(label = D96), align = T, offset = 0.5) +
    geom_cladelab(mapping = aes(subset = !is.na(clade_label_mrca), node = node, label = clade_label_mrca), align = T, offset = 0.6) +
    geom_treescale(width = 0.1) +
    theme(legend.position = "bottom")

ggsave(output_file, p, width = 3.5, height = 3.5)
