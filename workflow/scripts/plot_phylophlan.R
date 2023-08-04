snakemake@source("tree_functions.R")

with(snakemake@input, {
    gtdb_file <<- gtdb
    lanclos_file <<- lanclos
    oceandna_file <<- oceandna
    gene_files <<- genes
    tree_file <<- tree
})
output_file <- unlist(snakemake@output)
with(snakemake@params, {
    genera <<- c(genera, outgroup) %>%
        {setNames(names(.), gsub("g__", "", .))}
    outgroup <<- outgroup
    show_labels <<- show_labels
})

if (interactive()) {
    gene_files <- Sys.glob("analysis/diamond_collect/*.tsv")
    tree_file <- "analysis/IIIa_phylophlan/IIIa.tre.treefile"
    gtdb_file <- "analysis/metadata/gtdb_filtered.tsv"
    lanclos_file <- "analysis/metadata/lanclos.tsv"
    oceandna_file <- "analysis/metadata/oceandna_filtered.tsv"
    output_file <- "tmp.pdf"
    genera <- c(HIMB114 = "IIIa.1", `GCA-002704185` = "IIIa.2", IMCC9063 = "IIIa.3", Fonsibacter = "IIIb")
    outgroup <- c(IIIb = "Fonsibacter")
    show_labels <- "LSUCC|HIMB114|IMCC9063"
}
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
outgroup_tips <- filter(taxonomy, clade %in% names(outgroup), grepl("^GC", genome)) %>%
    pull(genome)

genes <- lapply(gene_files, read.table, sep = "\t", header = T, na.strings = "") %>%
    bind_rows %>%
    distinct(genome, gene) %>%
    mutate(x = as.numeric(as.factor(gene)))

tree <- read.tree(tree_file) %>%
    ape::root(outgroup_tips, edgelabel = T) %>%
    drop.tip(outgroup_tips) %>%
    as_tibble %>%
    left_join(taxonomy, by = c(label = "genome")) %>%
    mutate(support = ifelse(node %in% parent & node != parent, label, NA)) %>%
    separate(support, into = c("SH_aLRT", "UFboot"), sep = "/", convert = T) %>%
    mutate(name_show = ifelse(grepl(show_labels, Name), Name, NA)) %>%
    add_mrca(clade_label) %>%
    add_mrca(species)

p <- ggtree(to_treedata(tree), layout = "rectangular") +
    geom_facet(aes(x = x, color = gene), genes, geom_tile, "genes") +
    geom_point2(aes(subset = !is.na(UFboot) & UFboot >= 90, x = branch), color = "gray", size = 1) +
    geom_tiplab(aes(subset = !is.na(name_show), label = name_show, offset = 1)) +
    geom_cladelab(mapping = aes(subset = !is.na(clade_label_mrca), node = node, label = clade_label_mrca), offset = 0.05) +
    #geom_cladelab(mapping = aes(subset = !is.na(species_mrca), node = node, label = species_mrca), offset = 0.05) +
    geom_treescale(width = 0.2) +
    theme(
        strip.background = element_blank(),
        strip.text.x = element_blank()
    )

ggsave(output_file, facet_widths(p, widths = c(10, 1)), width = 3, height = 7)
