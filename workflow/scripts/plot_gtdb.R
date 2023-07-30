library(dplyr)
library(tidyr)
library(ggtree)
library(treeio)
library(phangorn)
library(ggplot2)
library(tibble)
library(igraph)

tsv_file <- snakemake@input["tsv"]
tree_file <- snakemake@input["tree"]
output_file <- unlist(snakemake@output)

# dist_threshold <- 0.1

#genes_file <- "analysis/diamond_collect/gtdb.tsv"
#tree_file <- "analysis/metadata/gtdb_clade.nwk"

#gtdb_file <- "analysis/metadata/gtdb_filtered.tsv"
#lanclos_file <- "analysis/metadata/lanclos.tsv"
#output_file <- "tmp.pdf"

lanclos <- read.table(lanclos_file, sep = "\t", header = T, fill = T, comment.char = "") %>%
    extract(GTDB_accession, into = "acc", regex = "_(\\d+)") %>%
    filter(!is.na(acc))
taxonomy <- read.table(gtdb_file, sep = "\t", col.names = c("genome", "taxonomy", "biosamples")) %>%
    extract(genome, into = "acc", regex = "_(\\d+)", remove = F) %>%
    left_join(lanclos, by = "acc") %>%
    extract(taxonomy, into = c("family", "genus"), regex = "f__([^;]+);g__([^;]+);") %>%
    select(acc, genome, family, genus, clade = SubcladeA) %>%
    mutate(clade = case_when(grepl("III", clade) ~ gsub("[.].+", "", clade), T ~ gsub("[a-z]", "", clade))) %>%
    group_by(genus) %>%
    mutate(clade = first(na.omit(clade))) %>%
    group_by(family) %>%
    mutate(clade = case_when(n_distinct(clade, na.rm = T) == 0 & is.na(clade) ~ paste0("f__", family), T ~ clade))

genes <- read.table(genes_file, sep = "\t", header = T, na.strings = "") %>%
    distinct(genome, gene) %>%
    mutate(x = as.numeric(as.factor(gene)))

get_mrca <- function(phylo, tips) {
    getMRCA(phylo, tips) %>%
        replace(is.null(.), NA)
}

add_mrca <- function(tree, colname) {
    colname <- deparse(substitute(colname))
    treedata <- to_treedata(tree)
    mrca <- mutate(tree, my_column = !!as.name(colname)) %>%
        group_by(my_column) %>%
        mutate(is.tip = label %in% treedata@phylo$tip.label) %>%
        mutate(no_data = all(is.na(my_column))) %>%
        mutate(mrca = get_mrca(treedata@phylo, node[is.tip])) %>%
        mutate(mrca = ifelse(no_data | is.na(mrca), node, mrca)) %>%
        group_by(mrca) %>%
        mutate(enough_tips = sum(is.tip) > 1) %>%
        mutate(ifelse(node == mrca & enough_tips, first(na.omit(my_column)), NA)) %>%
        pull
    tree[[paste0(colname, "_mrca")]] <- mrca
    return(tree)
}

add_ancestor <- function(tree, colname) {
    colname <- deparse(substitute(colname))
    treedata <- to_treedata(tree)
    ancestors <- mutate(tree, my_column = !!as.name(colname)) %>%
        filter(!is.na(my_column)) %>%
        with(setNames(my_column, node))
    ancestor_nodes <- as.numeric(names(ancestors))
    descs <- Descendants(treedata@phylo, ancestor_nodes) %>%
        lapply(data.frame) %>%
        setNames(ancestor_nodes) %>%
        bind_rows(.id = "ancestor") %>%
        setNames(c("ancestor", "node")) %>%
        distinct(node, .keep_all = T) %>%
        mutate(my_column = ancestors[ancestor]) %>%
        with(setNames(my_column, node))
    tree[[paste0(colname, "_ancestor")]] <- descs[as.character(tree$node)]
    return(tree)
}

to_treedata <- function(tree) {
    class(tree) <- c("tbl_tree", "tbl_df", "tbl", "data.frame")
    as.treedata(tree)
}

full_tree <- read.tree(tree_file)
#keep_labels <- cophenetic(full_tree) %>%
#    as.data.frame %>%
#    rownames_to_column(var = "A") %>%
#    pivot_longer(cols = -A, names_to = "B", values_to = "distance") %>%
#    filter(A <= B, distance <= dist_threshold) %>%
#    graph_from_data_frame(directed = F) %>%
#    clusters %>%
#    with(data.frame(label = names(membership), cluster = membership)) %>%
#    distinct(cluster, .keep_all = T) %>%
#    pull(label)

#tree <- keep.tip(full_tree, keep_labels) %>%
tree <- full_tree %>%
    as_tibble %>%
    mutate(label = gsub("'", "", label)) %>%
    separate(label, into = c("support", "label"), sep = "-", extra = "merge", fill = "left") %>%
    mutate(label = ifelse(node %in% parent, label, substring(label, 4))) %>%
    #mutate(genus = ifelse(grepl("g__", label), label, NA)) %>%
    extract(label, into = "acc", regex = "_(\\d+)", remove = F) %>%
    left_join(taxonomy, by = "acc") %>%
    add_mrca(clade)

p <- ggtree(to_treedata(tree), layout = "rectangular") +
    geom_facet(aes(x = x, color = gene), genes, geom_tile, "genes") +
    geom_cladelab(mapping = aes(subset = !is.na(clade_mrca), node = node, label = clade_mrca), offset = 0.1) +
    geom_treescale() +
    theme(
        strip.background = element_blank(),
        strip.text.x = element_blank()
    )

ggsave(output_file, facet_widths(p, widths = c(30, 1)), width = 3, height = 7)
