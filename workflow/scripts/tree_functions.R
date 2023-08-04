library(dplyr)
library(tidyr)
library(ggtree)
library(treeio)
library(phangorn)
library(ggplot2)
library(tibble)
library(igraph)

get_mrca <- function(phylo, tips) {
    getMRCA(phylo, tips) %>%
        replace(is.null(.), NA)
}

add_mrca <- function(tree, colname) {
    colname <- deparse(substitute(colname))
    treedata <- to_treedata(tree)
    mrca <- mutate(tree, my_column = !!as.name(colname)) %>%
        mutate(my_column = ifelse(my_column == "", NA, my_column)) %>%
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
