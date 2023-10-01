
library(dplyr)
library(tidyr)

if (interactive()) {
    tree_file <- "beast/rootAnnotator_annotatedMCCTree.nexus_fixed"
    gtdb_file <- "analysis/metadata/gtdb_filtered.tsv"
    lanclos_file <- "analysis/metadata/lanclos.tsv"
    oceandna_file <- "analysis/metadata/oceandna_filtered.tsv"
    output_file <- "tmp.pdf"
    genera <- c(HIMB114 = "IIIa.1", `GCA-002704185` = "IIIa.2", IMCC9063 = "IIIa.3")
    lineage <- "d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Pelagibacterales;f__Pelagibacteraceae"
    aln_file <- "analysis/TAT/rhodopsins.mafft"
    pos_file <- "input/pos.txt"
    ingroup_file <- "input/ingroup.tsv"
    ingroup_cds_file <- "input/ingroup.fna"
    genome_PR_cds_file <- "analysis/PR/genomes.cds"
    genome_TAT_cds_file <- "analysis/TAT/genomes.cds"
    habitats_file <- "input/habitats.tsv"
} else {
    with(snakemake@input, {
        tree_file <<- tree
        gtdb_file <<- gtdb
        lanclos_file <<- lanclos
        oceandna_file <<- oceandna
        aln_file <<- aln
        pos_file <<- pos
        ingroup_file <<- ingroup
    })
    with(snakemake@params, {
        lineage <<- lineage
        genera <<- genera %>%
            {setNames(names(.), gsub("g__", "", .))}
    })
    output_file <- unlist(snakemake@output)
}

gen_taxonomy <- sprintf("%s;g__%s;s__", lineage, names(genera)) %>%
    setNames(genera)

oceandna <- read.table(oceandna_file, sep = "\t", col.names = c("genome", "taxonomy")) %>%
    mutate(source = "OceanDNA")
lanclos <- read.table(lanclos_file, sep = "\t", header = T, fill = T, comment.char = "") %>%
    extract(GTDB_accession, into = "acc", regex = "_(\\d+)", remove = F) %>%
    mutate(lanclos_source = "Lanclos23")
ingroup <- read.table(ingroup_file, sep = "\t", header = T, na.strings = "") %>%
    select(label = Accession, alias = Alias, expressed = Expressed, ingroup_taxonomy = taxonomy, ingroup_source = source)
genomes <- read.table(gtdb_file, sep = "\t", col.names = c("genome", "taxonomy", "biosamples")) %>%
    mutate(source = "GTDB/NCBI") %>%
    extract(genome, into = "acc", regex = "_(\\d+)", remove = F) %>%
    full_join(lanclos, by = "acc") %>%
    bind_rows(oceandna) %>%
    mutate(source = ifelse(is.na(source), lanclos_source, source)) %>%
    mutate(genome = ifelse(is.na(genome), Name, genome))

pos <- readLines(pos_file)
ref <- pos[1]
dte <- as.numeric(pos[-1])

aln <- read.fasta(aln_file, type = "AA") %>%
    as.character %>%
    lapply(toupper)

read_cds <- function(filename) {
    read.fasta(filename) %>%
        as.character %>%
        lapply(toupper) %>%
        lapply(paste, collapse = "") %>%
        {data.frame(label = sub(" .+", "", names(.)), cds = unname(unlist(.)))}
}

cds <- list(TwR = read_cds(ingroup_cds_file), TwR = read_cds(genome_TAT_cds_file), PR = read_cds(genome_PR_cds_file)) %>%
    bind_rows(.id = "family")

ref_pos <- data.frame(res = aln[[ref]]) %>%
    mutate(aln_pos = 1:n()) %>%
    filter(res != "-") %>%
    filter(1:n() %in% dte) %>%
    pull(aln_pos)

motifs <- lapply(aln, `[`, ref_pos) %>%
    lapply(paste, collapse = "") %>%
    {data.frame(label = sub(" .+", "", names(.)), motif = unname(unlist(.)))}

fst <- function(x) first(na.omit(x))

genes <- lapply(aln, paste, collapse = "") %>%
    {data.frame(description = names(.), aa = unname(unlist(.)))} %>%
    mutate(label = sub(" .+", "", description)) %>%
    mutate(aa = gsub("-", "", aa)) %>%
    separate(label, into = c("genome", "gene"), sep = "@", fill = "left", remove = F) %>%
    extract(description, into = c("start", "end", "strand"), regex = "# (\\d+) # (\\d+) # (-?1)", convert = T, remove = F) %>%
    extract(description, into = c("left_partial", "right_partial"), regex = "partial=(\\d)(\\d);", convert = T) %>%
    left_join(motifs, by = "label") %>%
    left_join(cds, by = "label") %>%
    left_join(ingroup, by = "label") %>%
    left_join(genomes, by = "genome") %>%
    mutate(source = ifelse(is.na(ingroup_source), source, ingroup_source)) %>%
    mutate(taxonomy = ifelse(is.na(taxonomy), ingroup_taxonomy, taxonomy)) %>%
    extract(taxonomy, into = "genus", regex = "g__([^;]+);", remove = F) %>%
    mutate(IIIa_subgroup = recode(genus, !!!as.list(genera), .default = NA_character_)) %>%
    mutate(lanclos_taxonomy = recode(SubcladeA, !!!gen_taxonomy, .default = NA_character_)) %>%
    mutate(genome = ifelse(is.na(genome), GTDB_accession, genome)) %>%
    mutate(taxonomy = ifelse(is.na(taxonomy), lanclos_taxonomy, taxonomy)) %>%
    mutate(subgroup = ifelse(is.na(SubcladeA), IIIa_subgroup, SubcladeA)) %>%
    mutate(subgroup = ifelse(subgroup %in% genera, subgroup, NA_character_)) %>%
    replace_na(list(left_partial = 0, right_partial = 0)) %>%
    filter(!is.na(family)) %>%
    group_by(genome, aa, cds, family, left_partial, right_partial, motif) %>%
    arrange(is.na(alias)) %>%
    summarize(label = first(label), alias = fst(alias), expressed = fst(expressed), genome_name = fst(Name), taxonomy = fst(taxonomy), subgroup = fst(subgroup), source = fst(source)) %>%
    select(label, family, alias, expressed, left_partial, right_partial, motif, genome, genome_name, taxonomy, source, subgroup, aa, cds)
write.table(genes, file = output_file, sep = "\t", row.names = F, col.names = T, na = "")
