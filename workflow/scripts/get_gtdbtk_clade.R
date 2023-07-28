
library(ape)

tree_file <- unlist(snakemake@input)
with(snakemake@output, {
    clade_file <<- clade
    tips_file <<- tips
})
taxon <- snakemake@params["taxon"]

tree <- read.tree(tree_file)
parent <- with(tree, node.label[grepl(taxon, node.label)])
clade <- extract.clade(tree, parent)

write.tree(clade, file = clade_file)
cat(clade$tip.label, sep = "\n", file = tips_file)
