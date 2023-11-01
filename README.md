# TAT-like or Twin-peaked rhodopsins

This repository contains code ana analysis results used in the paper Mannen et al (2023) [Multiple roles of a conserved glutamate residue for unique biophysical properties in a new group of microbial rhodopsins homologous to TAT rhodopsin](https://doi.org/10.1016/j.jmb.2023.168331).

The workflow is divided into two parts: most of the analyses are covered by the default [snakemake](https://snakemake.readthedocs.io/) file `workflow/Snakefile` while the [beast2](https://www.beast2.org/) analysis is in `workflow/Beast2.snakefile`.

## Dependencies

Most of the dependencies of the main workflow pipilene are taken care of with [conda](https://docs.conda.io/), but in addition to snakemake and conda, the following dependencies have to be installed manually:

* [usearch](https://www.drive5.com/usearch/) is expected to be available from the `PATH`;
* [RootAnnotator](https://sourceforge.net/projects/rootannotator/) is expected to be in the directory named `RootAnnotator` in the current directory;
* [mad](https://www.mikrobio.uni-kiel.de/de/ag-dagan/ressourcen) is expected to be available from the `PATH`.

## Files in this repository

The workflow files are located in `workflow`. Input files to run the pipeline(s) from scretch are in `input`.

Final output files are in the folder `output`. Immediate analysis results needed to produce them are included as well:

* `analysis/IIIa_phylophlan/IIIa.tre.treefile` -- results of phylogenetic analysis of Pelagibacterales subclade IIIa
* `analysis/diamond_collect/{gtdb,lanclos,oceandna}.tsv` -- tsv files summarizing presence of rhodopsins in Pelagibacterales genomes obtained from three sources
* `analysis/TAT/rhodopsins.mafft` -- fasta files with alignment of all of the collected rhodopsins
* `analysis/metadata/{gtdb_filtered,lanclos,oceandna_filtered}.tsv` -- metadata for the analyzed Pelagibacterales genomes
* `analysis/metadata/gtdb_clade.nwk` -- tree in newick format corresponding to the o\_\_Pelagibacterales clade in GTDB r. 214.1
* `analysis/beast2/{beast_linked_models-codon12.trees,beast_linked_models.xml.state,beast_linked_models.log}` -- results of the beast2 run
* `analysis/beast2/beast_linked_models-rootAnnotator_annotatedMCCTree.nexus_fixed` -- (fixed) output of rootAnnotator, tree in nexus format
* `analysis/lazarus` -- lazarus analysis

Use [Issue tracker](https://github.com/BejaLab/TAT/issues) for questions/requests.
