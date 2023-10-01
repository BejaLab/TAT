
configfile: "config.yaml"

xmls ,= glob_wildcards("analysis/beast2/{xml}.xml")
xmls = [ "beast_linked_models" ]

rule all:
    input:
        expand("output/rhodopsin_tree-{xml}.pdf", xml = xmls),
        "output/lazarus.pdf"

rule beast:
    input:
        "analysis/beast2/{xml}.xml"
    output:
        "analysis/beast2/{xml}-codon12.trees",
        "analysis/beast2/{xml}.xml.state"
    log:
        "analysis/beast2/{xml}.log"
    threads:
        10
    shell:
        "beast -threads {threads} -working {input}"

rule root_annot:
    input:
        trees = "analysis/beast2/{xml}-codon12.trees",
        bindir = "RootAnnotator/dependencies/unix"
    output:
        "analysis/beast2/{xml}-rootAnnotator_annotatedMCCTree.nexus"
    params:
        burnin = 3000
    shadow:
        "minimal"
    shell:
        "{input.bindir}/rootAnnotator -s {input.trees} -a {input.bindir}/treeAnnotator/bin/treeAnnotator --computeMccTree --burnin {params.burnin} && mv rootAnnotator_annotatedMCCTree.nexus {output}"

rule fix_nexus:
    input:
        "analysis/beast2/{xml}-rootAnnotator_annotatedMCCTree.nexus"
    output:
        "analysis/beast2/{xml}-rootAnnotator_annotatedMCCTree.nexus_fixed"
    shell:
        "sed -E 's/Root_Probability=[0-9]*: */Root_Probability=/g' {input} | cat - <(echo) > {output}"

rule nexus_to_newick:
    input:
        "analysis/beast2/{xml}-rootAnnotator_annotatedMCCTree.nexus_fixed"
    output:
        "analysis/beast2/{xml}-rootAnnotator_annotatedMCCTree.nwk"
    conda:
        "envs/dendropy.yaml"
    script:
        "scripts/nexus_to_newick.py"

rule cp_fasta:
    input:
        "analysis/TAT/TwR.mafft"
    output:
        "analysis/lazarus/TwR.mafft.fasta"
    shell:
        "cp {input} {output}"

rule lazarus:
    input:
        fasta = "analysis/lazarus/TwR.mafft.fasta",
        tree = "analysis/beast2/beast_linked_models.raxml.bestTree.rooted"
    output:
        directory("analysis/lazarus/tree1/")
    params:
        asrv = 4,
        model = "lg"
    log:
        "analysis/lazarus/TwR.log"
    conda:
        "envs/lazarus.yaml"
    shell:
        """
        lazarus.py --alignment {input.fasta} --tree {input.tree} --codeml --model "$(realpath "$CONDA_PREFIX/dat/{params.model}.dat")" --asrv {params.asrv} --fix_asrv False --branch_lengths fixed --outputdir $(dirname {output}) --verbose {log}
        """

rule plot_lazarus:
    input:
        "analysis/lazarus/tree1/"
    output:
        "output/lazarus.pdf"
    params:
        positions = [ 82, 86, 93 ]
    conda:
        "envs/ggseqlogo.yaml"
    script:
        "scripts/plot_lazarus.R"

rule raxml_evaluate:
    input:
        fasta = "analysis/TAT/TwR.mafft",
        tree = "analysis/beast2/{xml}-rootAnnotator_annotatedMCCTree.nwk"
    output:
        "analysis/beast2/{xml}.raxml.bestTree",
        "analysis/beast2/{xml}.raxml.bestModel"
    params:
        prefix = "analysis/beast2/{xml}",
        seed = 123,
        model = "LG+G"
    log:
        "analysis/beast2/{xml}.raxml.log"
    conda:
        "envs/raxml-ng.yaml"
    shell:
        "raxml-ng --redo --evaluate --msa {input.fasta} --tree {input.tree} --model {params.model} --seed {params.seed} --prefix {params.prefix} &> {log}"

rule raxml_mad:
    input:
        "analysis/beast2/{xml}.raxml.bestTree"
    output:
        "analysis/beast2/{xml}.raxml.bestTree.rooted"
    shell:
        "mad {input}"

rule plot_rhodopsins:
    input:
        tree = "analysis/beast2/{xml}-rootAnnotator_annotatedMCCTree.nexus_fixed",
        gtdb = "analysis/metadata/gtdb_filtered.tsv",
        lanclos = "analysis/metadata/lanclos.tsv",
        oceandna = "analysis/metadata/oceandna_filtered.tsv",
        aln = "analysis/TAT/rhodopsins.mafft",
        pos = "input/pos.txt",
        ingroup = "input/ingroup.tsv",
        habitats = "input/habitats.tsv"
    output:
        "output/rhodopsin_tree-{xml}.pdf"
    params:
        genera = config["genera"]
    conda:
        "envs/ggtree.yaml"
    script:
        "scripts/plot_rhodopsins.R"
