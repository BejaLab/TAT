from Bio import SeqIO

# Input and output file paths
fasta_file = str(snakemake.input)
nexus_file = str(snakemake.output)
codons = snakemake.params["codons"]

# Read the FASTA file
records = list(SeqIO.parse(fasta_file, "fasta"))

num_records = len(records)
alignment_length = len(records[0].seq)

# Create a Nexus file
with open(nexus_file, "w") as nexus_out:
    nexus_out.write("#NEXUS\n\n")
    nexus_out.write("begin data;\n")
    nexus_out.write("dimensions ntax={} nchar={};\n".format(num_records, alignment_length))
    nexus_out.write("format datatype=dna missing=n gap=-;\n")
    nexus_out.write("matrix\n")

    # Write the sequences putting record IDs in single quotes
    for record in records:
        nexus_out.write("    '{}'\t{}\n".format(record.id, record.seq))

    nexus_out.write(";\n")
    nexus_out.write("end;\n")

    # Add SETS block with partition definitions, using "codon" as the prefix
    nexus_out.write("begin sets;\n")
    for codon in codons:
        nexus_out.write("    charset codon{} = {}-{}\\3;\n".format(codon, codon, alignment_length))
    nexus_out.write("end;\n")

