from pandas import read_csv

oceandna_file = snakemake.input['oceandna']
gtdb_file = snakemake.input['gtdb']
genera = snakemake.params['genera']
output_file = str(snakemake.output)

oceandna = read_csv(oceandna_file, sep = "\t")

biosamples = set()

genera_vals = []
for genus in genera.values():
    genera_vals.append(f"{genus};")

def belongs_to_genus(tax):
    for genus in genera_vals:
        if genus in tax:
            return True
    return False

with open(gtdb_file) as gtdb:
    for line in gtdb:
        acc, taxonomy, samples = line.rstrip().split('\t')
        if belongs_to_genus(taxonomy):
            biosamples.update(samples.split(','))

with open(output_file, 'w') as file:
    for i, line in oceandna.iterrows():
        if belongs_to_genus(line['gtdb_classification']) and line['original_biosample'] not in biosamples:
            file.write(f"{line['genome']}\t{line['gtdb_classification']}\n")
