import re
import json
from pandas import read_csv

lanclos_file = str(snakemake.input)
output_file = str(snakemake.output)

genera = snakemake.params['genera']

lanclos = read_csv(lanclos_file, sep = "\t")

def flatten(dlist):
    return [ id for sublist in dlist for id in sublist ]

with open(output_file, 'w') as file:
    for i, line in lanclos.iterrows():
        source = line['GTDB_accession']
        subclade = line['SubcladeA']
        if (source == 'This study' or source == 'SFB') and subclade in genera:
            file.write(f"{line['Name']}\t{subclade}\n")
