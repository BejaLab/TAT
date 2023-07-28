from Bio import Entrez
import xml.etree.ElementTree as ET
import re
from os import path, makedirs
import json

Entrez.email = str(snakemake.config['email'])

tsv_file = snakemake.input['tsv']
txt_file = snakemake.input['txt']
output_file = snakemake.output['tsv']
output_dir = snakemake.output['dir']
genera = snakemake.params['genera']

makedirs(output_dir)

with open(txt_file) as file:
    representatives = [ line.rstrip() for line in file ]

def esearch(acc, db, field):
    with Entrez.esearch(db = db, term = f"{acc}[{field}]") as handle:
        record = Entrez.read(handle)
    assert len(record['IdList']) > 0, f"Nothing found for {acc}"
    return record['IdList'][0]

def elink(id, dbfrom, db):
    with Entrez.elink(id = str(id), dbfrom = dbfrom, db = db) as handle:
        record = Entrez.read(handle)
    links = []
    for linkset in record:
        for link in linkset['LinkSetDb']:
            for linkid in link['Link']:
                links.append(linkid['Id'])
    return links

def esummary(id, db):
    with Entrez.esummary(id = id, db = db) as handle:
        record = Entrez.read(handle)
    if db == "biosample":
        record = record['DocumentSummarySet']['DocumentSummary'][0]
        record['SampleData'] = re.sub(r'\s+', ' ', record['SampleData'].strip())
    if db == "assembly":
        record = record['DocumentSummarySet']['DocumentSummary'][0]
        record['Meta'] = re.sub(r'\s+', ' ', record['Meta'].strip())
    record['id'] = int(id)
    return dict(record)

def get_derived(biosample):
    data = biosample['SampleData']
    data_lower = data.lower()
    if 'derived' in data_lower and 'biosample' in data_lower:
        pattern = r'SAM[A-Z]+\d+'
        samples = []
        for sample in re.findall(pattern, data):
            if sample not in samples and sample != biosample['Accession']:
                samples.append(sample)
        return samples
    return []

with open(tsv_file) as file:
    with open(output_file, 'w') as outfile:
        for line in file:
            gtdb_acc, taxonomy = line.strip().split("\t")
            ingroup = gtdb_acc in representatives
            if not ingroup:
                for genus in genera:
                    if ingroup := f"{genus};" in taxonomy:
                        break
            if ingroup:
                acc = gtdb_acc[3:]
                print(acc)
                assembly_id = esearch(acc, "assembly", "asac")
                assembly = esummary(assembly_id, "assembly")
                samples = {}
                more_biosamples = []
                for biosample_id in elink(assembly_id, "assembly", "biosample"):
                    biosample = esummary(biosample_id, "biosample")
                    samples[biosample['Accession']] = biosample
                    for derived in get_derived(biosample):
                        more_biosamples.append(derived)
                for biosample_acc in more_biosamples:
                    biosample_id = esearch(biosample_acc, "biosample", "accession")
                    biosample = esummary(biosample_id, "biosample")
                    samples[biosample_acc] = biosample
                outfile.write(f"{acc}\t{taxonomy}\t{','.join(samples.keys())}\n")
                json_file = path.join(output_dir, acc + '_assembly.json')
                with open(json_file, 'w') as fh:
                    fh.write(json.dumps(assembly, separators=(',', ':')) + '\n')
                json_file = path.join(output_dir, acc + '_biosamples.json')
                with open(json_file, 'w') as fh:
                    fh.write(json.dumps(samples, separators=(',', ':')) + '\n')
