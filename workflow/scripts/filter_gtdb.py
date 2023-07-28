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
genera_gtdb = genera.values()

makedirs(output_dir)

with open(txt_file) as file:
    representatives = [ line.rstrip() for line in file ]

def esearch(acc, db, field, check = True):
    with Entrez.esearch(db = db, term = f"{acc}[{field}]") as handle:
        record = Entrez.read(handle)
    assert len(record['IdList']) < 2, f"Multiple matches for {acc}: {record['IdList']}"
    if len(record['IdList']) > 0:
        return record['IdList'][0]
    else:
        return None

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
        record['Attributes'] = biosample_attributes(record)
    if db == "assembly":
        record = record['DocumentSummarySet']['DocumentSummary'][0]
        record['Meta'] = re.sub(r'\s+', ' ', record['Meta'].strip())
    if db == "sra":
        record = record[0]
    record['id'] = int(id)
    return dict(record)

def biosample_attributes(summary):
    attrs = {}
    root = ET.fromstring(summary['SampleData'])
    for child in root:
        if child.tag == 'Attributes':
            for attr in child:
                if attr.text:
                    attrs[attr.get("attribute_name")] = attr.text
    return attrs

def runs_from_biosample(biosample):
    pattern = r'[A-Z]RR\d+'
    runs = set()
    for run in re.findall(pattern, biosample['SampleData']):
        runs.add(run)
    return list(runs)

def biosamples_from_biosample(biosample):
    pattern = r'SAM[A-Z]+\d+'
    samples = set()
    for sample in re.findall(pattern, biosample['SampleData']):
        if sample != biosample['Accession']:
            samples.add(sample)
    return list(samples)

def sras_from_biosample(biosample):
    pattern = r'[A-Z]RS\d+'
    sras = set()
    for sra in re.findall(pattern, biosample['SampleData']):
        sras.add(sra)
    return list(sras)

#deprecated
def runs_from_biosample_attrs(attrs):
    sra_note = ''
    for key, val in summary['Attributes'].items():
        if 'sequence read archive' in val.lower() or 'SRA' in val:
            sra_note += ' ' + val
    pattern = r'[A-Z]RR\d+'
    return re.findall(pattern, sra_note)

def runs_from_sra(summary):
    runs = []
    root = ET.fromstring(f"<root>{summary['Runs']}</root>")
    for child in root:
        runs.append(child.get("acc"))
    return runs

def get_derived_attrs(attrs):
    if 'derived_from' in attrs and 'biosample' in attrs['derived_from'].lower():
        pattern = r'SAM[A-Z]\d+'
        return re.findall(pattern, attrs['derived_from'])[0]
    return None

with open(tsv_file) as file:
    with open(output_file, 'w') as outfile:
        for line in file:
            gtdb_acc, taxonomy = line.strip().split("\t")
            ingroup = gtdb_acc in representatives
            if not ingroup:
                for genus in genera_gtdb:
                    if ingroup := f"{genus};" in taxonomy:
                        break
            if ingroup:
                acc = gtdb_acc[3:]
                print(acc)
                assembly_id = esearch(acc, "assembly", "asac")
                assembly = esummary(assembly_id, "assembly")
                biosamples = []
                more_biosamples = []
                for biosample_id in elink(assembly_id, "assembly", "biosample"):
                    biosample = esummary(biosample_id, "biosample")
                    biosamples.append(biosample)
                    more_biosamples += biosamples_from_biosample(biosample)
                for biosample_acc in more_biosamples:
                    biosample_id = esearch(biosample_acc, "biosample", "accession")
                    biosample = esummary(biosample_id, "biosample")
                    biosamples.append(biosample)
                data = []
                for biosample in biosamples:
                    #runs = runs_from_biosample(biosample)
                    #sra_ids = elink(biosample['id'], "biosample", "sra")
                    #for sra_acc in sras_from_biosample(biosample):
                    #    sra_id = esearch(sra_acc, "sra", "Accession")
                    #    if sra_id:
                    #        sra_ids.append(sra_id)
                    #for sra_id in sra_ids:
                    #    sra = esummary(sra_id, "sra")
                    #    runs += runs_from_sra(sra)
                    #biosample['Runs'] = runs
                    data.append(biosample['Accession'])
                assert len(data) > 0, f"No data for {acc}"
                outfile.write(f"{acc}\t{taxonomy}\t{','.join(data)}\n")
                with open(path.join(output_dir, acc + '_assembly.json'), 'w') as fh:
                    fh.write(json.dumps(assembly, separators=(',', ':')) + '\n')
                with open(path.join(output_dir, acc + '_biosamples.json'), 'w') as fh:
                    fh.write(json.dumps(biosamples, separators=(',', ':')) + '\n')
