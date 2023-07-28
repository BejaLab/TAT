import urllib.request
from os.path import basename
import gzip
import shutil
import json

json_file = snakemake.params['json']
output_file = str(snakemake.output)

source = 'GenBank' if snakemake.wildcards['acc'].startswith('GCA') else 'RefSeq'

with open(json_file) as file:
    assembly = json.load(file)
    url = assembly['FtpPath_' + source]

file = basename(url) + '_genomic.fna.gz'
response = urllib.request.urlopen(f"{url}/{file}".replace("ftp://", "https://"))
with gzip.GzipFile(fileobj = response) as gzipped_file:
    with open(output_file, 'wb') as decompressed_file:
        shutil.copyfileobj(gzipped_file, decompressed_file)
