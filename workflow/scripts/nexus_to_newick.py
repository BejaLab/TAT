from dendropy import TreeList
import re

tree = TreeList.get(
    path = str(snakemake.input),
    schema = "nexus",
    preserve_underscores = True
).as_string(
    schema = "newick"
)

tree = re.sub('["\']|\[&R\] |:[0-9.]+', '', tree)

with open(str(snakemake.output), 'w') as file:
    file.write(tree)
