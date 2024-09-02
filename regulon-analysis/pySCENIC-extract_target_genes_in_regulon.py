# Checking the target gene list from pySCENIC regulon (ctx output)

import pandas as pd
import os
from pyscenic.transform import df2regulons as df2regs
import json

# use function from pyscenic for parsing the 'pyscenic ctx' output
def _df2regulons(fname):
  ext = os.path.splitext(fname,)[1]
  df = pd.read_csv(fname, sep=',', index_col=[0,1], header=[0,1], skipinitialspace=True)
  df[('Enrichment', 'Context')] = df[('Enrichment', 'Context')].apply(lambda s: eval(s))
  df[('Enrichment', 'TargetGenes')] = df[('Enrichment', 'TargetGenes')].apply(lambda s: eval(s))
  return df2regs(df)

# Form dictionary of { TF : Target } pairs from 'pyscenic ctx' output.
os.chdir("pySCENIC/output")
rdict = {}
regulons = _df2regulons('reg.csv')
for reg in regulons:
  targets = [ target for target in reg.gene2weight ]
  rdict[reg.name] = targets

# Write to JSON for import to R
rjson = json.dumps(rdict)
f = open("regulons.json","w")
f.write(rjson)
f.close()

library(jsonlite)
regs <- read_json('regulons.json', simplify = T)
