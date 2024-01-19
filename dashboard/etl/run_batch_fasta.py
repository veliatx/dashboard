import pathlib
import torch
import esm
from esm.esmfold.v1.esmfold import ESMFold
import os
from Bio import SeqIO
import sys
from tqdm import tqdm
INPUT = sys.argv[1]
OUTPUT = sys.argv[2]
print(OUTPUT)

sequences = SeqIO.parse(INPUT, 'fasta')
sequences = {k.id: str(k.seq) for k in sequences}

model = esm.pretrained.esmfold_v1()
model = model.eval().cuda()

import json
import numpy as np
# Multimer prediction can be done with chains separated by ':'


with torch.no_grad():
    with open(OUTPUT, 'w') as fout:

        for name, seq in tqdm(sequences.items()):
            seq = str(seq)
            try:
                output = model.infer_pdb(seq)
                plddt = model.infer(seq)
                fout.write(json.dumps({'id': name,
                                    'sequence': seq,
                                    'pdb': output,
                                    'mean_plddt': float(plddt['mean_plddt'].detach().cpu().numpy()[0]),
                                    'ptm': float(plddt['ptm'].detach().cpu().numpy()[0]),
                                    'plddt': [float(i) for i in plddt['plddt'].detach().cpu().numpy().mean(axis=(0, 2))],
                                    })
                        )
                fout.write('\n')
            except Exception as e:
                fout.write(json.dumps({'id': name,
                                    'sequence': seq,
                                    'pdb': 'failed',
                                    'mean_plddt': -1,#float(plddt['mean_plddt'].detach().cpu().numpy()[0]),
                                    'ptm': -1, #float(plddt['ptm'].detach().cpu().numpy()[0]),
                                    'plddt': [-1]#[float(i) for i in plddt['plddt'].detach().cpu().numpy().mean(axis=(0, 2))],
                                    }))
                print(name, seq, 'failed to infer pdb')
