import ast
import hashlib
import logging
import pandas as pd
import boto3
import smart_open

import pyfaidx
import fsspec
from tqdm import tqdm
import os

from pathlib import Path
from Bio.Seq import Seq

# from conservation import phylocsf
# from seqmap import genomic, utils

from sqlalchemy import and_, or_

from veliadb import base, settings, annotation_loading
from veliadb.base import (Assembly, Gene, Transcript, Protein, Orf, OrfXref,
                          TranscriptOrf, Cds, Dataset, SequenceRegionXref, Exon)
import veliadb.util as vdb_util

pd.options.display.max_columns = 100
pd.options.display.max_rows = 100
pd.options.display.max_colwidth = 200

def extract_nucleotide_sequence_broken_psl_starts(orf, reference):
# orf = [i for i in orfs if ';' in i.block_sizes][2]
    if orf.start == -1 or orf.end == -1:
        return None
    chrom = orf.assembly.ucsc_style_name
    blocks = list(zip(orf.chrom_starts.split(';'), orf.block_sizes.split(';'), orf.phases.split(';')))
    strand = orf.strand
    nucl = []
    for ix, (start, size, phase) in enumerate(blocks):
        start = int(start)
        size = int(size)
        phase = int(phase)
        if strand == '-':
            s = reference[chrom][start+2 : start+size+1]
            s = s.complement.seq[::-1]
        else:
            s = reference[chrom][start-1 : start+size-3]
            s = s.seq
        nucl.append(s)
    if strand == '-' and len(blocks)>1:
        nucl[-1] = reference[chrom][start+3 : start+size].complement.seq[::-1]
    elif strand == '-' and len(blocks)==1:
        nucl[0] = nucl[0][1:]
    if strand == '-':
        nucl = nucl[::-1]
        nucl = ' '.join(nucl)
    else:
        nucl = ' '.join(nucl)
    return nucl, str(Seq.Seq(nucl.replace(' ', '')).translate())