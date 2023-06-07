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
    return nucl, str(Seq(nucl.replace(' ', '')).translate())

import pyfaidx
import fsspec
from tqdm import tqdm
import os

# deduped_conservation = pd.read_parquet('deduped_phylocsf_length_filtered.parq')
transcripts = pyfaidx.Fasta('./gencode.v42.transcripts.fa')

def find_seq_substring(query, target_dict):
    return [k for k, s in target_dict.items() if query.lower() in s[:].seq.lower()]

if not os.path.exists('reference.fa'):
    with smart_open.open(GENOME_REFERENCE_PATH) as fhandle:
        with open('reference.fa', 'w') as f:
            for line in tqdm(fhandle.readlines()):
                f.write(line)
def query_overlapping_transcripts(o, session):
    results = session.query(Transcript).filter(and_(o.start >= Transcript.start, 
                                      o.end <= Transcript.end, 
                                      o.strand == Transcript.strand,
                                      o.assembly_id == Transcript.assembly_id)).all()
    return results

def compute_exact_transcript_matches(o):
    nt, aa = extract_nucleotide_sequence_broken_psl_starts(current_orf, reference)
    r = find_seq_substring(nt, transcripts)
    return o.id, nt, aa, [i.split('|')[0].split('.')[0] for i in r]

def run_id_mapping_parallel(orfs):
    import multiprocessing as mp
    NCPU = mp.cpu_count()
    with mp.Pool(NCPU) as ppool:
        results = {}
        for r in tqdm(ppool.imap(compute_exact_transcript_matches, orfs), total=len(orfs)):
            results[r[0]] = r[1:]
    return results
                
GENOME_REFERENCE_PATH = 's3://velia-annotation-dev/genomes/hg38/GRCh38.p13.genome.fa.gz'
reference = pyfaidx.Fasta('reference.fa')

# with jsonlines.open('sorf_table.jsonlines', mode = 'w') as fh:
#     for current_orf in tqdm(orfs):
#         nt, aa = extract_nucleotide_sequence_broken_psl_starts(current_orf, reference)
#         tids = find_seq_substring(nt, transcripts)
#         # if len(tids)>0:
#         tids = [i.split('|')[0].split('.')[0] for i in tids]
#         attributes = {
#             'chr': current_orf.assembly.ucsc_style_name,
#             'vtx': f"VTX-{str(current_orf.id).zfill(7)}",
#             'strand': current_orf.strand,
#             'start': current_orf.start,
#             'end': current_orf.end,
#             'nucl': nt,
#             'aa': aa,
#             'transcripts': tids
#         }
#         fh.write(attributes)