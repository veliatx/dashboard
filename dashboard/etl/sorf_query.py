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
from types import SimpleNamespace

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

import pathlib
current_folder = pathlib.Path(__file__).parents[0]

# class wrapper for orf query results to allow multiprocessing to iterate over orf objects
class OrfData(object):
    def __init__(self, orf_query_object):
        self.id = orf_query_object.id
        self.start = orf_query_object.start
        self.end = orf_query_object.end
        self.strand = orf_query_object.strand
        self.assembly = SimpleNamespace(ucsc_style_name=orf_query_object.assembly.ucsc_style_name)
        self.chrom_starts = orf_query_object.chrom_starts
        self.block_sizes = orf_query_object.block_sizes
        self.phases = orf_query_object.phases
        self.assembly_id = orf_query_object.assembly_id

def extract_nucleotide_sequence_veliadb(orf, reference):
# orf = [i for i in orfs if ';' in i.block_sizes][2]
    if orf.start == -1 or orf.end == -1:
        return '', ''
    chrom = orf.assembly.ucsc_style_name
    blocks = list(zip(orf.chrom_starts.split(';'), orf.block_sizes.split(';'), orf.phases.split(';')))
    strand = orf.strand
    nucl = []
    for ix, (start, size, phase) in enumerate(blocks):
        start = int(start) - 1
        size = int(size)
        phase = int(phase)
        if strand == '-':
            if ix == 0:
                s = reference[chrom][start : start+(size - 3)]
            else:
                s = reference[chrom][start : start+size]
            s = s.complement.seq[::-1]
        else:
            if ix == len(blocks)-1:
                s = reference[chrom][start : start+(size-3)]
            else:
                s = reference[chrom][start : start+size]
            s = s.seq
        nucl.append(s)
    if strand == '-':
        nucl = nucl[::-1]
        nucl = ' '.join(nucl)
    else:
        nucl = ' '.join(nucl)
    return nucl, str(Seq(nucl.replace(' ', '')).translate())

def extract_nucleotide_sequence_broken_psl_starts(orf, reference):
# orf = [i for i in orfs if ';' in i.block_sizes][2]
    if orf.start == -1 or orf.end == -1:
        return '', ''
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
transcripts = pyfaidx.Fasta(os.path.join(current_folder, './gencode.v42.transcripts.fa'))

def find_seq_substring(query, target_dict):
    if query == '':
        return []
    else:
        return [k for k, s in target_dict.items() if query.lower() in s[:].seq.lower()]

if not os.path.exists(os.path.join(current_folder, 'reference.fa')):
    with smart_open.open(GENOME_REFERENCE_PATH) as fhandle:
        with os.path.join(current_folder, 'reference.fa') as f:
            for line in tqdm(fhandle.readlines()):
                f.write(line)
def query_overlapping_transcripts(o, session):
    results = session.query(Transcript).filter(and_(o.start >= Transcript.start, 
                                      o.end <= Transcript.end, 
                                      o.strand == Transcript.strand,
                                      o.assembly_id == Transcript.assembly_id)).all()
    return results

def compute_exact_transcript_matches(o):
    nt, aa = extract_nucleotide_sequence_veliadb(o, reference)
    r_internal = find_seq_substring(nt, transcripts)
    return o.id, nt, aa, [i.split('|')[0].split('.')[0] for i in r_internal]

def run_id_mapping_parallel(orfs, NCPU = 1):
    import multiprocessing as mp
    results = {}
    if NCPU is None:
        NCPU = mp.cpu_count()
    if NCPU > 1:
        with mp.Pool(NCPU) as ppool:            
            for r0 in tqdm(ppool.imap(compute_exact_transcript_matches, orfs), total=len(orfs)):
                results[r0[0]] = r0[1:]
    else:
        for o in tqdm(orfs):
            r = compute_exact_transcript_matches(o)
            results[r[0]] = r[1:]
    return results

import jsonlines
import json
from tqdm import tqdm
def load_jsonlines_table(path_to_file, index_col = None):
    with open(path_to_file) as fh:
        results = []
        for line in tqdm(fh.readlines()):
            results.append(json.loads(line))
    df = pd.DataFrame(results)
    if index_col is not None:
        df.index = df[index_col]
    return df
                
GENOME_REFERENCE_PATH = 's3://velia-annotation-dev/genomes/hg38/GRCh38.p13.genome.fa.gz'
reference = pyfaidx.Fasta(os.path.join(current_folder, 'reference.fa'))

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