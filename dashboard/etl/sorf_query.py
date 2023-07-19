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
import re
from Bio import SeqIO
import pathlib
current_folder = pathlib.Path(__file__).parents[0]

transcripts = SeqIO.parse(os.path.join(current_folder, './gencode.v42.transcripts.fa'), 'fasta')
transcripts = {r.id: str(r.seq).lower() for r in transcripts}

# GENOME_REFERENCE_PATH = 's3://velia-annotation-dev/genomes/hg38/GRCh38.p13.genome.fa.gz'
# reference = SeqIO.parse(os.path.join(current_folder, 'reference.fa'))
# reference = {r.id: str(r.seq) for r in reference}


# from conservation import phylocsf
# from seqmap import genomic, utils

from sqlalchemy import and_, or_

from veliadb import base, settings, annotation_loading
from veliadb.base import (Assembly, Gene, Transcript, Protein, Orf, OrfXref, Protein,
                          TranscriptOrf, Cds, Dataset, SequenceRegionXref, Exon, ProteinXref)
import veliadb.util as vdb_util

pd.options.display.max_columns = 100
pd.options.display.max_rows = 100
pd.options.display.max_colwidth = 200



from copy import deepcopy
def parallel_sorf_query(vtx_id):
# Query DB
    session = base.Session() # connect to db
    orfs = session.query(Orf).filter(Orf.id == vtx_id).all()
    if len(orfs) == 0:
        print(f"{vtx_id} not found in veliadb")
        session.close()
        # return {}
    elif len(orfs) > 1:
        print(f"{vtx_id} had multiple entries found in veliadb")
        session.close()
        # return {}
    else:
        current_orf = orfs[0]
    # Loop over orfs, and populate sorf_table file with attributes of interest
    # nt, aa = extract_nucleotide_sequence_veliadb(current_orf, reference)
    nt = current_orf.nt_seq
    aa = current_orf.aa_seq
    r_internal = find_seq_substring(nt, transcripts)
    transcripts_exact = [i.split('|')[0].split('.')[0] for i in r_internal if i.startswith('ENST')]
    overlapping_tids = query_overlapping_transcripts(current_orf, session)
    overlapping_tids = [[i.split('.')[0] for i in [t.ensembl_id, t.refseq_id, t.chess_id] if i][0] for t in overlapping_tids]
    attributes = parse_orf(current_orf, session)
    attributes['transcripts_exact'] = transcripts_exact
    attributes['transcripts_overlapping'] = overlapping_tids
    # if attributes['aa'] == '':
    attributes['aa'] = aa
    # if attributes['nucl'] == '':
    attributes['nucl'] = nt
    session.close()
    return attributes

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
        self.secondary_orf_id = orf_query_object.secondary_orf_id
        
def parse_orf(orf, session):
    orf_xrefs = ';'.join([ox.xref for ox in session.query(OrfXref).filter(OrfXref.orf_id == orf.id).all() if ox.xref])
    pattern = r'U\w{9}-[0-9]*'
    match = re.search(pattern, orf_xrefs)
    if match:
        genscript_id = match.group()
    else:
        genscript_id = ""
    vtx_id = f'VTX-{orf.id:07d}'
    chrom = orf.assembly.ucsc_style_name
    strand = orf.strand
    start = str(orf.start)
    end = str(orf.end)
    chrom_starts = orf.chrom_starts
    block_sizes = orf.block_sizes
    phases = orf.phases
    seqs = orf.aa_seq
    source = ';'.join(set([x.orf_data_source.name for x in orf.xref]))
    
    
        
    ucsc_track = f'{orf.assembly.ucsc_style_name}:{orf.start}-{orf.end}'

    gene_ids = [g.id for g in session.query(Gene).filter(and_(Gene.assembly_id == orf.assembly_id,
                                                              Gene.strand == orf.strand,
                                                              Gene.start < orf.start,
                                                              Gene.end > orf.end)).all()]
    
    gene_xrefs = ';'.join([sx.xref for sx in session.query(SequenceRegionXref).filter(SequenceRegionXref.sequence_region_id.in_(gene_ids)).all()])
    
    
    transcript_ids = [t.id for t in session.query(Transcript).filter(and_(Transcript.assembly_id == orf.assembly_id,
                                                                          Transcript.strand == orf.strand,
                                                                          Transcript.start < orf.start,
                                                                          Transcript.end > orf.end)).all()]
    
    transcript_xrefs = ';'.join([sx.xref for sx in session.query(SequenceRegionXref).filter(SequenceRegionXref.sequence_region_id.in_(transcript_ids)).all()])
    protein_xrefs = ';'.join([str(px.xref) for px in \
                         session.query(ProteinXref)\
                                .join(Protein, Protein.id == ProteinXref.protein_id)\
                                .filter(Protein.aa_seq == seqs).all()])
    return {
        'vtx_id': vtx_id,
        'genscript_id': genscript_id,
        'velia_id': orf.velia_id,
        'secondary_id': orf.secondary_orf_id,
        'ucsc_track': ucsc_track,
        'chr': chrom,
        'strand': strand,
        'start': start,
        'end': end,
        'chrom_starts': chrom_starts,
        'block_sizes': block_sizes,
        'phases': phases,
        'aa': seqs,
        'nucl': orf.nt_seq,
        'orf_xrefs': orf_xrefs,
        'gene_xrefs': gene_xrefs,
        'transcript_xrefs': transcript_xrefs,
        'protein_xrefs': protein_xrefs,
        'source': source
    }
    
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
# transcripts = pyfaidx.Fasta(os.path.join(current_folder, './gencode.v42.transcripts.fa'))

def find_seq_substring(query, target_dict):
    if query == '':
        return []
    else:
        return [k for k, s in target_dict.items() if query.lower() in s]

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