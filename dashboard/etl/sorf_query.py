import json
import jsonlines
import os
import pathlib
import re
import smart_open

import pandas as pd
from tqdm import tqdm

from Bio.Seq import Seq
from Bio import SeqIO
from types import SimpleNamespace
from io import StringIO
from copy import deepcopy
from sqlalchemy import and_, or_

from veliadb import base
from veliadb.base import (Assembly, Gene, Transcript, Protein, Orf, OrfXref, Protein,
                          TranscriptOrf, Cds, Dataset, SequenceRegionXref, Exon, ProteinXref)
from dashboard.etl import CACHE_DIR, DATA_DIR

pd.options.display.max_columns = 100
pd.options.display.max_rows = 100
pd.options.display.max_colwidth = 200
    
def extract_nucleotide_sequence_veliadb(o, genome_seq):
    seqs = []
    blocks = [int(x) for x in o.block_sizes.split(';')]
    chrom_starts = [int(x)-1 for x in o.chrom_starts.split(';')]

    if '_' in o.assembly.ucsc_style_name:
        chrom = o.assembly.ucsc_style_name.split('_')[1].replace('v', '.')
    else:
        chrom = o.assembly.ucsc_style_name
   
    for i, block in enumerate(blocks):
        if o.strand == '-':
            seqs.append(str(genome_seq[chrom][chrom_starts[i]:(chrom_starts[i] + block + 1)].reverse_complement().seq).upper())
        else:
            seqs.append(str(genome_seq[chrom][chrom_starts[i]:(chrom_starts[i] + block + 1)].seq).upper())
    
    if o.strand == '-':    
        seqs.reverse()
    seq = ''.join(seqs)     
    return seq, str(Seq(seq).translate())


current_folder = pathlib.Path(__file__).parents[0]
# Function to parse the FASTA content
def parse_fasta(fasta_content):
    fasta_io = StringIO(fasta_content)
    records = list(SeqIO.parse(fasta_io, "fasta"))
    return records


with smart_open.open('s3://velia-annotation-dev/gencode/v42/gencode.v42.transcripts.fa.gz') as f:
    transcripts = f.read()
    transcripts = parse_fasta(transcripts)
    transcripts = {r.id: str(r.seq).lower() for r in transcripts}


def parallel_sorf_query(vtx_id, genome_reference):
    # Query DB
    session = base.Session() # connect to db
    try:
        orf = session.query(Orf).filter(Orf.id == vtx_id).one()
    except Exception as e:
        print(e)
    
    if orf.start != -1 and orf.nt_seq == '' or orf.aa_seq == '':
        try:
            nt_reconstructed, aa_reconstructed = extract_nucleotide_sequence_veliadb(orf, genome_reference)
            aa_reconstructed = aa_reconstructed.strip('*')
        except:
            print(f'Could not extract sequence for {vtx_id}')
            if orf.aa_seq == '':
                print(f'No registered AA seq')
                return
    else:
        nt_reconstructed = orf.nt_seq
        aa_reconstructed = orf.aa_seq
        
    if orf.nt_seq == nt_reconstructed:
        nt = orf.nt_seq
    elif orf.nt_seq == '':
        nt = nt_reconstructed
    else:
        nt = orf.nt_seq

    if orf.aa_seq == aa_reconstructed:
        aa = orf.aa_seq
    elif orf.aa_seq == '':
        aa = aa_reconstructed
    else:
        print(f'warning: {orf.id} has inconsistent coordinates in veliadb')
        aa = orf.aa_seq

    aa = aa.strip('*')

    r_internal = find_seq_substring(nt, transcripts)
    transcripts_exact = [i.split('|')[0] for i in r_internal if i.startswith('ENST')]
    # overlapping_tids = query_overlapping_transcripts(current_orf, session)
    # overlapping_tids = [[i.split('.')[0] for i in [t.ensembl_id, t.refseq_id, t.chess_id] if i][0] for t in overlapping_tids]
    attributes = parse_orf(orf, session)
    attributes['transcripts_exact'] = transcripts_exact
    # attributes['transcripts_overlapping'] = overlapping_tids
    # if attributes['aa'] == '':
    attributes['aa'] = aa.replace('*', '')
    # if attributes['nucl'] == '':
    attributes['nucl'] = nt
    session.close()
    return attributes


class OrfData(object):
    """
    Class wrapper for orf query results to allow multiprocessing to iterate over orf objects
    """
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


def find_seq_substring(query, target_dict):
    if query == '':
        return []
    else:
        return [k for k, s in target_dict.items() if query.lower() in s]

#if not os.path.exists(os.path.join(current_folder, 'reference.fa')):
#    with smart_open.open(GENOME_REFERENCE_PATH) as fhandle:
#        with open(os.path.join(current_folder, 'reference.fa'), 'w') as f:
#            for line in tqdm(fhandle.readlines()):
#                f.write(line)
# def query_overlapping_transcripts(o, session):
#     results = session.query(Transcript).filter(and_(o.start >= Transcript.start, 
#                                       o.end <= Transcript.end, 
#                                       o.strand == Transcript.strand,
#                                       o.assembly_id == Transcript.assembly_id)).all()
#     return results

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


def load_jsonlines_table(path_to_file, index_col = None):
    with open(path_to_file) as fh:
        results = []
        for line in tqdm(fh.readlines()):
            res = json.loads(line)
            if res:
                results.append(res)
    df = pd.DataFrame(results)
    if index_col is not None:
        df.index = df[index_col]
    return df


def parse_sorf_phase(sorf_df, session):
    phase_ids = []
    phase_entries = []
    for row in sorf_df.itertuples():
        if row.velia_id.startswith('Phase'):
            phase_ids.append(row.velia_id)
            phase_entries.append(f'Phase {row.velia_id[6]}')
        elif row.secondary_id.startswith('Phase'):
            phase_ids.append(row.secondary_id)
            phase_entries.append(f'Phase {row.secondary_id[6]}')
        elif any(['Phase' in xref for xref in row.orf_xrefs]):
            phase_id = [x for x in row.orf_xrefs if x.startswith('Phase')][0]
            phase_ids.append(phase_id)
            phase_entries.append(f'Phase {phase_id[6]}')
        elif row.velia_id.startswith('smORF'):
            phase_ids.append(row.velia_id)
            phase_entries.append('Phase 1')
        else:
            orf = session.query(Orf).filter(Orf.id == int(row.vtx_id.split('-')[1])).one()
            
            if orf.velia_id.startswith('Phase'):
                phase_entries.append(f'Phase {orf.velia_id[6]}')
                phase_ids.append(orf.velia_id)
            else:
                phase_entries.append('-1')
                phase_ids.append('-1')
    return phase_ids, phase_entries

def fix_missing_phase_ids(sorf_df):
    sorf_df = sorf_df.copy()
    missing_phase_ids = sorf_df[sorf_df['screening_phase']=='-1'].index
    interim_sheet = pd.read_csv(os.path.join(DATA_DIR, 'interim_phase1to7_all_20230717.csv'), index_col=0)
    # source = interim_sheet.loc[sorf_df[sorf_df['screening_phase']=='-1'].index]['source']
    new_phase_ids = []
    for vtx in missing_phase_ids:
        if vtx in interim_sheet.index:
            s = interim_sheet.loc[vtx, 'source']
        else:
            s = 'Not Screened'
        if s.startswith('velia_phase1'):
            i = 'Phase 1'
        elif s.startswith('velia_phase2'):
            i = 'Phase 2'
        elif s.startswith('velia_phase3'):
            i = 'Phase 3'
        else:
            i = s
        new_phase_ids.append(i)
    sorf_df.loc[missing_phase_ids, 'screening_phase'] = new_phase_ids
    return sorf_df
