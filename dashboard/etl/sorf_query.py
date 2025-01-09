"""Module for querying and processing sORF (small Open Reading Frame) data.

This module provides functions and utilities for:
- Loading and parsing sORF sequences from databases
- Querying sORF metadata and annotations
- Processing transcript and gene mappings
- Handling phase/screening information

The module interfaces with VeliaDB and external sequence databases to extract
and process sORF-related data.
"""

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
                         TranscriptOrf, Cds, Dataset, SequenceRegionXref, Exon, ProteinXref,
                         TranscriptXref)
from dashboard.etl import CACHE_DIR, DATA_DIR

pd.options.display.max_columns = 100
pd.options.display.max_rows = 100
pd.options.display.max_colwidth = 200

current_folder = pathlib.Path(__file__).parents[0]


def parse_fasta(fasta_content: str) -> list:
    """Parse FASTA format content into SeqIO records.

    Args:
        fasta_content: String containing FASTA formatted sequences

    Returns:
        List of SeqIO records
    """
    fasta_io = StringIO(fasta_content)
    records = list(SeqIO.parse(fasta_io, "fasta"))
    return records


# Load transcript sequences
with smart_open.open('s3://velia-annotation-dev/gencode/v42/gencode.v42.transcripts.fa.gz') as f:
    transcripts = f.read()
    transcripts = parse_fasta(transcripts)
    transcripts = {r.id: str(r.seq).lower() for r in transcripts}


def parallel_sorf_query(vtx_id: str) -> dict:
    """Query sORF data for a given VTX ID.

    Args:
        vtx_id: VTX identifier for the sORF

    Returns:
        Dictionary containing sORF attributes and annotations
    """
    session = base.Session()  # connect to db
    try:
        orf = session.query(Orf).filter(Orf.vtx_id == vtx_id).one()
    except Exception as e:
        print(e)

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

    attributes = parse_orf(orf, session)
    attributes['aa'] = aa.replace('*', '')
    attributes['nucl'] = nt
    session.close()
    return attributes


class OrfData:
    """Class wrapper for ORF query results to allow multiprocessing.

    Attributes:
        id: ORF identifier
        start: Start position
        end: End position
        strand: Strand orientation
        assembly: Assembly information
        chrom_starts: Chromosome start positions
        block_sizes: Exon block sizes
        phases: Phase information
        assembly_id: Assembly identifier
        secondary_orf_id: Secondary ORF identifier
    """

    def __init__(self, orf_query_object):
        """Initialize OrfData object from database query result."""
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


def parse_orf(orf: Orf, session) -> dict:
    """Parse ORF object and extract relevant attributes.

    Args:
        orf: ORF database object
        session: Database session

    Returns:
        Dictionary containing parsed ORF attributes
    """
    orf_xrefs = ';'.join([ox.xref for ox in session.query(OrfXref).filter(OrfXref.orf_id == orf.id).all() if ox.xref])
    pattern = r'U\w{9}-[0-9]*'
    match = re.search(pattern, orf_xrefs)
    genscript_id = match.group() if match else ""

    vtx_id = orf.vtx_id
    chrom = orf.assembly.ucsc_style_name
    strand = orf.strand
    start = str(orf.start)
    end = str(orf.end)
    chrom_starts = orf.chrom_starts
    block_sizes = orf.block_sizes
    phases = orf.phases
    seqs = orf.aa_seq.rstrip('*')

    source = ';'.join(set([x.orf_data_source.name for x in orf.xref]))
    ucsc_track = f'{orf.assembly.ucsc_style_name}:{orf.start}-{orf.end}'

    transcript_ids = list(set([t.transcript_id for t in session.query(TranscriptOrf).filter(TranscriptOrf.orf_id == orf.id).all()]))
    transcript_xrefs = ';'.join([sx.xref for sx in session.query(TranscriptXref).filter(TranscriptXref.transcript_id.in_(transcript_ids)).all()])
    transcript_objects = [t for t in session.query(Transcript).filter(Transcript.id.in_(transcript_ids)).all()]
    transcripts_exact = set([get_first_non_empty_property(t) for t in transcript_objects])
    gene_ids = set([t.gene.id for t in transcript_objects])

    protein_xrefs = ';'.join([str(px.xref) for px in
                             session.query(ProteinXref)
                             .join(Protein, Protein.id == ProteinXref.protein_id)
                             .filter(Protein.aa_seq == seqs).all()])

    gene_xrefs = ';'.join([sx.xref for sx in session.query(SequenceRegionXref).filter(SequenceRegionXref.sequence_region_id.in_(gene_ids)).all()])
    exact_transcripts = [t for t in session.query(Transcript).filter(Transcript.id.in_(transcript_ids)).all() if t.ensembl_id != '']
    exact_transcript_ids = [t.ensembl_id.split('.')[0] for t in exact_transcripts]

    return {
        'vtx_id': vtx_id,
        'genscript_id': genscript_id,
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
        'transcripts_exact': exact_transcript_ids,
        'protein_xrefs': protein_xrefs,
        'source': source
    }


def get_first_non_empty_property(obj: Transcript) -> str:
    """Get first non-empty identifier from transcript object.

    Args:
        obj: Transcript object

    Returns:
        First available identifier (ensembl, refseq, chess or internal id)
    """
    if obj.ensembl_id != '':
        id = obj.ensembl_id
    elif obj.refseq_id != '':
        id = obj.refseq_id
    elif obj.chess_id != '':
        id = obj.chess_id
    else:
        id = str(obj.id)
    return id


def load_jsonlines_table(path_to_file: str, index_col: str = None) -> pd.DataFrame:
    """Load JSON Lines file into pandas DataFrame.

    Args:
        path_to_file: Path to JSON Lines file
        index_col: Column to use as index

    Returns:
        DataFrame containing parsed JSON Lines data
    """
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


def parse_sorf_phase(sorf_df: pd.DataFrame, session) -> tuple:
    """Parse screening phase information from sORF data.

    Args:
        sorf_df: DataFrame containing sORF data
        session: Database session

    Returns:
        Tuple of (phase_ids, phase_entries)
    """
    phase_ids = []
    phase_entries = []
    for row in sorf_df.itertuples():
        if row.vtx_id.startswith('Phase'):
            phase_ids.append(row.vtx_id)
            phase_entries.append(f'Phase {row.vtx_id[6]}')
        elif row.secondary_id.startswith('Phase'):
            phase_ids.append(row.secondary_id)
            phase_entries.append(f'Phase {row.secondary_id[6]}')
        elif any(['Phase' in xref for xref in row.orf_xrefs]):
            phase_id = [x for x in row.orf_xrefs if x.startswith('Phase')][0]
            phase_ids.append(phase_id)
            phase_entries.append(f'Phase {phase_id[6]}')
        elif row.secondary_id.startswith('smORF'):
            phase_ids.append(row.vtx_id)
            phase_entries.append('Phase 1')
        else:
            orf = session.query(Orf).filter(Orf.vtx_id == row.vtx_id).one()

            if orf.secondary_orf_id.startswith('Phase'):
                phase_entries.append(f'Phase {orf.vtx_id[6]}')
                phase_ids.append(orf.vtx_id)
            else:
                phase_entries.append('-1')
                phase_ids.append('-1')
    return phase_ids, phase_entries


def fix_missing_phase_ids(sorf_df: pd.DataFrame) -> pd.DataFrame:
    """Fix missing phase IDs in sORF DataFrame.

    Args:
        sorf_df: DataFrame containing sORF data

    Returns:
        DataFrame with fixed phase IDs
    """
    sorf_df = sorf_df.copy()
    missing_phase_ids = sorf_df[sorf_df['screening_phase'] == '-1'].index
    interim_sheet = pd.read_csv(os.path.join(DATA_DIR, 'interim_phase1to7_all_20230717.csv'), index_col=0)

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
