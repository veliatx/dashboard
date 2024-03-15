"""Script to update the dashboard cache"""
from dashboard.etl.sorf_query import parallel_sorf_query, load_jsonlines_table, fix_missing_phase_ids, parse_sorf_phase
from dashboard.etl.transcript_features import load_xena_transcripts_with_metadata_from_s3, process_sums_dataframe_to_heatmap, create_comparison_groups_xena_tcga_vs_normal, read_tcga_de_from_s3, load_de_results
from dashboard.data_load import load_esmfold, load_mouse_blastp_results
from dashboard.etl.etl_utils import fasta_write_veliadb_protein_sequences, merge_sorf_df_blast, write_nonsignal_aa_sequences, parse_blastp_json
from dashboard.tabs.riboseq_atlas import get_average_coverage
from dashboard.etl import module_path
from dashboard.etl.run_protein_search_tools import run_protein_search_tools

from veliadb import base
from veliadb.base import Protein, ProteinXref, Orf

from Bio import SeqIO
from io import BytesIO
from pathlib import Path
from sqlalchemy import create_engine, text
from tqdm import tqdm

import multiprocessing as mp
import numpy as np
import pandas as pd

import boto3
import functools
import gzip
import click
import json
import logging
import pickle
import time
import shlex
import shutil
import subprocess
import sqlite3
import sys

def configure_logger(log_file=None, level=logging.INFO, overwrite_log=True,
                     format=logging.BASIC_FORMAT):
    """
    Helper to configure logging output.

    log_file: str
        Name of logfile
    level: str
        Density of log output
    overwrite_log: bool
        Flag to overwrite existing log
    format: str
        Type of logging output
    """

    if log_file is None:
        logging.basicConfig(stream=sys.stdout, level=level, format=format)
    else:
        logging.basicConfig(filename=log_file, level=level,
                            filemode=('w' if overwrite_log else 'a'),
                            format=format)
        console = logging.StreamHandler()
        console.setLevel(level)
        console.setFormatter(logging.Formatter(format))
        logging.getLogger('').addHandler(console)


def get_genome_reference(bucket_name, s3_file_key):
    """
    Helper function to retrieve genome reference

    bucket_name: str
        Valid Velia S3 bucket name
    """    
    s3_client = boto3.client('s3')

    file_obj = s3_client.get_object(Bucket=bucket_name, Key=s3_file_key)
    fasta_gzipped = file_obj['Body'].read()
    with gzip.open(BytesIO(fasta_gzipped), 'rt') as f:
        genome_reference = SeqIO.to_dict(SeqIO.parse(f, "fasta"))

    return genome_reference


def helper_function(arg, extra_arg):
    return pfunc(arg, extra_arg)


def pfunc(i, genome_reference):
        r = parallel_sorf_query(i, genome_reference)
        return r

@click.command()
@click.argument('vtx_ids_file', type=click.Path(exists=True, resolve_path=True))
@click.argument('cache_dir', type=click.Path(resolve_path=True))
@click.argument('data_dir', type=click.Path(exists=True, resolve_path=True))
# @click.argument('protein_tools_path', type=click.Path(exists=True, resolve_path=True))

@click.option("--overwrite", is_flag=True, show_default=True, 
              default=False, help="Whether to overwrite existing cache directory")
@click.option("--resume", is_flag=True, show_default=True, 
              default=False, help="Whether to overwrite existing cache directory")
@click.option("--run_protein_tools", is_flag=True, show_default=True, 
              default=False, help="Whether to overwrite protein_data in the cache directory")
@click.option("--number_threads", help="Number of threads to use for parallel computing.", default=8, type=int)

def update_cache(vtx_ids_file, cache_dir, data_dir, overwrite, resume, run_protein_tools, number_threads):
    """
    Primary function to update cache

    vtx_ids_file: Path to file containing a single VTX entry per line

    cache_dir: Path to desired output directory
    
    data_dir: Path to directory with data folder for dashboard
    
    protein_tools_path: Path to directory with protein_tools repo
    
    overwrite: Flag indicating whether to overwrite ALL existing results
    
    resume: Flag indicating whether to resume partial run on existing cache dir
    
    run_protein_tools: Flag indicating whether to run protein_tools (takes 12hrs+)
    """

    configure_logger(f"{time.strftime('%Y%m%d_%H%M%S')}_veliadb_load_db.log",
                     level=logging.INFO)

    cache_dir = Path(cache_dir).absolute()

    bucket_name = 'velia-annotation-dev'
    s3_file_key = 'genomes/hg38/GRCh38.p13.genome.fa.gz'
    logging.info('Loading hg38 genome')
    genome_reference = get_genome_reference(bucket_name, s3_file_key)

    if cache_dir.exists() and not overwrite and not resume:
        logging.info(f'Cache directory {cache_dir} exists and overwrite and resume are set to false')
        return
    elif cache_dir.exists() and resume:
        logging.info(f'Running job from existing cache {cache_dir}')
    elif cache_dir.exists() and overwrite:
        shutil.rmtree(cache_dir)
        cache_dir.mkdir()
        cache_dir.joinpath('protein_data').mkdir()
    elif not cache_dir.exists():
        cache_dir.mkdir()
        cache_dir.joinpath('protein_data').mkdir()

    with open(vtx_ids_file) as fhandle:
        ids = [int(i.replace('VTX-', '')) for i in fhandle.readlines()]
        ids = list(set(ids))

    sorfs_json_exists = cache_dir.joinpath('sorf_table.jsonlines').exists()
    if overwrite or not sorfs_json_exists:
        # Query DB
        logging.info("Query Orf Table")
        session = base.Session() 
        orfs = session.query(Orf).filter(Orf.id.in_(ids)).all()
        missing_orfs = set(ids) - set([i.id for i in orfs])
        
        if len(missing_orfs) > 0:
            logging.info('WARNING: some of the provided IDs were not found in veliadb.', *missing_orfs)
        with open(cache_dir.joinpath('sorf_table.jsonlines'), 'w') as fopen:
            with mp.Pool(number_threads) as ppool:
                func = functools.partial(helper_function, extra_arg=genome_reference)
                for r in tqdm(ppool.imap(func, ids, chunksize=number_threads), total=len(ids)):
                    fopen.write(json.dumps(r))
                    fopen.write('\n')

        session.close()

    sorf_df = load_jsonlines_table(cache_dir.joinpath('sorf_table.jsonlines'), index_col='vtx_id')
    
    del genome_reference

    # Format table to conform to standardized entries
    session = base.Session()
    sorf_df['index_copy'] = sorf_df.index
    sorf_df['show_details'] = False
    sorf_df['orf_xrefs'] = sorf_df.apply(lambda x: tuple(x.orf_xrefs.split(';')), axis=1)
    sorf_df['source'] = sorf_df.apply(lambda x: tuple(x.source.split(';')), axis=1)

    cols = list(sorf_df.columns)
    cols.insert(0, cols.pop(cols.index('show_details')))

    phase_ids = []
    phase_entries = []
    protein_xrefs = []
    for row in tqdm(sorf_df.itertuples()):
        protein_xrefs.append(tuple([str(px.xref) for px in \
                                session.query(ProteinXref)\
                                        .join(Protein, Protein.id == ProteinXref.protein_id)\
                                        .filter(Protein.aa_seq == row.aa).all()]))

    logging.info('parsing sorf_phase')
    phase_ids, phase_entries = parse_sorf_phase(sorf_df, session)
    sorf_df['screening_phase_id'] = phase_ids
    sorf_df['screening_phase'] = phase_entries
    sorf_df['protein_xrefs'] = protein_xrefs
    sorf_df['aa_length'] = sorf_df.apply(lambda x: len(x.aa), axis=1)

    sorf_df = fix_missing_phase_ids(sorf_df)

    cols = ['show_details', 'vtx_id', 'aa_length', 'screening_phase_id', 'screening_phase', 'ucsc_track', 
            'source', 'orf_xrefs', 'protein_xrefs', 'gene_xrefs', 'transcript_xrefs',  
            'transcripts_exact', 'aa', 'nucl', 
            'index_copy', 'genscript_id', 'chr', 'strand', 'start', 'end', 
            'chrom_starts', 'block_sizes', 'phases',]

    sorf_df = sorf_df[cols]
    session.close()
    # Add ribo-seq coverage
    ribo_df = get_average_coverage()
    vtx_with_any_support = ribo_df[(ribo_df.sum(axis=1)>50) & (ribo_df.max(axis=1)>10)].index
    array_to_add = ['True' if i in vtx_with_any_support else 'False' for i in sorf_df.index]
    sorf_df['Ribo-Seq RPKM Support'] = array_to_add
    
    # Generate fasta file for protein_tools run
    with open(cache_dir.joinpath('protein_data', 'protein_tools_input.fasta'), 'w') as fopen:
       for ix, row in sorf_df.iterrows():
           fopen.write(f">{row['vtx_id']}\n{row['aa'].replace('*', '')}\n")

    if run_protein_tools:
        vtx_aa_fasta_path = cache_dir.joinpath('protein_data', 'protein_tools_input.fasta')
        output_prefix_path = cache_dir.joinpath('protein_data')
        protein_etl_script = __file__.replace('update_cache.py', 'dashboard_etl.py')
        cmd = f"python {protein_etl_script} -i {vtx_aa_fasta_path} -o {output_prefix_path}"
        logging.info(f'Running protein feature prediction tools {cmd}')
        subprocess.run(shlex.split(cmd))
        write_nonsignal_aa_sequences(output_prefix_path)
        fasta_write_veliadb_protein_sequences(output_prefix_path) # Writes fastas with reference protein sequences in cache
        logging.info(f'Running protein search tools {cmd}')
        cmd = f"run_protein_search_tools {vtx_aa_fasta_path} {output_prefix_path} --number_threads {number_threads}"
        subprocess.run(shlex.split(cmd))
        hits_per_query, blast_sorf_table = parse_blastp_json(cache_dir.joinpath('protein_data', 'blastp.results.json'))
        with open(cache_dir.joinpath('protein_data', 'blastp.results.pkl'), 'wb') as fwrite:
            pickle.dump((hits_per_query, blast_sorf_table), fwrite)

        logging.info("Finished running protein tools.")
        
    # Massage table to standard format
    _, blastp_table = load_mouse_blastp_results(cache_dir)
    sorf_df = merge_sorf_df_blast(sorf_df, blastp_table, cache_dir.joinpath('protein_data'))

    
    logging.info(f'Adding ESMFold data')
    esmfold = load_esmfold(cache_dir.joinpath('protein_data', 'esmfold.jsonlines'))
    sorf_df['ESMFold plddt 90th percentile'] = [np.percentile(esmfold[s.replace('*', '')]['plddt'], 90) if s.replace('*', '') in esmfold else -1 for s in sorf_df['aa'].values]
    sorf_df.to_parquet(cache_dir.joinpath('sorf_df.parq'))
    
    logging.info(f'Start ETL expression data')
    transcripts_to_map = np.concatenate([*sorf_df['transcripts_exact']])
    transcripts_to_map = [str(i).split('.')[0] for i in transcripts_to_map]
    xena, metadata, tissue_pairs = load_xena_transcripts_with_metadata_from_s3(transcripts_to_map)
    groups = create_comparison_groups_xena_tcga_vs_normal(xena, tissue_pairs)
    rows = []
    for cancer, g in groups.items():
        n = xena.loc[g['normal_indices']][xena.columns[6:]].mean(axis=0)
        c = xena.loc[g['cancer_indices']][xena.columns[6:]].mean(axis=0)
        for t in n.index:
            rows.append([t, cancer, g['GTEx Tissue'], g['TCGA Cancer'], n.loc[t], c.loc[t]])
    
    tissue_pairs.to_parquet(cache_dir.joinpath('gtex_tcga_pairs.parq'))
    xena.to_parquet(cache_dir.joinpath('xena.parq'))
    de_genes = read_tcga_de_from_s3('velia-analyses-dev',
                     'VAP_20230329_tcga_differential_expression', output_dir = cache_dir)
    
    xena_expression = xena[xena.columns[6:]]
    xena_metadata = metadata
    
    logging.info(f'Sum expression over each VTX')
    all_transcripts = [i for i in transcripts_to_map if i.startswith('ENST')]
    xena_vtx_sums = xena_expression.T.copy()
    xena_vtx_sums = xena_vtx_sums.loc[xena_vtx_sums.index.intersection(all_transcripts)]
    xena_exact_vtx_sums = {}
    transcripts_in_xena = xena_expression.columns
    for vtx_id, row in tqdm(sorf_df.iterrows()):
        transcripts_parsed = [i.split('.')[0] if i.startswith('ENST') else i for i in row['transcripts_exact']]
        intersection_transcripts = transcripts_in_xena.intersection(transcripts_parsed)
        if len(intersection_transcripts) > 0:
            xena_exact_vtx_sums[vtx_id] = xena_expression[intersection_transcripts].sum(axis=1)
    
    xena_exact_vtx_sums = pd.DataFrame(xena_exact_vtx_sums)
    xena_exact_heatmap_data = process_sums_dataframe_to_heatmap(xena_exact_vtx_sums, xena_metadata)
    xena_metadata.to_parquet(cache_dir.joinpath('xena_metadata.parq'))
    xena_expression.to_parquet(cache_dir.joinpath('xena_app.parq'))
    pickle.dump(xena_exact_heatmap_data, open(cache_dir.joinpath('xena_exact_heatmap.pkl'), 'wb'))
    de_tables_dict, de_metadata = load_de_results(cache_dir, all_transcripts)
    de_tables_dict = {k.split('.')[0]:v for k, v in de_tables_dict.items()}
    tables = []
    for k, v in de_tables_dict.items():
        v['transcript_id'] = k
        tables.append(v)

    pd.concat(tables).rename({'log2FC': 'log2FoldChange', 
                            'Cancer Average': 'Cancer Mean', 
                            'GTEx Average': 'GTEx Mean'}, axis=1).to_sql('transcript_de', sqlite3.connect(cache_dir.joinpath('xena.db')))

    engine = create_engine(f"sqlite:///{cache_dir.joinpath('xena.db')}", future=True)
    with engine.connect() as conn:
        conn.execute(text(f"CREATE INDEX idx_transcript_id_de ON transcript_de (transcript_id);"))
        conn.commit()
    de_metadata.rename({'TCGA Cancer Type': 'TCGA'}, axis=1).to_sql('sample_metadata_de', sqlite3.connect(cache_dir.joinpath('xena.db')), index=False)
    pickle.dump(de_tables_dict, open(cache_dir.joinpath('xena_de_tables_dict.pkl'), 'wb'))
    de_metadata.to_parquet(cache_dir.joinpath('xena_de_metadata.parq'))
    
    

    
