from sorf_query import *
from transcript_features import *
from amino_acid_features import *
from dashboard.app import load_mouse_blastp_results
from dashboard.etl.sorf_query import parse_sorf_phase
from dashboard.etl.transcript_features import load_xena_transcripts_with_metadata_from_s3, process_sums_dataframe_to_heatmap, create_comparison_groups_xena_tcga_vs_normal, read_tcga_de_from_s3, load_de_results
from dashboard.data_load import load_esmfold
import numpy as np

import pandas as pd
import sys
import multiprocessing as mp
import subprocess
import shlex
import json
import pickle
import os

from dashboard.etl import CACHE_DIR

NCPU = 16
pd.options.display.max_columns = 100
pd.options.display.max_rows = 100
pd.options.display.max_colwidth = 200

if __name__ == '__main__':
    # OUTPUT_DIR = '../../cache_update'
    # CACHE_DIR = OUTPUT_DIR
    OUTPUT_DIR = CACHE_DIR
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    os.makedirs(os.path.join(OUTPUT_DIR, 'protein_data'), exist_ok=True)

    input_file = sys.argv[1]
    with open(input_file) as fhandle:
        ids = [int(i.replace('VTX-', '')) for i in fhandle.readlines()]
    # Query DB
    session = base.Session() # connect to db
    orfs = session.query(Orf).filter(Orf.id.in_(ids)).all()
    missing_orfs = set(ids) - set([i.id for i in orfs])
    if len(missing_orfs) > 0:
        print('WARNING: some of the provided IDs were not found in veliadb.', *missing_orfs)
    with open(os.path.join(OUTPUT_DIR, 'sorf_table.jsonlines'), 'w') as fopen:
        with mp.Pool(NCPU) as ppool:
            for r in tqdm(ppool.imap(parallel_sorf_query, ids), total=len(ids)):
                fopen.write(json.dumps(r))
                fopen.write('\n')
            
    sorf_df = load_jsonlines_table(os.path.join(OUTPUT_DIR, 'sorf_table.jsonlines'), index_col='vtx_id')
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


    phase_ids, phase_entries = parse_sorf_phase(sorf_df, session)
    sorf_df['screening_phase_id'] = phase_ids
    sorf_df['screening_phase'] = phase_entries
    sorf_df['protein_xrefs'] = protein_xrefs
    sorf_df['aa_length'] = sorf_df.apply(lambda x: len(x.aa), axis=1)

    sorf_df = fix_missing_phase_ids(sorf_df)
    # sorf_df = sorf_df[sorf_df['screening_phase'] != '-1']

    cols = ['show_details', 'vtx_id', 'aa_length', 'screening_phase_id', 'screening_phase', 'ucsc_track', 
            'source', 'orf_xrefs', 'protein_xrefs', 'gene_xrefs', 'transcript_xrefs',  
            'transcripts_exact', 'aa', 'nucl', 
            'index_copy', 'genscript_id', 'chr', 'strand', 'start', 'end', 
            'chrom_starts', 'block_sizes', 'phases',]

    sorf_df = sorf_df[cols]
    session.close()
    
    # Munge sorf_df into format consistent with historical app.
    with open(os.path.join(OUTPUT_DIR, 'protein_data', 'protein_tools_input.fasta'), 'w') as fopen:
       for ix, row in sorf_df.iterrows():
           fopen.write(f">{row['vtx_id']}\n{row['aa'].replace('*', '')}\n")
    python_executable = '/home/ec2-user/anaconda/envs/protein_tools/bin/python'
    # subprocess.run(shlex.split(f"{python_executable} /home/ec2-user/repos/protein_tools/dashboard_etl.py -i {os.path.abspath(os.path.join(OUTPUT_DIR, 'protein_data', 'protein_tools_input.fasta'))} -o {os.path.abspath(os.path.join(OUTPUT_DIR, 'protein_data'))}"))
    # Massage table to standard format
    _, blastp_table = load_mouse_blastp_results(CACHE_DIR = OUTPUT_DIR)
    sorf_df = sorf_df.merge(pd.DataFrame(blastp_table).T, left_index=True, right_index=True, how='left')
    protein_scores = pd.read_csv(os.path.join(OUTPUT_DIR, 'protein_data', 'sequence_features_scores.csv'), index_col=0)

    sorf_df = sorf_df.merge(protein_scores[['Deepsig', 'SignalP 6slow', 'SignalP 5b', 'SignalP 4.1']],
                    left_index=True, right_index=True, how='left')
    protein_strings = pd.read_csv(os.path.join(OUTPUT_DIR, 'protein_data', 'sequence_features_strings.csv'), index_col=0)
    protein_cutsite = protein_strings.apply(lambda x: x.str.find('SO')+1).replace(0, -1).drop('Sequence', axis=1)
    sorf_df = sorf_df.merge(protein_cutsite,
                    left_index=True, right_index=True, how='left', suffixes=('_score', '_cut'))

    id_data = pd.read_csv('s3://velia-data-dev/VDC_001_screening_collections/all_phases/interim_phase1to7_non-sigp_20230723.csv')

    id_data.index = id_data['vtx_id']
    sorf_df = sorf_df.merge(id_data[['trans1',
            'trans2', 'trans3', 'sec1', 'sec2', 'sec3', 'translated_mean',
            'secreted_mean', 'translated', 'swissprot_isoform', 'ensembl_isoform',
            'refseq_isoform', 'phylocsf_58m_avg', 'phylocsf_58m_max', 'phylocsf_58m_min',
            'phylocsf_vals']], left_index=True, right_index=True, how='left')

    with open(DATA_DIR / 'all_secreted_phase1to7.txt', 'r') as f:
        secreted_ids = [i.strip() for i in f.readlines()]

    sorf_df.insert(int(np.where(sorf_df.columns=='translated')[0][0]), 'secreted', [True if i in secreted_ids else False for i in sorf_df.index])
    sorf_df['transcripts_exact'] = [tuple(i) for i in sorf_df['transcripts_exact']]
    # add esmfold data
    esmfold = load_esmfold()
    sorf_df['ESMFold plddt 90th percentile'] = [np.percentile(esmfold[s.replace('*', '')]['plddt'], 90) if s.replace('*', '') in esmfold else -1 for s in sorf_df['aa'].values]
    sorf_df.to_parquet(os.path.join(OUTPUT_DIR, 'sorf_df.parq'))
    # Done formatting
    
    # Start ETL expression data
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
    
    tissue_pairs.to_parquet(os.path.join(OUTPUT_DIR, 'gtex_tcga_pairs.parq'))
    xena.to_parquet(os.path.join(OUTPUT_DIR, 'xena.parq'))
    de_genes = read_tcga_de_from_s3('velia-analyses-dev',
                     'VAP_20230329_tcga_differential_expression', output_dir = OUTPUT_DIR)
    
    xena_expression = xena[xena.columns[6:]]
    xena_metadata = metadata
    
    all_transcripts = [i for i in transcripts_to_map if i.startswith('ENST')]
    # Sum expression over each VTX            
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
    xena_metadata.to_parquet(os.path.join(OUTPUT_DIR, 'xena_metadata.parq'))
    xena_expression.to_parquet(os.path.join(OUTPUT_DIR, 'xena_app.parq'))
    pickle.dump(xena_exact_heatmap_data, open(os.path.join(OUTPUT_DIR, 'xena_exact_heatmap.pkl'), 'wb'))
    de_tables_dict, de_metadata = load_de_results(all_transcripts)
    de_tables_dict = {k.split('.')[0]:v for k, v in de_tables_dict.items()}
    pickle.dump(de_tables_dict, open(os.path.join(OUTPUT_DIR, 'xena_de_tables_dict.pkl'), 'wb'))
    de_metadata.to_parquet(os.path.join(OUTPUT_DIR, 'xena_de_metadata.parq'))
    
    

    
