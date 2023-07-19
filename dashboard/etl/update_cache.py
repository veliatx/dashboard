from sorf_query import *
from transcript_features import *
from amino_acid_features import *
from dashboard.app import load_mouse_blastp_results
import numpy as np

import jsonlines
import pandas as pd
import sys
import subprocess, shlex
import multiprocessing as mp


NCPU = 16
pd.options.display.max_columns = 100
pd.options.display.max_rows = 100
pd.options.display.max_colwidth = 200

if __name__ == '__main__':
    OUTPUT_DIR = '../../cache'
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
            
    sorf_table = load_jsonlines_table(os.path.join(OUTPUT_DIR, 'sorf_table.jsonlines'), index_col='vtx_id')
    
    transcripts_to_map = np.concatenate((*sorf_table['transcripts_exact'], *sorf_table['transcripts_overlapping']))
    transcripts_to_map = [str(i) for i in transcripts_to_map]
    xena, metadata, tissue_pairs = load_xena_transcripts_with_metadata_from_s3(transcripts_to_map)
    groups = create_comparison_groups_xena_tcga_vs_normal(xena, tissue_pairs)
    rows = []
    for cancer, g in groups.items():
        n = xena.loc[g['normal_indices']][xena.columns[6:]].mean(axis=0)
        c = xena.loc[g['cancer_indices']][xena.columns[6:]].mean(axis=0)
        for t in n.index:
            rows.append([t, cancer, g['GTEx Tissue'], g['TCGA Cancer'], n.loc[t], c.loc[t]])
    
    normal_vs_gtex_expression = pd.DataFrame(rows, columns = ['Transcript', 'TCGA', 'GTEx', 'Description', 'Normal', 'Cancer'])
    normal_vs_gtex_expression.to_parquet(os.path.join(OUTPUT_DIR, 'gtex_tcga_pairs.parq'))
    xena.to_parquet(os.path.join(OUTPUT_DIR, 'xena.parq'))
    de_genes = read_tcga_de_from_s3('velia-analyses-dev',
                     'VAP_20230329_tcga_differential_expression', output_dir = OUTPUT_DIR)
    tissue_pairs.to_parquet(os.path.join(OUTPUT_DIR, 'gtex_tcga_pairs.parq'))
    
    with open(os.path.join(OUTPUT_DIR, 'protein_data', 'protein_tools_input.fasta'), 'w') as fopen:
       for ix, row in sorf_table.iterrows():
           fopen.write(f">{row['vtx_id']}\n{row['aa'].replace('*', '')}\n")
    subprocess.run(shlex.split(f"/home/ec2-user/anaconda/envs/protein_tools/bin/python /home/ec2-user/repos/protein_tools/dashboard_etl.py -i {os.path.abspath(os.path.join(OUTPUT_DIR, 'protein_data', 'protein_tools_input.fasta'))} -o {os.path.abspath(os.path.join(OUTPUT_DIR, 'protein_data'))}"))
    
    

    
