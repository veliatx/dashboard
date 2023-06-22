from sorf_query import *
from transcript_features import *
from amino_acid_features import *
import numpy as np

import jsonlines
import pandas as pd
import sys
NCPU = None
pd.options.display.max_columns = 100
pd.options.display.max_rows = 100
pd.options.display.max_colwidth = 200

if __name__ == '__main__':
    OUTPUT_DIR = '../../cache'
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    input_file = sys.argv[1]
    with open(input_file) as fhandle:
        ids = [int(i.replace('VTX-', '')) for i in fhandle.readlines()]
    # Query DB
    session = base.Session() # connect to db
    orfs = session.query(Orf).filter(Orf.id.in_(ids)).all()
    orfs = [OrfData(i) for i in orfs]
    missing_orfs = set(ids) - set([i.id for i in orfs])
    if len(missing_orfs) > 0:
        print('WARNING: some of the provided IDs were not found in veliadb.', *missing_orfs)

    transcript_matching_results = run_id_mapping_parallel(orfs, NCPU = 1) # implemented for multithreading but had issues with db access in multiprocessing

    # Loop over orfs, and populate sorf_table file with attributes of interest
    transcripts_to_map = []
    with jsonlines.open(os.path.join(OUTPUT_DIR, 'sorf_table.jsonlines'), mode = 'w') as fh:
        for current_orf in tqdm(orfs):
            # nt, aa = extract_nucleotide_sequence_broken_psl_starts(current_orf, reference)
            overlapping_tids = query_overlapping_transcripts(current_orf, session)
            overlapping_tids = [[i.split('.')[0] for i in [t.ensembl_id, t.refseq_id, t.chess_id] if i][0] for t in overlapping_tids]
            # exact_tids = find_seq_substring(nt, transcripts)
            # if len(tids)>0:
            # exact_tids = [i.split('|')[0].split('.')[0] for i in exact_tids]
            exact_tids = transcript_matching_results[current_orf.id][-1]
            nt = transcript_matching_results[current_orf.id][0]
            aa = transcript_matching_results[current_orf.id][1]
            attributes = {
                'chr': current_orf.assembly.ucsc_style_name,
                'vtx': f"VTX-{str(current_orf.id).zfill(7)}",
                'strand': current_orf.strand,
                'start': current_orf.start,
                'end': current_orf.end,
                'nucl': nt,
                'aa': aa,
                'transcripts_exact': exact_tids,
                'transcripts_overlapping': overlapping_tids,
                
            }
            transcripts_to_map+=exact_tids
            transcripts_to_map+=[i.split('.')[0] for i in overlapping_tids if i.startswith('ENST')]
            
            fh.write(attributes)

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
    
    # tissue_pairs.to_parquet(os.path.join(OUTPUT_DIR, 'gtex_tcga_pairs.parq'))
    

    
