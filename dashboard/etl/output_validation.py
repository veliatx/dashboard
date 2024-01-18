import pandas as pd
import numpy as np

def check_esmfold_completeness(esmfold, sorf_df):
    missing_vtx = []
    for ix, row in sorf_df.iterrows():
        if row['aa'] not in esmfold:
            missing_vtx.append(ix)
    if len(missing_vtx) == 0:
        return True, None
    else:
        return False, missing_vtx
    
def check_blast_csv(blast_csv, sorf_df):
    sorf_table_order = list(sorf_df['vtx_id'])
    indexes_of_blast_results = []
    for ix, row in blast_csv.iterrows():
        indexes_of_blast_results.append(sorf_table_order.index(ix))
    fraction_covered = max(indexes_of_blast_results)/len(sorf_table_order)
    largest_gap = max(np.diff(sorted(indexes_of_blast_results)))
    return fraction_covered, largest_gap

def check_sorf_df_duplicates(sorf_df):
    return [(i, c) for i, c in sorf_df.index.value_counts().items() if c>1]
