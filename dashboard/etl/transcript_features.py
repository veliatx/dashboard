import boto3
import smart_open
import pandas as pd
import numpy as np

def load_xena_transcripts_with_metadata_from_s3(transcripts_to_load = None):
    """
    Wrapper to load data from s3.

    Potentially this will be moved to redshift, but for now this handles loading Xena data.

    Loads as TPM + 0.001
    
    """
    xena_transcripts = pd.read_csv('s3://velia-data-dev/VDC_004_annotation/20230329_tcga_data/xena_transcript_names.csv')
    if transcripts_to_load is None:
        transcripts_to_map = xena_transcripts['Xena Transcripts']
    else:
        overlapping_xena_sorf_transcripts = set(xena_transcripts['Xena Transcripts']).intersection(set(transcripts_to_load))
    xena = pd.read_feather('s3://velia-data-dev/VDC_004_annotation/20230329_tcga_data/xena.feather', columns = ['index']+list(overlapping_xena_sorf_transcripts))
    xena.index = xena.pop('index')
    metadata = pd.read_table('s3://velia-data-dev/VDC_004_annotation/20230329_tcga_data/Xena/TcgaTargetGTEX_phenotype.txt',
                             encoding = "ISO-8859-1", index_col=0)
    metadata = metadata[~metadata['_primary_site'].isna()]
    ordr = xena.index.intersection(metadata.index)
    xena = xena.loc[ordr].copy()
    metadata = metadata.loc[ordr].copy()

    tissue_pairs = pd.read_excel('s3://velia-data-dev/VDC_004_annotation/20230329_tcga_data/Xena/GTEx_TCGA_comparisons.xlsx')
    tissue_pairs['GTEx Tissue Type'] = tissue_pairs['GTEx Tissue Type'].str.replace('whole ', '')
    tissue_pairs = tissue_pairs[tissue_pairs['GTEx Tissue Type']!='--']

    xena = np.exp2(xena)
    return metadata.merge(xena, left_index=True, right_index=True), metadata, tissue_pairs

def create_comparison_groups_xena_tcga_vs_normal(xena, tissue_pairs):
    groups = {}
    for ix, row in tissue_pairs.iterrows():
        # Normal sample indexes for specific tissues
        
        normal = xena.index[(xena['_primary_site'].str.lower() == row['GTEx Tissue Type'].lower())
                                & (xena['_sample_type'] == 'Normal Tissue')]
        # Cancer sample indexes
        if row['TCGA Cancer Type'] == 'LAML':
            cancer = xena.index[(xena['detailed_category'].str.lower() == row['Description'].lower())
                                    & (xena['_sample_type'] == 'Primary Blood Derived Cancer - Peripheral Blood')]
        else:
            cancer = xena.index[(xena['detailed_category'].str.lower() == row['Description'].lower())
                                    & (xena['_sample_type'] == 'Primary Tumor')]
        groups[row['TCGA Cancer Type']] = {'normal_indices': normal, 'cancer_indices': cancer,
                                       'GTEx Tissue': row['GTEx Tissue Type'], 'TCGA Cancer': row['Description']}
    return groups

def load_ccle_from_s3(path_to_ccle):
    pass
    