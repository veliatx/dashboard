import boto3
import smart_open
import pandas as pd

def load_xena_transcripts_with_metadata_from_s3(transcripts_to_load = None):
    """
    Wrapper to load data from s3.

    Potentially this will be moved to redshift, but for now this handles loading Xena data.
    
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
    return xena, metadata, tissue_pairs

def load_ccle_from_s3(path_to_ccle):
    pass
    