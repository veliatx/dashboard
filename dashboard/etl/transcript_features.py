"""Module for loading and processing transcript expression data from various sources.

This module provides functions to:
- Load transcript expression data from TCGA/GTEx/CCLE
- Process differential expression results
- Calculate tissue specificity scores
- Handle data transformations and metadata merging

The module interfaces with S3 to load data and processes it into formats suitable
for downstream analysis.
"""

import boto3
import smart_open
import pandas as pd
import numpy as np
from tqdm import tqdm
import os
from dashboard.etl import CACHE_DIR, TPM_DESEQ2_FACTOR
from collections import defaultdict
from typing import Dict, List, Optional, Tuple, Any


def load_CCLE_data_from_s3() -> pd.DataFrame:
    """Load CCLE expression data from S3.
    
    Returns:
        DataFrame containing CCLE expression data
    """
    exp = pd.read_csv('s3://velia-analyses-dev/VAP_20230324_pull_ccle/data/OmicsExpressionProteinCodingGenesTPMLogp1_22Q4.csv')
    return exp


def load_xena_transcripts_with_metadata_from_s3(
    transcripts_to_load: Optional[List[str]] = None
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Load TCGA/GTEx transcript expression data and metadata from S3.

    Args:
        transcripts_to_load: Optional list of transcript IDs to load. If None, loads all transcripts.

    Returns:
        Tuple containing:
        - DataFrame with merged expression and metadata
        - DataFrame with just metadata
        - DataFrame with tissue type mapping between GTEx and TCGA
    """
    xena_transcripts = pd.read_csv('s3://velia-data-dev/VDC_004_annotation/tcga/20230329_tcga_data/xena_transcript_names.csv')
    if transcripts_to_load is None:
        selected_transcripts = [str(t[0]) for t in xena_transcripts.values]
    else:
        overlapping_xena_sorf_transcripts = set(xena_transcripts['Xena Transcripts']).intersection(set(transcripts_to_load))
        selected_transcripts = [str(t) for t in overlapping_xena_sorf_transcripts]

    xena = pd.read_feather('s3://velia-data-dev/VDC_004_annotation/tcga/20230329_tcga_data/xena.feather', columns=['index'] + selected_transcripts)
    xena.index = xena.pop('index')
    metadata = pd.read_table('s3://velia-data-dev/VDC_004_annotation/tcga/20230329_tcga_data/Xena/TcgaTargetGTEX_phenotype.txt',
                           encoding="ISO-8859-1", index_col=0)
    metadata = metadata[~metadata['_primary_site'].isna()]
    ordr = xena.index.intersection(metadata.index)
    xena = xena.loc[ordr].copy()
    metadata = metadata.loc[ordr].copy()

    tissue_pairs = pd.read_excel('s3://velia-data-dev/VDC_004_annotation/tcga/20230329_tcga_data/Xena/GTEx_TCGA_comparisons.xlsx')
    tissue_pairs['GTEx Tissue Type'] = tissue_pairs['GTEx Tissue Type'].str.replace('whole ', '')
    tissue_pairs = tissue_pairs[tissue_pairs['GTEx Tissue Type'] != '--']

    xena = np.exp2(xena)
    return metadata.merge(xena, left_index=True, right_index=True), metadata, tissue_pairs


def read_tcga_de_from_s3(
    bucket: str,
    object_prefix: str,
    output_dir: Optional[str] = None,
    tcga_cancer_codes: Optional[List[str]] = None
) -> Dict[str, pd.DataFrame]:
    """Read differential expression results from S3 for TCGA cancer types.

    Args:
        bucket: S3 bucket name
        object_prefix: Prefix path in bucket to search
        output_dir: Optional directory to save parquet files
        tcga_cancer_codes: Optional list of TCGA cancer type codes to load

    Returns:
        Dictionary mapping cancer types to their differential expression results
    """
    session = boto3.Session()
    s3 = session.resource('s3')
    my_bucket = s3.Bucket(bucket)
    
    cancer_de_results = {}
    for f in my_bucket.objects.filter(Prefix=object_prefix):
        fname = f.key
        if fname.endswith('.csv'):
            cancer_de_results[fname.split('/')[-1].split('_')[0]] = fname

    if tcga_cancer_codes is None:
        tcga_cancer_codes = list(cancer_de_results.keys())

    tables = {}
    for c in tqdm(tcga_cancer_codes):
        if c in cancer_de_results:
            table = pd.read_csv(f"s3://{bucket}/{cancer_de_results[c]}", index_col=0)
            if output_dir is not None:
                table.to_parquet(f"{output_dir}/{c}_de.parq")
            tables[c] = table    
        else:
            print(f"{c} not found in results")
    return tables


def create_comparison_groups_xena_tcga_vs_normal(
    xena: pd.DataFrame,
    tissue_pairs: pd.DataFrame
) -> Dict[str, Dict[str, Any]]:
    """Create sample groupings for TCGA tumor vs GTEx normal comparisons.

    Args:
        xena: DataFrame containing sample metadata
        tissue_pairs: DataFrame mapping GTEx tissues to TCGA cancer types

    Returns:
        Dictionary mapping cancer types to their sample groupings
    """
    groups = {}
    for ix, row in tissue_pairs.iterrows():
        normal = xena.index[(xena['_primary_site'].str.lower() == row['GTEx Tissue Type'].lower())
                           & (xena['_sample_type'] == 'Normal Tissue')]
        
        if row['TCGA Cancer Type'] == 'LAML':
            cancer = xena.index[(xena['detailed_category'].str.lower() == row['Description'].lower())
                               & (xena['_sample_type'] == 'Primary Blood Derived Cancer - Peripheral Blood')]
        else:
            cancer = xena.index[(xena['detailed_category'].str.lower() == row['Description'].lower())
                               & (xena['_sample_type'] == 'Primary Tumor')]
            
        groups[row['TCGA Cancer Type']] = {
            'normal_indices': normal,
            'cancer_indices': cancer,
            'GTEx Tissue': row['GTEx Tissue Type'],
            'TCGA Cancer': row['Description']
        }
    return groups


def process_sums_dataframe_to_heatmap(
    xena_vtx_sum_df: pd.DataFrame,
    xena_metadata_df: pd.DataFrame
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Process expression data into heatmap format and calculate tissue specificity.

    Args:
        xena_vtx_sum_df: DataFrame containing expression values
        xena_metadata_df: DataFrame containing sample metadata

    Returns:
        Tuple containing:
        - DataFrame with tau scores
        - DataFrame with processed expression values
        - DataFrame with merged expression and metadata
    """
    xena_vtx_sum_df = np.log2(xena_vtx_sum_df + 1)

    xena_vtx_exp_df = xena_metadata_df.merge(xena_vtx_sum_df, left_index=True, right_index=True)

    xena_vtx_sum_df = xena_vtx_sum_df.merge(xena_metadata_df, left_index=True, right_index=True)
    transcript_col_names = [i for i in xena_vtx_sum_df.columns if i not in xena_metadata_df.columns]
    groups = xena_metadata_df.loc[xena_vtx_sum_df.index][['_primary_site', '_study']].apply(lambda x: '-'.join(x), axis=1)
    xena_vtx_sum_df = xena_vtx_sum_df[transcript_col_names].groupby(groups).aggregate(np.mean)

    threshold = 0.1
    mean_vals = xena_vtx_sum_df.max()
    cols_to_remove = mean_vals[mean_vals < threshold].index
    xena_vtx_sum_df = xena_vtx_sum_df.drop(cols_to_remove, axis=1)
    
    tau_df = xena_vtx_sum_df/xena_vtx_sum_df.max()
    tau = ((1-tau_df).sum())/(tau_df.shape[0]-1)
    tau.name = 'tau'
    xena_tau_df = xena_vtx_sum_df.T.merge(tau, left_index=True, right_index=True)

    return xena_tau_df, xena_vtx_sum_df, xena_vtx_exp_df


def load_de_results(
    cache_directory: str,
    transcripts: List[str]
) -> Tuple[Dict[str, pd.DataFrame], pd.DataFrame]:
    """Load differential expression results from cache.

    Args:
        cache_directory: Path to cache directory
        transcripts: List of transcript IDs to load

    Returns:
        Tuple containing:
        - Dictionary mapping transcripts to their DE results
        - DataFrame with tissue type mappings
    """
    cache_filez = os.listdir(cache_directory)
    temp_dict = {}
    for f in cache_filez:
        if f.endswith('_de.parq') and not (f == 'expression_de.parq'):
            df = pd.read_parquet(os.path.join(cache_directory, f))
            df['transcript'] = df.apply(lambda x: x.name.split('.')[0], axis=1)
            df = df[df['transcript'].isin(transcripts)].copy()
            temp_dict[f.split('_')[0]] = df
            
    de_tables_dict = defaultdict(dict)
    for c, df in tqdm(temp_dict.items()):
        for row in df.itertuples():
            de_tables_dict[row[0]][c] = {
                'Cancer Average': row._7/TPM_DESEQ2_FACTOR,
                'GTEx Average': row._8/TPM_DESEQ2_FACTOR,
                'log2FC': row.log2FoldChange,
                'padj': row.padj
            }
    for t, d in de_tables_dict.items():
        de_tables_dict[t] = pd.DataFrame(d).T

    tcga_gtex_tissue_metadata = pd.read_parquet(os.path.join(cache_directory, 'gtex_tcga_pairs.parq'))
    tcga_gtex_tissue_metadata = tcga_gtex_tissue_metadata.drop_duplicates(['TCGA Cancer Type', 'GTEx Tissue Type']).copy()
    tcga_gtex_tissue_metadata.index = tcga_gtex_tissue_metadata['TCGA Cancer Type']
    
    return de_tables_dict, tcga_gtex_tissue_metadata