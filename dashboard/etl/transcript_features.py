import boto3
import smart_open
import pandas as pd
import numpy as np
from tqdm import tqdm
import os
from dashboard.etl import CACHE_DIR, TPM_DESEQ2_FACTOR
from collections import defaultdict

# def load_s3_transcript_DE_stats

def load_CCLE_data_from_s3():
    exp = pd.read_csv('s3://velia-analyses-dev/VAP_20230324_pull_ccle/data/OmicsExpressionProteinCodingGenesTPMLogp1_22Q4.csv')


def load_xena_transcripts_with_metadata_from_s3(transcripts_to_load = None):
    """
    Wrapper to load data from s3.

    Potentially this will be moved to redshift, but for now this handles loading Xena data.

    Loads as TPM + 0.001
    
    """
    xena_transcripts = pd.read_csv('s3://velia-data-dev/VDC_004_annotation/tcga/20230329_tcga_data/xena_transcript_names.csv')
    if transcripts_to_load is None:
        selected_transcripts = [str(t[0]) for t in xena_transcripts.values]
    else:
        overlapping_xena_sorf_transcripts = set(xena_transcripts['Xena Transcripts']).intersection(set(transcripts_to_load))
        selected_transcripts = [str(t) for t in overlapping_xena_sorf_transcripts]

    xena = pd.read_feather('s3://velia-data-dev/VDC_004_annotation/tcga/20230329_tcga_data/xena.feather', columns = ['index'] + selected_transcripts)
    xena.index = xena.pop('index')
    metadata = pd.read_table('s3://velia-data-dev/VDC_004_annotation/tcga/20230329_tcga_data/Xena/TcgaTargetGTEX_phenotype.txt',
                             encoding = "ISO-8859-1", index_col=0)
    metadata = metadata[~metadata['_primary_site'].isna()]
    ordr = xena.index.intersection(metadata.index)
    xena = xena.loc[ordr].copy()
    metadata = metadata.loc[ordr].copy()

    tissue_pairs = pd.read_excel('s3://velia-data-dev/VDC_004_annotation/tcga/20230329_tcga_data/Xena/GTEx_TCGA_comparisons.xlsx')
    tissue_pairs['GTEx Tissue Type'] = tissue_pairs['GTEx Tissue Type'].str.replace('whole ', '')
    tissue_pairs = tissue_pairs[tissue_pairs['GTEx Tissue Type']!='--']

    xena = np.exp2(xena)
    return metadata.merge(xena, left_index=True, right_index=True), metadata, tissue_pairs


def read_tcga_de_from_s3(bucket, object_prefix, output_dir = None,
                         tcga_cancer_codes = None):
    session = boto3.Session()
    s3 = session.resource('s3')
    my_bucket = s3.Bucket(bucket)
    
    cancer_de_results = {}
    for f in my_bucket.objects.filter(Prefix=object_prefix):
        fname = f.key
        if fname.endswith('.csv'):
            cancer_de_results[fname.split('/')[-1].split('_')[0]] = fname
    if tcga_cancer_codes is None:
        tcga_cancer_codes = cancer_de_results.keys()
    tables = {}
    for c in tqdm(tcga_cancer_codes):
        if c in cancer_de_results.keys():
            table = pd.read_csv(f"s3://{bucket}/{cancer_de_results[c]}", index_col=0)
            if not output_dir is None:
                table.to_parquet(f"{output_dir}/{c}_de.parq")
            tables[c] = table    
        else:
            print(f"{c} not found in results")
    return tables


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

def process_sums_dataframe_to_heatmap(xena_vtx_sum_df, xena_metadata_df):
    xena_vtx_sum_df = np.log2(xena_vtx_sum_df + 1)

    xena_vtx_exp_df = xena_metadata_df.merge(xena_vtx_sum_df, left_index=True, right_index=True)

    xena_vtx_sum_df = xena_vtx_sum_df.merge(xena_metadata_df, left_index=True, right_index=True)
    transcript_col_names = [i for i in xena_vtx_sum_df.columns if i not in xena_metadata_df.columns]
    groups = xena_metadata_df.loc[xena_vtx_sum_df.index][['_primary_site', '_study']].apply(lambda x: '-'.join(x), axis=1)
    xena_vtx_sum_df = xena_vtx_sum_df[transcript_col_names].groupby(groups).aggregate(np.mean)

    threshold = .1
    mean_vals = xena_vtx_sum_df.max()
    cols_to_remove = mean_vals[mean_vals < threshold].index
    xena_vtx_sum_df = xena_vtx_sum_df.drop(cols_to_remove, axis=1)
    
    tau_df = xena_vtx_sum_df/xena_vtx_sum_df.max()
    tau = ((1-tau_df).sum())/(tau_df.shape[0]-1)
    tau.name = 'tau'
    xena_tau_df = xena_vtx_sum_df.T.merge(tau, left_index=True, right_index=True)

    return xena_tau_df, xena_vtx_sum_df, xena_vtx_exp_df

def load_de_results(cache_directory, transcripts):
    cache_filez = os.listdir(cache_directory)
    temp_dict = {}
    for f in cache_filez:
        if f.endswith('_de.parq') and not (f=='expression_de.parq'):
            df = pd.read_parquet(os.path.join(cache_directory, f))
            df['transcript'] = df.apply(lambda x: x.name.split('.')[0], axis=1)
            df = df[df['transcript'].isin(transcripts)].copy()
            temp_dict[f.split('_')[0]] = df
            
    de_tables_dict = defaultdict(dict)
    for c, df in tqdm(temp_dict.items()):
        for row in df.itertuples():
            de_tables_dict[row[0]][c] = {'Cancer Average': row._7/TPM_DESEQ2_FACTOR, 'GTEx Average': row._8/TPM_DESEQ2_FACTOR, 
                                         'log2FC': row.log2FoldChange, 'padj': row.padj}
    for t, d in de_tables_dict.items():
        de_tables_dict[t] = pd.DataFrame(d).T
    tcga_gtex_tissue_metadata = pd.read_parquet(os.path.join(cache_directory, 'gtex_tcga_pairs.parq'))
    tcga_gtex_tissue_metadata = tcga_gtex_tissue_metadata.drop_duplicates(['TCGA Cancer Type', 'GTEx Tissue Type']).copy()
    tcga_gtex_tissue_metadata.index = tcga_gtex_tissue_metadata['TCGA Cancer Type']
    return de_tables_dict, tcga_gtex_tissue_metadata