import json
import jsonlines
import pickle
import numpy as np

import pandas as pd
import pyarrow.parquet as pq
import streamlit as st
import sqlite3

from dashboard.util import convert_list_string
from dashboard.etl import CACHE_DIR, DATA_DIR, NOTEBOOK_DATA_DIR

from ast import literal_eval
from collections import defaultdict
from pathlib import Path

@st.cache_data()
def load_autoimmune_atlas():
    """
    """
    meta = pd.read_parquet("s3://velia-athena-dev/expression_atlas_v1_metadata.parq")
    meta['dashboard_group'] = meta['atlas_group']
    return meta

@st.cache_data()
def load_autoimmune_metadata():

    """
    Temporary table for joining contrasts from the autoimmune studies to specific samples.
    To be replaced by table in autoimmune sqlitedb.
    """
    with sqlite3.connect(DATA_DIR.joinpath('autoimmune_expression_atlas_v1.db')) as sqliteConnection:
        contrast_samples = pd.read_sql("SELECT * FROM transcript_contrast", sqliteConnection)#pd.read_csv(CACHE_DIR.joinpath('de_contrast_table.csv'))
        sample_meta_df = pd.read_sql(
                                f"""SELECT
                                        sample_id, 
                                        sample_condition_1,
                                        sample_condition_2,
                                        sample_condition_3,
                                        sample_type_1,
                                        sample_type_2
                                    FROM sample_metadata""",
                                sqliteConnection,
                            )
            
    contrast_samples.index = contrast_samples.apply(lambda x: f"{x['velia_study']} -- {x['contrast']}", axis=1)
    contrast_samples.index.name = 'study -- contrast'
    contrast_samples['display'] = contrast_samples.apply(
                                    lambda x: f"{x.velia_study} -- " +
                                        f"{x.contrast.upper().split('VS')[0 if x.contrast_side == 'left' else 1].strip('_')}", 
                                    axis=1,
                                    )
    return contrast_samples, sample_meta_df

@st.cache_data(ttl='24h')
def pull_s3_normed_counts(
                    study:str,
                    contrast:str,
                    ):
    """
    """
    normed_counts_df = pd.read_csv(
        f's3://velia-piperuns-dev/expression_atlas/v1/{study}/de_results/{study}_transcript_{contrast}.csv',
        index_col=0,
        )
    normed_counts_df = normed_counts_df.loc[:,
                            (normed_counts_df.columns.str.startswith('condition') & 
                            (~normed_counts_df.columns.str.contains('meannormed'))) |
                            normed_counts_df.columns.isin(('padj','log2FoldChange'))
                            ]
    normed_counts_df.rename(lambda x: x.split('_')[-1], inplace=True, axis=1)
    normed_counts_df.loc['velia_study',:] = study
    normed_counts_df.loc['contrast',:] = contrast
    normed_counts_df['padj'].fillna(1.0, inplace=True)
                        
    return normed_counts_df

@st.cache_data()
def load_sorf_df_conformed():
    """
    """
    df = pd.read_parquet(CACHE_DIR.joinpath('sorf_df.parq'))
    df.drop('nonsignal_seqs', axis=1, inplace=True)
    #df = df[df['aa_length'] <= 150].copy()

    df = add_temp_ms_ribo_info(df) # Fine it's not cache it's data until ribo-seq info is in veliadb or something
    
    df = add_temp_isoform_info(df) # Converted to use ETL format

    df = add_temp_tblastn_info(df) # Uses ETL Version
    
    df = filter_riboseq(df) # Fine?? doesn't rely on any data/cache info, but could perform this in cache update too

    df = add_temp_gwas(df)

    # df = add_temp_feature_info(df) # Columns added directly to the sequence_features_strings.csv file
    feature_df = pd.read_csv(CACHE_DIR.joinpath('protein_data', 'sequence_features_strings.csv'), index_col=0)
    feature_cols = ['nonsignal_seqs', 'DeepTMHMM_prediction', 'DeepTMHMM_length']
    df = df.merge(feature_df[feature_cols], left_index=True, right_index=True, how='left')
    df.index.name = 'vtx_id'
    isoform_cols = ['swissprot_isoform', 'ensembl_isoform', 'refseq_isoform']
    df[isoform_cols] = df[isoform_cols].apply(lambda x: [str(i) if isinstance(i, np.ndarray) else 'None' for i in x], axis=0)
    df['vtx_id'] = df.index
    
    df = add_temp_nonsig_cons_info(df) # generates core output in run_protein_search_tools.py

    df = reorder_table_cols(df)
    df.rename({'secreted': 'secreted_hibit',
               'translated': 'translated_hibit'}, axis=1, inplace=True)
    
    df[['start', 'end']] = df[['start', 'end']].astype(int)
     
    return df


def reorder_table_cols(df):
    """
    """
    view_cols = [
        'show_details', 'vtx_id', 'aa_length', 'ucsc_track', 'source', 
        'screening_phase_id', 'screening_phase', 'genscript_id',
        'orf_xrefs', 'protein_xrefs', 'gene_xrefs', 'transcript_xrefs', 
        'transcripts_exact', 'aa', 'nucl', 'chr', 'strand', 'start', 'end',
        'chrom_starts', 'block_sizes', 'phases', 
        'blastp_subject', 'blastp_hit_description',
        'blastp_align_length', 'blastp_align_identity', 
        'tblastn_hit_id', 'tblastn_description',
        'tblastn_align_length', 'tblastn_align_identity', 
        'Deepsig_score', 'SignalP 6slow_score', 'SignalP 5b_score', 'SignalP 4.1_score', 
        'Deepsig_cut', 'SignalP 6slow_cut', 'SignalP 5b_cut', 'SignalP 4.1_cut', 
        'Phobius', 'DeepTMHMM', 
        'translated_mean', 'secreted_mean', 'secreted', 'translated', 
        'phylocsf_58m_avg', 'phylocsf_58m_max', 'phylocsf_58m_min', 'ESMFold plddt 90th percentile',
        'MS_evidence', 'swissprot_isoform', 'ensembl_isoform', 'refseq_isoform', 
        'Ribo-Seq RPKM Support', 'Ribo-Seq sORF',
        'nonsignal_seqs', 'DeepTMHMM_prediction', 'DeepTMHMM_length',
        'nonsig_blastp_align_identity', 'nonsig_tblastn_align_identity',
        'SNPS', 'MAPPED_TRAIT', 'P-VALUE']
    
    return df[view_cols]


def add_temp_gwas(df):
    """
    """
    sorf_gwas_df = pd.read_csv(DATA_DIR.joinpath('embl_nhgri_gwas_v1.0.2.dashboard_sorfs.tsv'), sep='\t')
    agg_df = sorf_gwas_df[['vtx_id', 'SNPS', 'MAPPED_TRAIT', 'P-VALUE']].groupby('vtx_id').aggregate(list)
    df = df.merge(agg_df, left_index=True, right_index=True, how='left')
    
    return df

def add_temp_nonsig_cons_info(df):
    ""
    ""
    tdf = pd.read_parquet(CACHE_DIR.joinpath('protein_data', 'nonsignal_seq_blast_tblastn.parq'))

    df = df.merge(tdf, left_index=True, right_index=True, how='left')

    return df


def add_temp_feature_info(df):
    """
    """
    feature_df = pd.read_csv(CACHE_DIR.joinpath('protein_data', 'sequence_features_strings.csv'), index_col=0)

    signal_cols = ['SignalP 6slow', 'SignalP 4.1', 'SignalP 5b', 'Deepsig']
    nonsignal_seqs = []

    for i, row in feature_df.iterrows():
        nonsignal_aa = ''
        for col in signal_cols:
            if row[col].startswith('S'):
                signal_len = row[col].count('S')
                nonsignal_aa = row['Sequence'][signal_len:-1]
                break
        nonsignal_seqs.append(nonsignal_aa)

    feature_df['nonsignal_seqs'] = nonsignal_seqs

    feature_df['DeepTMHMM_prediction'] = feature_df.apply(lambda x: 'M' in x['DeepTMHMM'] or 'B' in x['DeepTMHMM'], axis=1)
    feature_df['DeepTMHMM_length'] = feature_df.apply(lambda x: x['DeepTMHMM'].count('M'), axis=1)

    feature_cols = ['nonsignal_seqs', 'DeepTMHMM_prediction', 'DeepTMHMM_length']
    df = df.merge(feature_df[feature_cols], left_index=True, right_index=True, how='left')
    
    df.index.name = 'vtx_id'
    df['vtx_id'] = df.index

    return df


def add_temp_tblastn_info(df):
    """
    """
    tblastn_df = pd.read_parquet(CACHE_DIR.joinpath('protein_data', 'mouse_tblastn.parq'))
    df = df.merge(tblastn_df[['tblastn_hit_id', 'tblastn_description',
                              'tblastn_align_length', 'tblastn_align_identity', 'tblastn_evalue']], how='left', left_index=True, right_index=True)               
    df.drop('phylocsf_vals', axis=1, inplace=True)
    
    return df


def add_temp_isoform_info(df):
    """
    """
    df.drop(columns=['swissprot_isoform', 
                        'ensembl_isoform', 
                        'refseq_isoform'], inplace=True)

    isoforms = pd.read_parquet(CACHE_DIR.joinpath('protein_data', 'isoforms_search.parq'))

    isoform_cols = ['swissprot_isoform', 'ensembl_isoform', 'refseq_isoform'] 
    df = df.merge(isoforms, left_index=True, right_index=True, how='left')
    # df[isoform_cols] = df[isoform_cols].apply(lambda x: [literal_eval(y) for y in x])

    return df


def add_temp_ms_ribo_info(df):
    """
    """
    ribo_df = pd.read_excel(DATA_DIR.joinpath('Secreted_mP_Riboseq_SAF.xlsx'))
    ribo_vtx = set(ribo_df[ribo_df['manual_check'] == 1]['vtx_id'])
    ccle_df = pd.read_excel(DATA_DIR.joinpath('SummaryIdentification_CCLE_strongerConfidence.xlsx'), index_col=0)
    gtex_df = pd.read_excel(DATA_DIR.joinpath('SummaryIdentification_GTEX_strongerConfidence.xlsx'), index_col=0)

    ccle_vtx = set(ccle_df['vtx_id'])
    gtex_vtx = set(gtex_df['vtx_id'])

    ms_vtx = gtex_vtx.union(ccle_vtx)
    support_vtx = ms_vtx.union(ribo_vtx)

    df['manual_riboseq'] = df.apply(lambda x: True if x.vtx_id in ribo_vtx else False, axis=1)
    df['MS_evidence'] = df.apply(lambda x: True if x.vtx_id in ms_vtx else False, axis=1)
    df['MS or Riboseq'] = df.apply(lambda x: True if x.vtx_id in support_vtx else False, axis=1)

    return df


def filter_riboseq(df):
    """
    Temporary function to enforce ribo-seq filtering
    for primary entries in collection
    """
    df['Ribo-Seq sORF'] = (
        (df['source'].apply(lambda x: 'gencode_riboseq' in x)) | \
        (df['source'].apply(lambda x: 'velia_phase1_Bona fide' in x)) | \
        (df['source'].apply(lambda x: 'velia_phase1_Chen' in x)) | \
        (df['source'].apply(lambda x: 'velia_phase1_Prensner' in x)) | \
        (df['source'].apply(lambda x: 'velia_phase2_Chang_Saghatelian' in x)) | \
        (df['source'].apply(lambda x: 'velia_phase2_Chothani2022_SignalP' in x)) | \
        (df['source'].apply(lambda x: 'velia_phase2_Bianca_Chen' in x)) | \
        (df['source'].apply(lambda x: 'velia_phase2_Bonafide_Bianca' in x)) | \
        (df['source'].apply(lambda x: 'velia_phase2_Cao_Slavoff_MINAS60' in x)) | \
        (df['source'].apply(lambda x: 'velia_phase2_Rat_Cardiac_Huang' in x)) | \
        (df['source'].apply(lambda x: 'velia_phase2_Mudge2022_SignalP' in x)) | \
        (df['source'].apply(lambda x: 'velia_phase5_Blume_Mudge' in x)) | \
        (df['source'].apply(lambda x: 'velia_phase5_bona fide' in x)) | \
        (df['source'].apply(lambda x: 'velia_phase6_plasma_mass_spec' in x)) | \
        (df['source'].apply(lambda x: 'velia_phase6_public_mass_spec' in x)) | \
        (df['source'].apply(lambda x: 'velia_phase7_Ribo-seq_PBMC_LPS_R848' in x)) | \
        (df['source'].apply(lambda x: 'velia_phase9_orfrater' in x)) | \
        (df['source'].apply(lambda x: 'velia_phase9_Olsen' in x)) | \
        (df['source'].apply(lambda x: 'velia_phase7_Ribo-seq_PBMC_LPS_R848' in x)) | \
        (df['source'].apply(lambda x: 'velia_phase10_riboseq_230114' in x)) | \
        (df['screening_phase'] == 'Not Screened') |
        (df['orf_xrefs'].astype(str).str.contains('RibORF')))
    
    ribo_df = df[df['Ribo-Seq sORF']].copy()
    x = ribo_df.groupby('aa').aggregate(list)

    vtx_to_keep = []

    for i, row in x.iterrows():
        vtx_id = ''
        
        if len(row['vtx_id']) > 1:
                        
            for j, phase in enumerate(row['screening_phase']):
                if 'phase' in phase.lower():
                    vtx_id = row['vtx_id'][j]
                    
            if vtx_id == '':
                vtx_id = row['vtx_id'][0]
        else:
            vtx_id = row['vtx_id'][0]
            
        vtx_to_keep.append(vtx_id)
        
    ribo_df = ribo_df[ribo_df['vtx_id'].isin(vtx_to_keep)].copy()

    ribo_aa = set(ribo_df['aa'])

    non_ribo_df = df[~df['Ribo-Seq sORF']].copy()
    non_ribo_df = non_ribo_df[~non_ribo_df['aa'].isin(ribo_aa)]

    df = pd.concat([ribo_df, non_ribo_df])

    return df


@st.cache_data()
def load_protein_feature_string_representations():
    """
    """
    df = pd.read_csv(CACHE_DIR.joinpath('protein_data', 'sequence_features_strings.csv'), index_col=0).T
    return df


@st.cache_data()
def load_xena_transcripts(transcripts):
    """
    """
    xena_expression = pd.read_parquet(CACHE_DIR.joinpath('xena_app.parq'), columns = transcripts)
    return xena_expression


@st.cache_data()
def load_xena_metadata():
    """
    """
    xena_metadata = pd.read_parquet(CACHE_DIR.joinpath('xena_metadata.parq'))
    xena_metadata['dashboard_group'] = list(map(lambda x: '-'.join(map(str, x)), xena_metadata[['_primary_site', '_study']].values))
    parq_metadata = pq.read_metadata(CACHE_DIR.joinpath('xena_app.parq'))
    xena_transcripts = set(parq_metadata.schema.names)
    return xena_metadata, xena_transcripts


@st.cache_data()
def load_xena_heatmap_data():
    with open(CACHE_DIR.joinpath('xena_exact_heatmap.pkl'), 'rb') as infile:
        xena_tau_df, xena_vtx_sum_df, xena_vtx_exp_df = pickle.load(infile)

    non_vtx_cols = [x for x in list(xena_vtx_exp_df.columns) if x[0:3] != 'VTX']
    vtx_cols = [x for x in list(xena_vtx_exp_df.columns) if x[0:3] == 'VTX']

    index_cols = xena_vtx_exp_df[non_vtx_cols].T.drop_duplicates().T
    data_cols = xena_vtx_exp_df[vtx_cols]

    xena_vtx_exp_df = index_cols.merge(data_cols, left_index=True, right_index=True)

    return xena_tau_df, xena_vtx_sum_df, xena_vtx_exp_df


@st.cache_data()
def load_esmfold(file_path):
    """
    """
    esmfold = {}
    with open(file_path) as fopen:
        j_reader = jsonlines.Reader(fopen)
        for l in j_reader:
            esmfold[l['sequence']] = l
    return esmfold


@st.cache_data()
def load_mouse_blastp_results(CACHE_DIR=CACHE_DIR):
    """
    """
    with open(CACHE_DIR.joinpath('protein_data', 'blastp.results.pkl'), 'rb') as fopen:
        hits_per_query, sorf_table_data = pickle.load(fopen)
    return hits_per_query, sorf_table_data

@st.cache_data()
def load_phylocsf_data():
    """
    """
    pcsf = pd.read_csv(DATA_DIR.joinpath(f"interim_phase1to7_all_phylocsf-vals_20230628.csv"), index_col=0)
    pcsf['phylocsf_vals'] = pcsf['phylocsf_vals'].apply(convert_list_string)
    pcsf = pcsf[['phylocsf_58m_avg', 'phylocsf_58m_max',
           'phylocsf_58m_min', 'phylocsf_58m_std', 'phylocsf_vals']]
    return pcsf
