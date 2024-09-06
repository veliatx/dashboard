import json
import jsonlines
import pickle
import numpy as np

import pandas as pd
import pyarrow.parquet as pq
import streamlit as st
import sqlite3

from dashboard.util import convert_list_string
from dashboard.etl import CACHE_DIR, DATA_DIR, HMMER_S3_LOC, NOTEBOOK_DATA_DIR

from ast import literal_eval
from collections import defaultdict
from pathlib import Path
from sqlalchemy import and_
from veliadb import base
from veliadb.base import Protein, Orf, OrfXref, Dataset


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

@st.cache_data()
def load_sorf_df_conformed():
    """
    """
    df = pd.read_parquet(CACHE_DIR.joinpath('sorf_df.parq'))
    df.drop('nonsignal_seqs', axis=1, inplace=True)
    
    df = add_temp_ms_ribo_info(df) # Fine it's not cache it's data until ribo-seq info is in veliadb
    
    df = add_temp_isoform_info(df) # Converted to use ETL format
    
    df = add_temp_tblastn_info(df) # Uses ETL Version

    df = add_temp_source(df)

    df = add_temp_hibit(df)

    # TODO: These are getting filtered because Ribo-Seq sORF is not flagged for these transcripts. 
    keep_vtx_no_translated = \
    """VTX-0850289
    VTX-0087278
    VTX-0774612
    VTX-0850643
    VTX-0851353
    VTX-0060798
    VTX-0850465
    VTX-0069174
    VTX-0699909
    VTX-0852737
    VTX-0851971
    VTX-0851127
    VTX-0850174
    VTX-0850841
    VTX-0852738
    VTX-0851309
    VTX-0826742
    VTX-0852041
    VTX-0015094
    VTX-0851455
    VTX-0087278
    VTX-0860277
    VTX-0079370"""
    keep_vtx_no_translated = df['vtx_id'].isin([v.strip() for v in keep_vtx_no_translated.split('\n')])
    df.loc[keep_vtx_no_translated, 'screening_phase'] = 'TEMPORARY_KEEP'

    df = filter_riboseq(df) # Fine?? doesn't rely on any data/cache info, but could perform this in cache update too

    df = add_temp_gwas(df)

    df = add_temp_rarevar(df)

    df = add_temp_metaorf_score(df)

    df = add_temp_uniprot_annotation(df)

    # df = add_temp_feature_info(df) # Columns added directly to the sequence_features_strings.csv file
    feature_df = pd.read_csv(CACHE_DIR.joinpath('protein_data', 'sequence_features_strings.csv'), index_col=0)
    feature_cols = ['nonsignal_seqs', 'DeepTMHMM_prediction', 'DeepTMHMM_length']
    df = df.merge(feature_df[feature_cols], left_index=True, right_index=True, how='left')
    df.index.name = 'vtx_id'
    isoform_cols = ['swissprot_isoform', 'ensembl_isoform', 'refseq_isoform']
    df[isoform_cols] = df[isoform_cols].apply(lambda x: [str(i) if isinstance(i, np.ndarray) else 'None' for i in x], axis=0)
    df['vtx_id'] = df.index
    
    df = add_temp_nonsig_cons_info(df) # generates core output in run_protein_search_tools.py

    df = add_temp_transcript_exact(df)

    df = reorder_table_cols(df)
    df.rename({'secreted': 'secreted_hibit',
               'translated': 'translated_hibit'}, axis=1, inplace=True)
    
    df[['start', 'end']] = df[['start', 'end']].astype(int)

    signal_cols = ['SignalP 4.1_cut', 'SignalP 5b_cut', 'SignalP 6slow_cut', 'Deepsig_cut']
    measured_secreted_or_predicted_secreted = df['secreted_hibit'] | (df[signal_cols] > -1).any(axis=1)
    df = df[(measured_secreted_or_predicted_secreted) | (df['DeepTMHMM_prediction'])]

    return df


def reorder_table_cols(df):
    """
    """
    view_cols = [
        'show_details', 'vtx_id', 'aa_length', 'ucsc_track', 'source', 
        'protein_xrefs', 'gene_xrefs', 'transcripts_exact', 
        'screening_phase_id', 'uniprot_annotation_score', 'MetaORF v1.0 Score',
        'aa', 'nonsignal_seqs', 
        'blastp_subject', 'blastp_hit_description',
        'blastp_align_length', 'blastp_align_identity', 'nonsig_blastp_align_identity',
        'tblastn_hit_id', 'tblastn_description',
        'tblastn_align_length', 'tblastn_align_identity', 'nonsig_tblastn_align_identity',
        'Deepsig_score', 'SignalP 6slow_score', 'SignalP 5b_score', 'SignalP 4.1_score', 
        'Deepsig_cut', 'SignalP 6slow_cut', 'SignalP 5b_cut', 'SignalP 4.1_cut', 
        'DeepTMHMM', 'DeepTMHMM_prediction', 'DeepTMHMM_length',
        'translated_mean', 'secreted_mean', 'secreted', 'translated', 
        'phylocsf_58m_avg', 'phylocsf_58m_max', 'ESMFold plddt 90th percentile',
        'swissprot_isoform', 'ensembl_isoform', 'refseq_isoform', 
        'spdis_ot', 'consequences_ot', 'trait_ot', 'coding_variant_ot',
        'spdis_gb', 'consequences_gb', 'trait_gb', 'coding_variant_gb',
        'chr', 'start', 'end', 'Ribo-Seq sORF'
        ]
    
    return df[view_cols]


def add_temp_gwas(df):
    """
    """
    session = base.Session()

    opentargets_df = pd.read_parquet(DATA_DIR.joinpath('genetics', 'opentargets_variant_summary_table.parq'))
    opentargets_df = opentargets_df[~opentargets_df['orf_idx_str'].str.startswith('ENSG')]

    group_cols = ['orf_idx_str', 'consequences', 'coding_change', 'ot_spdis', 'ot_pvals', 'ot_betas', 'ot_traits', ]
    grouped_opentargets_df = opentargets_df[group_cols].groupby('orf_idx_str').aggregate(list)
    vtx_map = dict(session.query(Orf.orf_idx_str, Orf.vtx_id).filter(Orf.orf_idx_str.in_(grouped_opentargets_df.index)).all())
    grouped_opentargets_df.reset_index(inplace=True)
    grouped_opentargets_df['vtx_id'] = grouped_opentargets_df.apply(lambda x: vtx_map[x['orf_idx_str']] if x['orf_idx_str'] in vtx_map.keys() else '', axis=1)
    grouped_opentargets_df.set_index('vtx_id', inplace=True)
    df = df.merge(grouped_opentargets_df, left_index=True, right_index=True, how='left')
    df['coding_variant_ot'] = df.apply(lambda x: True in x['coding_change'] if isinstance(x['coding_change'], list) else False, axis=1)
    df.set_index('vtx_id', inplace=True)
    session.close()

    return df


def add_temp_rarevar(df):
    """
    """
    session = base.Session()

    genebass_df = pd.read_parquet(DATA_DIR.joinpath('genetics', 'genebass_variant_summary_table_20240612.parq'))
    genebass_df = genebass_df[~genebass_df['orf_idx_str'].str.startswith('ENSG')]

    group_cols = ['orf_idx_str', 'consequences', 'coding_change', 'gb_spdis', 'gb_pvals', 'gb_betas', 'gb_traits', 'gb_coding_descriptions']

    grouped_genebass_df = genebass_df[group_cols].groupby('orf_idx_str').aggregate(list)
    vtx_map = dict(session.query(Orf.orf_idx_str, Orf.vtx_id).filter(Orf.orf_idx_str.in_(grouped_genebass_df.index)).all())
    grouped_genebass_df.reset_index(inplace=True)
    grouped_genebass_df['vtx_id'] = grouped_genebass_df.apply(lambda x: vtx_map[x['orf_idx_str']] if x['orf_idx_str'] in vtx_map.keys() else '', axis=1)
    grouped_genebass_df = grouped_genebass_df[~pd.isna(grouped_genebass_df['vtx_id'])]
    grouped_genebass_df.set_index('vtx_id', inplace=True)
    df = df.merge(grouped_genebass_df, suffixes=('_ot', '_gb'), left_index=True, right_index=True, how='left')
    df['coding_variant_gb'] = df.apply(lambda x: True in x['coding_change_gb'] if isinstance(x['coding_change_gb'], list) else False, axis=1)
    #df.set_index('vtx_id', inplace=True)

    df['spdis_ot'] = df['ot_spdis'].apply(lambda x: ';'.join(x) if isinstance(x, list) else x).astype(str)
    df['spdis_gb'] = df['gb_spdis'].apply(lambda x: ';'.join(x) if isinstance(x, list) else x).astype(str)

    df['consequences_ot'] = df['consequences_ot'].apply(lambda x: ';'.join(x) if isinstance(x, list) else x)
    df['consequences_gb'] = df['consequences_gb'].apply(lambda x: ';'.join(x) if isinstance(x, list) else x)

    df['trait_ot'] = df['ot_traits'].apply(lambda x: ';'.join(x) if isinstance(x, list) else x)
    #df['trait_gb'] = df['gb_traits'].apply(lambda x: ';'.join(x) if isinstance(x, list) else x)

    df['pval_ot'] = df['ot_pvals'].apply(lambda x: ';'.join([str(y) for y in x]) if isinstance(x, list) else x)
    df['pval_gb'] = df['gb_pvals'].apply(lambda x: ';'.join([str(y) for y in x]) if isinstance(x, list) else x)
    
    df['trait_gb'] = df['gb_coding_descriptions'].apply(lambda x: ';'.join([str(y) for y in x]) if isinstance(x, list) else x)

    session.close()

    return df


def add_temp_hibit(df):
    """
    """
    df['nucl_no_stop'] = df.apply(lambda x: x.nucl[:-3], axis=1)

    ph1to8_hit_df = pd.read_csv(DATA_DIR.joinpath('interim_phase1to8_all_20231012.csv'))
    ph9_hit_df = pd.read_excel(DATA_DIR.joinpath('Ph9 hit list final.xlsx'), sheet_name='Secreted')

    secreted_ph1to8_vtx = list(ph1to8_hit_df[ph1to8_hit_df['secreted']]['vtx_id'].values)
    secreted_ph9_vtx = list(ph9_hit_df.merge(df, left_on='endogenous_nt_seq', right_on='nucl_no_stop')['vtx_id'].values)

    secreted_vtx = set(secreted_ph1to8_vtx + secreted_ph9_vtx)

    df['secreted'] = df.apply(lambda x: x.vtx_id in secreted_vtx or x['secreted'], axis=1)
    hibit = pd.read_excel(DATA_DIR.joinpath('phase12_standardized.xlsx'), sheet_name = 'Primary HTS')
    hibit_confirmation = pd.read_excel(DATA_DIR.joinpath('phase12_standardized.xlsx'), sheet_name = 'Confirmation HTS', index_col=0)

    for ix, row in hibit.iterrows():
        vtx_id = row['VTX ID']
        if vtx_id not in df.index:
            print(f'{vtx_id} not found in df please rerun etl.')
        else:
            for key in ['trans1', 'trans2', 'trans3', 'sec1', 'sec2', 'sec3']:
                df.at[vtx_id, key] = row[key]
            df.at[vtx_id, 'secreted_mean'] = np.nanmean(df.loc[vtx_id][['trans1', 'trans2', 'trans3']])
            df.at[vtx_id, 'translated_mean'] = np.nanmean(df.loc[vtx_id][['sec1', 'sec2', 'sec3']])
            df.at[vtx_id, 'secreted'] = True if hibit_confirmation.loc[vtx_id, 'secretion hit '] == 'Yes' else False
            df.at[vtx_id, 'translated'] = True if hibit_confirmation.loc[vtx_id, 'translation hit'] == 'Yes' else False
    return df


def add_temp_source(df):
    """
    """
    session = base.Session()
    source_df = pd.DataFrame(session.query(Orf.vtx_id, Dataset.name).\
                                     join(OrfXref, OrfXref.orf_id == Orf.id).\
                                     join(Dataset, Dataset.id == OrfXref.xref_dataset_id).\
                                     filter(Orf.vtx_id.in_(list(df['vtx_id']))).all())
    
    source_df = source_df.groupby('vtx_id').aggregate(list)
    source_df['source'] = source_df.apply(lambda x: list(set(x['name'])), axis=1)
    df.drop(columns='source', inplace=True)
    df = df.merge(source_df[['source']], left_index=True, right_index=True, how='left')
    session.close()

    df['source'] = df.apply(lambda x: x.source if isinstance(x.source, list) else [], axis=1)

    return df


def add_temp_uniprot_annotation(df):
    """
    """
    session = base.Session()
    seqs = list(df['aa'].values)
    aa_prot = {p.aa_seq: p for p in session.query(Protein).filter(Protein.aa_seq.in_(seqs)).all()}
    
    annotation_scores = []
    for i, row in df.iterrows():
        if row.aa in aa_prot.keys():
            prot = aa_prot[row.aa]
            try:
                annotation_scores.append(prot.attrs['uniprot_annotation_score'])
            except:
                annotation_scores.append(-1.0)
        else:
            annotation_scores.append(-1.0)
    
    df['uniprot_annotation_score'] = annotation_scores

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


def get_tx_synonyms(tx):
    """
    """
    tx_synonyms = []
    if tx.ensembl_id != '':
        tx_synonyms.append(tx.ensembl_id)
    if tx.refseq_id != '':
        tx_synonyms.append(tx.refseq_id)
    if tx.chess_id != '':
        tx_synonyms.append(tx.chess_id)

    return tx_synonyms


def add_temp_transcript_exact(df):
    """
    """
    session = base.Session()

    orfs_query = session.query(Orf).filter(Orf.vtx_id.in_(df['vtx_id']))
    vtx_orf_id_map = {o.id: o.vtx_id for o in orfs_query.all()}

    tx_orf_df = pd.DataFrame(session.query(base.TranscriptOrf.orf_id, base.TranscriptOrf.transcript_id).filter(base.TranscriptOrf.orf_id.in_(vtx_orf_id_map.keys())).all())
    transcripts = session.query(base.Transcript).filter(base.Transcript.id.in_(tx_orf_df['transcript_id'])).all()
    tx_synonyms = {tx.id: get_tx_synonyms(tx) for tx in transcripts}

    tx_orf_df['transcripts_exact'] = tx_orf_df.apply(lambda x: tx_synonyms[x.transcript_id], axis=1)
    tx_orf_df['vtx_id'] = tx_orf_df.apply(lambda x: vtx_orf_id_map[x.orf_id], axis=1)

    tx_group_df = tx_orf_df[['vtx_id', 'transcripts_exact']].groupby('vtx_id').aggregate(list)
    tx_group_df['transcripts_exact'] = tx_group_df.apply(lambda x: list(set([subitem for item in x.transcripts_exact for subitem in item])), axis=1)

    df.drop(columns='transcripts_exact', inplace=True)

    df = df.merge(tx_group_df, left_index=True, right_index=True, how='left')

    return df


def add_temp_metaorf_score(df):
    """
    """
    session = base.Session()

    metaorf_score_df = pd.DataFrame(session.query(OrfXref.orf_id, Orf.vtx_id, OrfXref.xref).\
                    join(Orf, Orf.id == OrfXref.orf_id).\
                    filter(Orf.vtx_id.in_(df.index)).\
                    filter(and_(OrfXref.xref_dataset_id == 145, OrfXref.type == 'score')).all())
    metaorf_score_df = metaorf_score_df.groupby('vtx_id').first()
    metaorf_score_df['MetaORF v1.0 Score'] = metaorf_score_df['xref'].astype(float)
    df = df.merge(metaorf_score_df, left_index=True, right_index=True, how='left')

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
        (df['source'].apply(lambda x: 'velia_phase9_orfrater' in x)) | \
        (df['source'].apply(lambda x: 'velia_phase9_Olsen' in x)) | \
        (df['source'].apply(lambda x: 'velia_phase9_Li et al VSMC' in x)) | \
        (df['source'].apply(lambda x: 'velia_phase7_Ribo-seq_PBMC_LPS_R848' in x)) | \
        (df['source'].apply(lambda x: 'velia_phase10_riboseq_230114' in x)) | \
        (df['source'].apply(lambda x: 'velia_phase11_riboseq_240214' in x)) | \
        (df['source'].apply(lambda x: 'swissprot' in x)) | \
        (df['source'].apply(lambda x: 'MetaORF v1.0' in x)) | \
        (df['source'].apply(lambda x: 'ENSEMBL' in x)) | \
        (df['source'].apply(lambda x: 'BestRefSeq' in x)) | \
        (df['source'].apply(lambda x: 'velia_phase5_uniprot-tremble' in x)) | \
        (df['screening_phase'] == 'Not Screened') | \
        (df['orf_xrefs'].astype(str).str.contains('RibORF')) | \
        (df['protein_xrefs'].astype(str).str.contains('RibORF')) | \
        (df['screening_phase'] == 'TEMPORARY_KEEP') | \
        (df['MS or Riboseq'])
    )
    
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


@st.cache_data()
def load_hmmer_results(vtx_id, hmmer_loc=HMMER_S3_LOC):
    """
    """
    try:
        return pd.read_parquet(f'{hmmer_loc}{vtx_id}.parq')
    except FileNotFoundError as e:
        return 
