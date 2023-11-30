import json
import jsonlines
import os
import pickle
import gzip

import pandas as pd
import streamlit as st

from collections import defaultdict
from streamlit_plotly_events import plotly_events

from streamlit_echarts import st_echarts
from scipy.cluster.hierarchy import linkage, leaves_list

from dashboard.util import filter_dataframe, convert_list_string
from dashboard.etl.sorf_query import load_jsonlines_table
from dashboard import tabs

import pyarrow.parquet as pq
from dashboard.etl import CACHE_DIR, TPM_DESEQ2_FACTOR


@st.cache_data()
def load_autoimmune_atlas():
    meta = pd.read_parquet("s3://velia-athena-dev/expression_atlas_v1_metadata.parq")
    meta['dashboard_group'] = meta['atlas_group']
    return meta

@st.cache_data()
def load_sorf_df_conformed():
    df = pd.read_parquet(os.path.join(CACHE_DIR, 'sorf_df.parq'))
    
    # TODO remove this as temp addition
    ribo_df = pd.read_excel('../data/Secreted_mP_Riboseq_SAF.xlsx')
    ribo_vtx = set(ribo_df[ribo_df['manual_check'] == 1]['vtx_id'])
    ccle_df = pd.read_excel('../data/SummaryIdentification_CCLE_strongerConfidence.xlsx', index_col=0)
    gtex_df = pd.read_excel('../data/SummaryIdentification_GTEX_strongerConfidence.xlsx', index_col=0)

    ccle_vtx = set(ccle_df['vtx_id'])
    gtex_vtx = set(gtex_df['vtx_id'])

    ms_vtx = gtex_vtx.union(ccle_vtx)
    support_vtx = ms_vtx.union(ribo_vtx)

    df['manual_riboseq'] = df.apply(lambda x: True if x.vtx_id in ribo_vtx else False, axis=1)
    df['MS_evidence'] = df.apply(lambda x: True if x.vtx_id in ms_vtx else False, axis=1)
    df['MS or Riboseq'] = df.apply(lambda x: True if x.vtx_id in support_vtx else False, axis=1)

    feat_df = pd.read_csv('../data/interim_phase1to8_all_20231012.csv')
    feat_df.set_index('vtx_id', inplace=True)
    df = df.merge(feat_df[['tblastn_hit_id', 'tblastn_description',
                           'tblastn_score', 'tblastn_query_coverage', 'tblastn_align_length',
                           'tblastn_align_identity', 'tblastn_gaps', 'tblastn_evalue']], how='left', left_index=True, right_index=True)
    
    df.index.name = 'vtx_id'
    df['vtx_id'] = df.index

    return df

@st.cache_data()
def load_protein_feature_string_representations():
    df = pd.read_csv(os.path.join(CACHE_DIR, 'protein_data', 'sequence_features_strings.csv'), index_col=0).T
    return df

@st.cache_data()
def load_kibby_results(sorf_table):
    kibby = pd.read_csv(os.path.join(CACHE_DIR, 'protein_data', 'kibby.out'), index_col=0)
    
    kibby['conservation'] = kibby.conservation.apply(lambda x: list(map(float, x.strip().split(' '))))
    subdf = sorf_table.loc[sorf_table.index.intersection(kibby.index)]
    kibby = kibby.loc[subdf.index]
    kibby.index = sorf_table.loc[kibby.index, 'vtx_id']
    return kibby

@st.cache_data()
def load_xena_transcripts(transcripts):
    xena_expression = pd.read_parquet(os.path.join(CACHE_DIR, 'xena_app.parq'), columns = transcripts)    
    return xena_expression

@st.cache_data()
def load_xena_metadata():
    xena_metadata = pd.read_parquet(os.path.join(CACHE_DIR, 'xena_metadata.parq'))
    xena_metadata['dashboard_group'] = list(map(lambda x: '-'.join(map(str, x)), xena_metadata[['_primary_site', '_study']].values))
    parq_metadata = pq.read_metadata(os.path.join(CACHE_DIR, 'xena_app.parq'))
    xena_transcripts = set(parq_metadata.schema.names)
    return xena_metadata, xena_transcripts

@st.cache_data()
def load_xena_heatmap_data():
    xena_exact_heatmap_data = pickle.load(open(os.path.join(CACHE_DIR, 'xena_exact_heatmap.pkl'), 'rb'))
    return xena_exact_heatmap_data

@st.cache_data()
def load_esmfold():
    """
    """
    esmfold = {}
    with gzip.open('../data/phase1to7_all_esmfold.jsonlines.gz') as fopen:
        j_reader = jsonlines.Reader(fopen)
        for l in j_reader:
            esmfold[l['sequence']] = l
    return esmfold

@st.cache_data()
def load_mouse_blastp_results(CACHE_DIR = '../cache'):
    hits_per_query = defaultdict(list)
    sorf_table_data = {}
    with open(os.path.join(CACHE_DIR, 'protein_data', 'blastp.results.json'), 'r') as fopen:
        blastp = json.load(fopen)
        blastp = blastp['BlastOutput2']
    for entry in blastp:
        entry = entry['report']['results']
        q = entry['search']['query_title']
        hits = entry['search']['hits']
        if len(hits) == 0:
            # print('No alignments found with mouse')
            pass
        else:
            for h in hits:
                ids = []
                for item in h['description']:
                    ids.append(item['accession'])
                alignment = h['hsps']
                alignment = alignment[0]
                align_str = '  \n'.join([h['description'][0]['title'], alignment['qseq'], alignment['midline'], alignment['hseq']])
                alignment['hit_ids'] = ';'.join(ids)
                alignment['alignment'] = align_str
                hits_per_query[q].append(alignment)
                if isinstance(alignment, dict):
                    best_hit = alignment
                else:
                    best_hit = pd.DataFrame(alignment).sort_values('score', ascending=False).iloc[0]
                best_hit_description = [h for h in hits if h['num'] == best_hit['num']][0]['description'][0]
                sorf_table_data[q] = {'blastp_score': best_hit['score'],
                 'blastp_query_coverage': best_hit['align_len']/len(best_hit['qseq']),
                 'blastp_align_length': best_hit['align_len'],
                 'blastp_gaps': best_hit['gaps'],
                 'blastp_align_identity': best_hit['identity']/best_hit['align_len'],
                'blastp_subject': best_hit_description['id'],
                'blastp_hit_description': best_hit_description['title']
                }
    return hits_per_query, sorf_table_data


@st.cache_data()
def load_phylocsf_data():
    pcsf = pd.read_csv(f"../data/interim_phase1to7_all_phylocsf-vals_20230628.csv", index_col=0)
    pcsf['phylocsf_vals'] = pcsf['phylocsf_vals'].apply(convert_list_string)
    pcsf = pcsf[['phylocsf_58m_avg', 'phylocsf_58m_max',
           'phylocsf_58m_min', 'phylocsf_58m_std', 'phylocsf_vals']]
    return pcsf
