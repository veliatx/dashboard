import json
import jsonlines
import os
import pickle
import gzip
import py3Dmol

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px
import streamlit as st
import streamlit.components.v1 as components
import streamlit_scrollable_textbox as stx

from collections import defaultdict
from streamlit_plotly_events import plotly_events
from tqdm import tqdm

from streamlit_echarts import st_echarts
from scipy.cluster.hierarchy import linkage, leaves_list

from dashboard.util import filter_dataframe, convert_list_string
from dashboard.etl.sorf_query import load_jsonlines_table
from dashboard import plotting, description
from dashboard import tabs

import pyarrow.parquet as pq

CACHE_DIR = '../cache'
TPM_DESEQ2_FACTOR = 80


@st.cache_data()
def load_autoimmune_atlas():
#     exp = pd.read_parquet("s3://velia-athena-dev/expression_atlas_v1_expression_per_group.parq")
#     exp = exp.groupby([x.split('.')[0] if x.startswith('ENST') else x for x in exp.index]).mean()
    meta = pd.read_parquet("s3://velia-athena-dev/expression_atlas_v1_metadata.parq")
    meta['dashboard_group'] = meta['atlas_group']
    return meta

# def load_autoimmune_expression_for_transcript(transcript_id):
#     log10padj_threshold = -2
#     minimum_expression = 2
#     con = sqlite3.connect('../data/autoimmune_expression_atlas_v1.db')
#     query = """SELECT *
#     FROM transcript_de
#     WHERE transcript_de.transcript_id = '{0}'
#     AND transcript_de.log10_padj <= {1}
#     AND (transcript_de.case_mean >= {2} OR transcript_de.control_mean >= {2})
#     """.format(transcript_id, log10padj_threshold, minimum_expression)
#     return pd.read_sql(query, con)

@st.cache_data()
def load_sorf_df_conformed():
    df = pd.read_parquet(os.path.join(CACHE_DIR, 'sorf_df.parq'))
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


# # @st.cache_data()
# def load_de_results(transcripts):
#     cache_filez = os.listdir(CACHE_DIR)
#     temp_dict = {}
#     for f in cache_filez:
#         if f.endswith('_de.parq') and not (f=='expression_de.parq'):
#             df = pd.read_parquet(os.path.join(CACHE_DIR, f))
#             df['transcript'] = df.apply(lambda x: x.name.split('.')[0], axis=1)
#             df = df[df['transcript'].isin(transcripts)].copy()
#             temp_dict[f.split('_')[0]] = df
            
#     de_tables_dict = defaultdict(dict)
#     for c, df in tqdm(temp_dict.items()):
#         for row in df.itertuples():
#             de_tables_dict[row[0]][c] = {'Cancer Average': row._7/TPM_DESEQ2_FACTOR, 'GTEx Average': row._8/TPM_DESEQ2_FACTOR, 
#                                          'log2FC': row.log2FoldChange, 'padj': row.padj}
#     for t, d in de_tables_dict.items():
#         de_tables_dict[t] = pd.DataFrame(d).T
#     tcga_gtex_tissue_metadata = pd.read_parquet(os.path.join(CACHE_DIR, 'gtex_tcga_pairs.parq'))
#     tcga_gtex_tissue_metadata = tcga_gtex_tissue_metadata.drop_duplicates(['TCGA Cancer Type', 'GTEx Tissue Type']).copy()
#     tcga_gtex_tissue_metadata.index = tcga_gtex_tissue_metadata['TCGA Cancer Type']
#     return de_tables_dict, tcga_gtex_tissue_metadata

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