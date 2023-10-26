import json
import jsonlines
import os
import pickle
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

from dashboard import plotting, description
from dashboard.util import filter_dataframe, convert_list_string, convert_df
from dashboard.data_load import *

def tcga_page(sorf_df, tcga_data):
    st.title("sORF Transcriptome Atlas (TCGA)")

    df = filter_dataframe(sorf_df, 'tcga_transcriptome')

    st.write(f'{df.shape[0]} sORF entries')

    xena_metadata, xena_expression, xena_exact_heatmap_data, xena_overlapping_heatmap_data, de_tables_dict, de_metadata = tcga_data
    with st.container():
        col1, col2, col3 = st.columns(3)
        with col1:
            tx_type = st.radio("sORF to transcript mapping",
                               ('Boundary Overlap', 'Exact Overlap'))
            
            if tx_type == 'Boundary Overlap':
                xena_tau_df, xena_vtx_sum_df, xena_vtx_exp_df = xena_overlapping_heatmap_data

            else:
                xena_tau_df, xena_vtx_sum_df, xena_vtx_exp_df = xena_exact_heatmap_data

        with col2:
            values = st.slider(
                'Select Tissue Specificity Tau',
                0.0, 1.0, (.8, 1.0),
                help='Higher values of [Tau](https://academic.oup.com/bib/article/18/2/205/2562739) indicate tissue specific expression of a given sORF')
            tissue_specific_vtx_ids = list(xena_tau_df[xena_tau_df['tau'].between(*values)].index)

        with col3:
            st.write(f'{len(tissue_specific_vtx_ids)}')
            
        option, events = plotting.expression_atlas_heatmap_plot(tissue_specific_vtx_ids, xena_vtx_sum_df)

        display_cols = ['vtx_id', 'screening_phase_id', 'screening_phase', 'orf_xrefs', 'protein_xrefs', 'gene_xrefs', 'transcript_xrefs', 'source']#, 'secreted_mean', 'translated_mean', 'isoform_of']
        value = st_echarts(option, height="1000px", events=events)
        if value:
            st.header('Selected sORF')
            st.dataframe(sorf_df[sorf_df['vtx_id'] == value][display_cols])
            
            fig = plotting.expression_vtx_boxplot(value, xena_vtx_exp_df)
            st.plotly_chart(fig, use_container_width=True)

        df = sorf_df[sorf_df['vtx_id'].isin(tissue_specific_vtx_ids)][display_cols]
        exp_df = xena_vtx_sum_df[tissue_specific_vtx_ids].copy()
        df = df.merge(exp_df.T, left_index=True, right_index=True)
        st.header('Tissue specific sORFs')
        st.dataframe(df)

        st.download_button(
            label="Download as CSV",
            data=convert_df(df),
            file_name='sorf_atlas_selection.csv',
            mime='text/csv',
        )
        
        
# def autoimmune_page(sorf_df, autoimmune_data):
#     st.title("sORF Transcriptome Atlas (Auto Immune)")

#     df = filter_dataframe(sorf_df, 'transcriptome')

#     st.write(f'{df.shape[0]} sORF entries')

#     autoimmune_expression, autoimmune_metadata = autoimmune_data
#     with st.container():
#         col1, col2, col3 = st.columns(3)
#         with col1:
#             tx_type = st.radio("sORF to transcript mapping",
#                                ('Boundary Overlap', 'Exact Overlap'))
            
#             if tx_type == 'Boundary Overlap':
#                 xena_tau_df, xena_vtx_sum_df, xena_vtx_exp_df = xena_overlapping_heatmap_data

#             else:
#                 xena_tau_df, xena_vtx_sum_df, xena_vtx_exp_df = xena_exact_heatmap_data

#         with col2:
#             values = st.slider(
#                 'Select Tissue Specificity Tau',
#                 0.0, 1.0, (.8, 1.0),
#                 help='Higher values of [Tau](https://academic.oup.com/bib/article/18/2/205/2562739) indicate tissue specific expression of a given sORF')
#             tissue_specific_vtx_ids = list(xena_tau_df[xena_tau_df['tau'].between(*values)].index)

#         with col3:
#             st.write(f'{len(tissue_specific_vtx_ids)}')
            
#         option, events = plotting.expression_atlas_heatmap_plot(tissue_specific_vtx_ids, xena_vtx_sum_df)

#         display_cols = ['vtx_id', 'screening_phase_id', 'screening_phase', 'orf_xrefs', 'protein_xrefs', 'gene_xrefs', 'transcript_xrefs', 'source']#, 'secreted_mean', 'translated_mean', 'isoform_of']
#         value = st_echarts(option, height="1000px", events=events)
#         if value:
#             st.header('Selected sORF')
#             st.dataframe(sorf_df[sorf_df['vtx_id'] == value][display_cols])
            
#             fig = plotting.expression_vtx_boxplot(value, xena_vtx_exp_df)
#             st.plotly_chart(fig, use_container_width=True)

#         df = sorf_df[sorf_df['vtx_id'].isin(tissue_specific_vtx_ids)][display_cols]
#         exp_df = xena_vtx_sum_df[tissue_specific_vtx_ids].copy()
#         df = df.merge(exp_df.T, left_index=True, right_index=True)
#         st.header('Tissue specific sORFs')
#         st.dataframe(df)

#         st.download_button(
#             label="Download as CSV",
#             data=convert_df(df),
#             file_name='sorf_atlas_selection.csv',
#             mime='text/csv',
#         )