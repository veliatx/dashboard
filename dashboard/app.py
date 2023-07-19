import json
import jsonlines
import os
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
from dashboard.util import filter_dataframe, convert_list_string
from dashboard.etl.sorf_query import load_jsonlines_table

from veliadb import base
from veliadb.base import Orf, Protein, ProteinXref
import gzip

CACHE_DIR = '../cache'
TPM_DESEQ2_FACTOR = 80


@st.cache_data()
def load_protein_feature_string_representations():
    df = pd.read_csv(os.path.join(CACHE_DIR, 'protein_data', 'sequence_features_strings.csv'), index_col=0).T
    return df


@st.cache_data()
def load_sorf_df():
    """
    """
    session = base.Session()

    sorf_df = pd.read_excel('../data/interim_phase1to6_secreted_hits_20230330.xlsx')
    secreted_ph7_hit_ids = pd.read_csv('../data/phase_7_secreted_hits.csv', index_col=0)

    secreted_vtx_ids = set(sorf_df['vtx_id']).union(set(secreted_ph7_hit_ids.index))

    sorf_df = load_jsonlines_table(os.path.join(CACHE_DIR, 'sorf_table.jsonlines'), index_col='vtx_id')

    # removes any sORFs with multiple VTX IDs (e.g. multi-mappers to the genome)
    sorf_df = sorf_df[sorf_df['vtx_id'].isin(secreted_vtx_ids)]
    sorf_df = sorf_df[~sorf_df['vtx_id'].str.contains('\|')]

    sorf_df['index_copy'] = sorf_df.index

    sorf_df['show_details'] = False
    sorf_df['orf_xrefs'] = sorf_df.apply(lambda x: x.orf_xrefs.split(';'), axis=1)
    sorf_df['source'] = sorf_df.apply(lambda x: x.source.split(';'), axis=1)

    cols = list(sorf_df.columns)
    cols.insert(0, cols.pop(cols.index('show_details')))

    phase_ids = []
    phase_entries = []
    protein_xrefs = []
    for i, row in sorf_df.iterrows():
        protein_xrefs.append([str(px.xref) for px in \
                                session.query(ProteinXref)\
                                        .join(Protein, Protein.id == ProteinXref.protein_id)\
                                        .filter(Protein.aa_seq == row.aa).all()])
        if row.velia_id.startswith('Phase'):
            phase_ids.append(row.velia_id)
            phase_entries.append(f'Phase {row.velia_id[6]}')
        elif row.secondary_id.startswith('Phase'):
            phase_ids.append(row.secondary_id)
            phase_entries.append(f'Phase {row.secondary_id[6]}')
        elif any(['Phase' in xref for xref in row.orf_xrefs]):
            phase_id = [x for x in row.orf_xrefs if x.startswith('Phase')][0]
            phase_ids.append(phase_id)
            phase_entries.append(f'Phase {phase_id[6]}')
        elif row.velia_id.startswith('smORF'):
            phase_ids.append(row.velia_id)
            phase_entries.append('Phase 1')
        else:
            orf = session.query(Orf).filter(Orf.id == int(row.vtx_id.split('-')[1])).one()
            phase_ids.append(orf.velia_id)
            if orf.velia_id.startswith('Phase'):
                phase_entries.append(f'Phase {orf.velia_id[6]}')
            else:
                phase_entries.append('-1')

    sorf_df['screening_phase_id'] = phase_ids
    sorf_df['screening_phase'] = phase_entries
    sorf_df['protein_xrefs'] = protein_xrefs

    sorf_df = sorf_df[sorf_df['screening_phase'] != '-1']

    cols = ['show_details', 'vtx_id', 'screening_phase_id', 'screening_phase', 'ucsc_track', 
            'source', 'orf_xrefs', 'protein_xrefs', 'gene_xrefs', 'transcript_xrefs',  
            'transcripts_exact', 'transcripts_overlapping', 'aa', 'nucl', 
            'index_copy', 'genscript_id', 'chr', 'strand', 'start', 'end', 
            'chrom_starts', 'block_sizes', 'phases',]

    sorf_df = sorf_df[cols]
    session.close()
    _, blastp_table = load_mouse_blastp_results()
    sorf_df = sorf_df.merge(pd.DataFrame(blastp_table).T, left_index=True, right_index=True, how='left')
    protein_scores = pd.read_csv(os.path.join(CACHE_DIR, 'protein_data', 'sequence_features_scores.csv'), index_col=0)

    sorf_df = sorf_df.merge(protein_scores[['Deepsig', 'SignalP 6slow', 'SignalP 5b', 'SignalP 4.1']],
                  left_index=True, right_index=True, how='left')
    protein_strings = pd.read_csv(os.path.join(CACHE_DIR, 'protein_data', 'sequence_features_strings.csv'), index_col=0)
    sorf_df.insert(4, 'aa_length', protein_strings['Sequence'].str.len())
    protein_cutsite = protein_strings.apply(lambda x: x.str.find('SO')+1).replace(0, -1).drop('Sequence', axis=1)
    sorf_df = sorf_df.merge(protein_cutsite,
                  left_index=True, right_index=True, how='left', suffixes=('_score', '_cut'))
    id_data = pd.read_csv('s3://velia-data-dev/VDC_001_screening_collections/all_phases/interim_phase1to7_non-sigp_20230718.csv')
    id_data.index = id_data['vtx_id']
    sorf_df = sorf_df.merge(id_data[['trans1',
           'trans2', 'trans3', 'sec1', 'sec2', 'sec3', 'translated_mean',
           'secreted_mean', 'translated', 'swissprot_isoform', 'ensembl_isoform',
           'refseq_isoform', 'phylocsf_58m_avg', 'phylocsf_58m_max', 'phylocsf_58m_min',
           'phylocsf_vals']], left_index=True, right_index=True, how='left')
    
    with open('../data/all_secreted_phase1to7.txt', 'r') as f:
        secreted_ids = [i.strip() for i in f.readlines()]
    sorf_df.insert(int(np.where(sorf_df.columns=='translated')[0][0]), 'secreted', [True if i in secreted_ids else False for i in sorf_df.index])

    
    return sorf_df


@st.cache_data()
def load_kibby_results(sorf_table):
    kibby = pd.read_csv(os.path.join(CACHE_DIR, 'protein_data', 'kibby.out'), index_col=0)
    
    kibby['conservation'] = kibby.conservation.apply(lambda x: list(map(float, x.strip().split(' '))))
    subdf = sorf_table.loc[sorf_table.index.intersection(kibby.index)]
    kibby = kibby.loc[subdf.index]
    kibby.index = sorf_table.loc[kibby.index, 'vtx_id']
    return kibby


@st.cache_data()
def load_de_results(transcripts):
    cache_filez = os.listdir(CACHE_DIR)
    temp_dict = {}
    for f in cache_filez:
        if f.endswith('_de.parq') and not (f=='expression_de.parq'):
            df = pd.read_parquet(os.path.join(CACHE_DIR, f))
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
    tcga_gtex_tissue_metadata = pd.read_parquet(os.path.join(CACHE_DIR, 'gtex_tcga_pairs.parq'))
    tcga_gtex_tissue_metadata = tcga_gtex_tissue_metadata.drop_duplicates(['TCGA Cancer Type', 'GTEx Tissue Type']).copy()
    tcga_gtex_tissue_metadata.index = tcga_gtex_tissue_metadata['TCGA Cancer Type']
    return de_tables_dict, tcga_gtex_tissue_metadata


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
def load_xena_tcga_gtex_target(vtx_id_to_transcripts):
    # Expression is saved as TPM + 0.001 (NOT LOGGED)
    xena_expression = pd.read_parquet(os.path.join(CACHE_DIR, 'xena.parq'))
    xena_metadata = xena_expression[xena_expression.columns[:6]]
    xena_expression = xena_expression[xena_expression.columns[6:]]

    # vtx_id_to_transcripts = load_jsonlines_table(os.path.join(CACHE_DIR, 'sorf_table.jsonlines'), index_col='vtx')
    all_transcripts = [i for i in np.concatenate((*vtx_id_to_transcripts['transcripts_exact'],
                                                  *vtx_id_to_transcripts['transcripts_overlapping']))
                       if i.startswith('ENST')]
    de_tables_dict, de_metadata = load_de_results(all_transcripts)
    de_tables_dict = {k.split('.')[0]:v for k, v in de_tables_dict.items()}
    # Sum expression over each VTX            
    xena_vtx_sums = xena_expression.T.copy()
    xena_vtx_sums = xena_vtx_sums.loc[xena_vtx_sums.index.intersection(all_transcripts)]
    xena_exact_vtx_sums = {}
    xena_overlapping_vtx_sums = {}
    transcripts_in_xena = xena_expression.columns
    for vtx_id, row in vtx_id_to_transcripts.iterrows():
        if len(transcripts_in_xena.intersection(row['transcripts_exact'])) > 0:
            xena_exact_vtx_sums[vtx_id] = xena_expression[transcripts_in_xena.intersection(row['transcripts_exact'])].sum(axis=1)
        if len(transcripts_in_xena.intersection(row['transcripts_overlapping'])) > 0:
            xena_overlapping_vtx_sums[vtx_id] = xena_expression[transcripts_in_xena.intersection(row['transcripts_overlapping'])].sum(axis=1)
    xena_exact_vtx_sums = pd.DataFrame(xena_exact_vtx_sums)
    xena_overlapping_vtx_sums = pd.DataFrame(xena_overlapping_vtx_sums)
    xena_exact_heatmap_data = process_sums_dataframe_to_heatmap(xena_exact_vtx_sums, xena_metadata)
    xena_overlapping_heatmap_data = process_sums_dataframe_to_heatmap(xena_overlapping_vtx_sums, xena_metadata)
    return xena_metadata, xena_expression, vtx_id_to_transcripts, xena_exact_heatmap_data, xena_overlapping_heatmap_data, de_tables_dict, de_metadata


@st.cache_data()
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


@st.cache_data()
def load_mouse_blastp_results():
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


@st.cache_data()
def convert_df(df):
    return df.to_csv().encode('utf-8')


def sorf_details(sorf_df):
    st.title('sORF Table')
    st.write('Table contains library of secreted sORFs.')

    df = filter_dataframe(sorf_df)

    if 'data_editor_prev' in st.session_state.keys():
        curr_rows = st.session_state['data_editor']['edited_rows']
        prev_rows = st.session_state['data_editor_prev']['edited_rows']

        update_rows = {key: curr_rows[key] for key in curr_rows if key not in prev_rows}

        if len(update_rows) > 0:
            row_idx = list(update_rows.keys())[0]
            
            df['show_details'] = False
            row_idx_name = df.iloc[row_idx].name
            df.at[row_idx_name, 'show_details'] = True
            st.session_state['curr_vtx_id'] = str(df.loc[row_idx_name]['vtx_id'])

    st.write(f'{df.shape[0]} sORF entries')
    update_df = st.data_editor(
        df,
        column_config={
            'vtx_id': st.column_config.TextColumn(disabled=True),
            'screening_phase_id': st.column_config.TextColumn(disabled=True),
            'screening_phase': st.column_config.TextColumn(disabled=True),
            'genscript_id': st.column_config.TextColumn(disabled=True),
            #'orf_xrefs': st.column_config.TextColumn(disabled=True),
            #'protein_xrefs': st.column_config.TextColumn(disabled=True),
            #'gene_xrefs': st.column_config.TextColumn(disabled=True),

        },
        key='data_editor',
        hide_index=True
    )

    st.download_button(
        label="Download as CSV",
        data=convert_df(update_df),
        file_name='sorf_selection.csv',
        mime='text/csv',
    )

    st.session_state['data_editor_prev'] = st.session_state['data_editor'].copy()

    # Load data
    tcga_data = load_xena_tcga_gtex_target(sorf_df)
    xena_metadata, xena_expression, vtx_id_to_transcripts, _, _, de_tables_dict, de_metadata = tcga_data
    esmfold = load_esmfold()
    blastp_mouse_hits, blastp_data_for_sorf_table = load_mouse_blastp_results()
    kibby = load_kibby_results(sorf_df)
    protein_features_df = load_protein_feature_string_representations()
    phylocsf_dataframe = load_phylocsf_data()
    xena_overlap = []

    if 'curr_vtx_id' in st.session_state.keys():

        vtx_id = st.session_state['curr_vtx_id']
        selected_row = df[df['vtx_id'] == st.session_state['curr_vtx_id']]

        st.divider()
        st.header('sORF Details')
        st.dataframe(selected_row[['vtx_id', 'screening_phase_id', 'orf_xrefs', 'protein_xrefs', 'gene_xrefs']])

        #link = ucsc_link(selected_row['chromosome'], selected_row['start'], selected_row['end'])
        #st.write(f"UCSC Browser [link]({link})")
        selected_transcripts_exact = vtx_id_to_transcripts.loc[vtx_id, 'transcripts_exact']
        selected_transcripts_overlapping = vtx_id_to_transcripts.loc[vtx_id, 'transcripts_overlapping']
        selected_transcripts = np.concatenate([selected_transcripts_exact, selected_transcripts_overlapping])        
        xena_overlap = xena_expression.columns.intersection(selected_transcripts)
        value = None

        with st.expander("Transcription Data", expanded=True):
            col1, col2 = st.columns(2)

            with col1:
                title = f'Transcript Specific Expression - {vtx_id}'
                option, events = plotting.expression_heatmap_plot(vtx_id, vtx_id_to_transcripts, xena_expression, xena_metadata, title)
                
                if option:
                    value = st_echarts(option, height="1000px", events=events, renderer='svg')
                else:
                    st.write('No transcripts in TCGA/GTEx/TARGET found containing this sORF')

            with col2:

                if len(xena_overlap) > 0:
                    if value:
                        if value.startswith('**'):
                            selected_transcript = [value[2:]]
                        else:
                            selected_transcript = [value]

                    elif selected_transcripts.shape[0]:
                        selected_transcript = [xena_overlap[0]]

                    chart_title = f'Differential Expression - {selected_transcript[0]}'

                    de_exact_echarts_options_b = plotting.plot_transcripts_differential_expression_barplot(selected_transcript, 
                                                                                                           de_tables_dict, de_metadata,
                                                                                                           chart_title)
                                                                                                    
                    st_echarts(options=de_exact_echarts_options_b, key='b', height='900px', width = '600px', renderer='svg')
                
            if (len(xena_overlap)>0) and value:
                st.write(value)

                xena_vtx_exp_df = xena_metadata.merge(xena_expression, left_index=True, right_index=True)
                fig = plotting.expression_vtx_boxplot(value, xena_vtx_exp_df)
                st.plotly_chart(fig, use_container_width=True)
                
        with st.expander("Protein Structure and Function", expanded=True):
            col3, col4 = st.columns(2)

            with col3:

                # Load esmfold data for selected sORF
                sorf_aa_seq = sorf_df[sorf_df['vtx_id']==vtx_id]['aa'].iloc[0]
                if sorf_aa_seq[-1] == '*':
                    sorf_aa_seq = sorf_aa_seq[:-1]

                plddt = esmfold[sorf_aa_seq]['plddt']
                # Plot plDDT, Phylocsf, and kibby
                achart = plotting.plot_sequence_line_plots_altair(vtx_id, sorf_aa_seq, phylocsf_dataframe, kibby, esmfold)
                col3.altair_chart(achart, use_container_width=False)
                        
            with col4:
                structure = esmfold[sorf_aa_seq]['pdb']
                modified_structure_colors = plotting.color_protein_terminal_ends(sorf_aa_seq, structure)
                view = py3Dmol.view(js='https://3dmol.org/build/3Dmol.js',)
                view.addModel(modified_structure_colors, 'pdb')
                view.setStyle({'cartoon': {'colorscheme': {'prop':'b','gradient': 'roygb','min':50,'max':90}}})
                view.zoomTo()
                st.header('sORF ESMfold', help="Red - Low Confidence  \nBlue - High Confidence  \nConfidence is based on plDDT score from ESMFold  \nN-term is blue and C-term is red")
                # st.text('Red - Low confidence')
                # st.text('Blue - High confidence')
                # st.text('N-term is blue and C-term is red')
                components.html(view._make_html(), height=500, width=600)
                
            f = protein_features_df[vtx_id]
            imdf = plotting.format_protein_feature_strings_for_altair_heatmap(f)
            altair_signal_features_fig = plotting.altair_protein_features_plot(imdf)
            st.header('AA Feature Predictions', help=description.amino_acid_features_hover_text)
            st.altair_chart(altair_signal_features_fig, use_container_width=True)

        
        with st.expander("BLASTp results", expanded=True):
            # Blastp Mouse
            # primary_id = sorf_df[sorf_df['vtx_id'] == vtx_id].iloc[0]['screening_phase_id']
            blastp_results_selected_sorf = blastp_mouse_hits[vtx_id]
            if len(blastp_results_selected_sorf) == 0:
                long_text = "No alignments with mouse found." #st.write('No alignments with mouse found.')
            else:
                long_text = ""
                for h in blastp_results_selected_sorf:
                    hit_text = f"Match IDs: {h['hit_ids']}  \nAlign Stats: Score - {h['score']}, Length - {h['align_len']}  \n"
                    long_text+=hit_text
                    long_text+= h['alignment'] + '  \n  \n'
    
            stx.scrollableTextbox(long_text, height = 300, fontFamily='Courier')


def sorf_transcriptome_atlas(sorf_df):
    st.title("sORF Transcriptome Atlas")
    tcga_data = load_xena_tcga_gtex_target(sorf_df)
    xena_metadata, xena_expression, vtx_id_to_transcripts, xena_exact_heatmap_data, xena_overlapping_heatmap_data, de_tables_dict, de_metadata = tcga_data
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
           

def selector(sorf_df):
    st.title("Select sORFs Interactively")
    st.session_state['x_feature'] = st.selectbox("Select X Feature", options = sorf_df.columns, index=11)
    st.session_state['y_feature'] = st.selectbox("Select Y Feature", options = sorf_df.columns, index=12)
    
    fig = px.scatter(sorf_df, x=st.session_state['x_feature'], y=st.session_state['y_feature'])
    selected_points = plotly_events(fig)
    st.write(selected_points)


def genome_browser():
    """
    """
    components.iframe("http://10.65.25.231:8080/velia_collections.html", height=1200, scrolling=True)


def main():
    
    st.set_page_config(layout="wide")
    st.markdown(""" <style>iframe[title="streamlit_echarts.st_echarts"]{ height: 1000px !important } """, unsafe_allow_html=True)
    with open("style.css") as f:
        st.markdown(f"<style>{f.read()}</style>", unsafe_allow_html=True)
    
    # Define a dictionary with the page names and corresponding functions
    pages = {
        "sORF Details": sorf_details,
        "sORF Transcriptome Atlas": sorf_transcriptome_atlas,
        "sORF Genome Browser": genome_browser,
        # "sORF Selector": selector,
    }
    sorf_df = load_sorf_df()
  
    tab1, tab2, tab3 = st.tabs(list(pages.keys()))

    with tab1:
        sorf_details(sorf_df)

    with tab2:
        sorf_transcriptome_atlas(sorf_df)

    with tab3:
        genome_browser()
        
    # with tab4:
    #     selector(sorf_df)


if __name__ == "__main__":
    main()
