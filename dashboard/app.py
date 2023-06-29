import streamlit as st
import numpy as np
from streamlit_plotly_events import plotly_events
import plotly.express as px
from plotly import graph_objects as go
import matplotlib.pyplot as plt
import pandas as pd
from collections import defaultdict
from st_aggrid import AgGrid, GridUpdateMode
from st_aggrid.grid_options_builder import GridOptionsBuilder
import streamlit_scrollable_textbox as stx
# import seaborn_altair as salt
import seaborn as sns
import json
import jsonlines
import py3Dmol

from streamlit_echarts import st_echarts
from scipy.cluster.hierarchy import linkage, leaves_list

from dashboard.util import filter_dataframe
from dashboard.etl.sorf_query import load_jsonlines_table
from dashboard import plotting

import os
import altair as alt
from tqdm import tqdm
import streamlit.components.v1 as components

from plotting import plot_structure_plddt, expression_de_to_echarts_data, bar_plot_expression_groups, features_to_int_arrays
CACHE_DIR = '../cache'
TPM_DESEQ2_FACTOR = 80


@st.cache_data()
def load_protein_feature_string_representations():
    df = pd.read_parquet('../data/sequence_features_strings.parq').T
    return df


@st.cache_data()
def load_sorf_excel():
    """
    """
    sorf_excel_df = pd.read_excel('../data/interim_phase1to6_secreted_hits_20230330.xlsx')
    
    # removes any sORFs with multiple VTX IDs (e.g. multi-mappers to the genome)
    sorf_excel_df = sorf_excel_df[~sorf_excel_df['vtx_id'].str.contains('\|')]
    sorf_excel_df['index_copy'] = sorf_excel_df.index

    sorf_excel_df['show_details'] = False

    cols = list(sorf_excel_df.columns)
    cols.insert(0, cols.pop(cols.index('show_details')))

    sorf_excel_df = sorf_excel_df[cols]
    return sorf_excel_df


@st.cache_data()
def load_kibby_results(sorf_excel):
    kibby = pd.read_csv('../data/phase1to6_kibby.csv', index_col=0)
    kibby['conservation'] = kibby.conservation.apply(lambda x: list(map(float, x.strip().split(' '))))
    kibby = kibby.loc[sorf_excel['primary_id']]
    kibby.index = sorf_excel['vtx_id']
    return kibby


@st.cache_data()
def load_de_results(transcripts):
    cache_filez = os.listdir(CACHE_DIR)
    de_tables_dict = {}
    for f in cache_filez:
        if f.endswith('_de.parq') and not (f=='expression_de.parq'):
            df = pd.read_parquet(os.path.join(CACHE_DIR, f))
            df = df.loc[df.index.intersection(transcripts)]
            de_tables_dict[f.split('_')[0]] = df
    de_tables_dict = defaultdict(dict)
    for c, df in tqdm(de_tables_dict.items()):
        for row in df.itertuples():
            de_tables_dict[row[0]][c] = {'Cancer Average': row._7, 'GTEx Average': row._8, 
                                         'log2FC': row.log2FoldChange, 'padj': row.padj}
    for t, d in de_tables_dict.items():
        de_tables_dict[t] = pd.DataFrame(d).T
    return de_tables_dict


@st.cache_data()
def load_esmfold():
    """
    """
    esmfold = {}
    with jsonlines.open('../data/phase1to6_secreted_esmfold.json') as fopen:
        for l in fopen.iter():
            esmfold[l['sequence']] = l
    return esmfold


@st.cache_data()
def load_xena_tcga_gtex_target(vtx_combination_type='transcripts_exact'):
    # Expression is saved as TPM + 0.001 (NOT LOGGED)
    xena_expression = pd.read_parquet(os.path.join(CACHE_DIR, 'xena.parq'))
    xena_metadata = xena_expression[xena_expression.columns[:6]]
    xena_expression = xena_expression[xena_expression.columns[6:]]

    vtx_id_to_transcripts = load_jsonlines_table(os.path.join(CACHE_DIR, 'sorf_table.jsonlines'), index_col='vtx')
    # Map VTX to transcript ids - some 
    transcript_to_vtx_id = {}
    for ix, row in vtx_id_to_transcripts.iterrows():
        for val in row[vtx_combination_type]:
            transcript_to_vtx_id[val] = ix

    # Sum transcripts for VTX id
    xena_vtx_sums = xena_expression.T.copy()
    xena_vtx_sums = xena_vtx_sums.loc[xena_vtx_sums.index.intersection(transcript_to_vtx_id.keys())]
    xena_vtx_sums['vtx_id'] = xena_vtx_sums.apply(lambda x: transcript_to_vtx_id[x.name], axis=1)
    xena_vtx_sums = xena_vtx_sums.groupby('vtx_id').aggregate(np.sum).T
    de_tables_dict = load_de_results(list(transcript_to_vtx_id.keys()))
    return xena_metadata, xena_expression, vtx_id_to_transcripts, xena_vtx_sums, de_tables_dict


@st.cache_data()
def load_xena_heatmap(vtx_combination_type='transcripts_overlapping'):
    xena_metadata_df, xena_exp_df, _, xena_vtx_sum_df, _ = load_xena_tcga_gtex_target(vtx_combination_type)
    
    xena_vtx_sum_df = np.log2(xena_vtx_sum_df + 1)

    xena_vtx_exp_df = xena_metadata_df.merge(xena_vtx_sum_df, left_index=True, right_index=True)

    xena_vtx_sum_df = xena_vtx_sum_df.merge(xena_metadata_df, left_index=True, right_index=True)
    transcript_col_names = [i for i in xena_vtx_sum_df.columns if i not in xena_metadata_df.columns]
    xena_vtx_sum_df = xena_vtx_sum_df[transcript_col_names].groupby(xena_vtx_sum_df['primary disease or tissue']).aggregate(np.mean)

    threshold = .1
    mean_vals = xena_vtx_sum_df.max()
    cols_to_remove = mean_vals[mean_vals < threshold].index
    xena_vtx_sum_df = xena_vtx_sum_df.drop(cols_to_remove, axis=1)
    
    tau_df = xena_vtx_sum_df/xena_vtx_sum_df.max()
    tau = ((1-tau_df).sum())/(tau_df.shape[0]-1)
    tau.name = 'tau'
    xena_tau_df = xena_vtx_sum_df.T.merge(tau, left_index=True, right_index=True)

    # compute the clusters
    row_clusters = linkage(xena_vtx_sum_df.T.values, method='complete', metric='euclidean')
    col_clusters = linkage(xena_vtx_sum_df.values, method='complete', metric='euclidean')

    # compute the leaves order
    row_leaves = leaves_list(row_clusters)
    col_leaves = leaves_list(col_clusters)

    # reorder the DataFrame according to the clusters
    plot_df = xena_vtx_sum_df.T.iloc[row_leaves, col_leaves]

    col_map = {x: i for i,x in enumerate(plot_df.columns)}
    row_map = {x: i for i,x in enumerate(plot_df.index)}
    data = [(row_map[k[0]], col_map[k[1]], v) for k,v in plot_df.stack().items()]
    
    col_names = list(col_map.keys())
    row_names = list(row_map.keys())

    return data, col_names, row_names, xena_tau_df, xena_vtx_exp_df


@st.cache_data()
def load_mouse_blastp_results():
    hits_per_query = defaultdict(list)
    with open('../data/phase1to6_secreted_mouse_blastp.json', 'r') as fopen:
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
                align_str = '  \n'.join([alignment['qseq'], alignment['midline'], alignment['hseq']])
                alignment['hit_ids'] = ';'.join(ids)
                alignment['alignment'] = align_str
                hits_per_query[q].append(alignment)
    return hits_per_query


@st.cache_data
def convert_df(df):
    return df.to_csv().encode('utf-8')


def sorf_table(sorf_excel_df):
    st.title('sORF Library')
    st.write('Table contains secreted sORFs.')

    df = filter_dataframe(sorf_excel_df)

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

    update_df = st.data_editor(
        df,
        column_config={
            'vtx_id': st.column_config.TextColumn(disabled=True),
        },
        key='data_editor'
    )

    st.download_button(
        label="Download as CSV",
        data=convert_df(update_df),
        file_name='sorf_selection.csv',
        mime='text/csv',
    )

    st.session_state['data_editor_prev'] = st.session_state['data_editor'].copy()

    xena_metadata, xena_expression, vtx_id_to_transcripts, xena_vtx_sums, de_tables_dict = load_xena_tcga_gtex_target()
    esmfold = load_esmfold()
    blastp_mouse_hits = load_mouse_blastp_results()
    kibby = load_kibby_results(sorf_excel_df)
    protein_features_df = load_protein_feature_string_representations()

    if 'curr_vtx_id' in st.session_state.keys():

        vtx_id = st.session_state['curr_vtx_id']
        selected_row = df[df['vtx_id'] == st.session_state['curr_vtx_id']]

        st.header('sORF Details')
        st.dataframe(selected_row[['vtx_id', 'primary_id', 'orf_xref', 'protein_xrefs', 'gene_xref']])

        #link = ucsc_link(selected_row['chromosome'], selected_row['start'], selected_row['end'])
        #st.write(f"UCSC Browser [link]({link})")
        
        col1, col2 = st.columns(2)
        with col1:
            # Plot transcript expression levels
            fig = plotting.expression_heatmap_plot(vtx_id, vtx_id_to_transcripts, xena_expression, xena_metadata)
            if fig:
                st.pyplot(fig)
            else:
                st.write('No transcripts in TCGA/GTEx/TARGET found containing this sORF')
        
        # st.write(xena_overlap)
        # # Barplot Tumor vs NAT
        # fig, ax = plt.subplots(figsize=(9, 3))
        with col2:
            fig, result = plotting.expression_de_plot(vtx_id, vtx_id_to_transcripts, de_tables_dict)
            st.plotly_chart(fig)
            #with st.container():
            #    st.write(result)
            #    option = bar_plot_expression_groups(result, 'TCGA', ['GTEx', 'Cancer'])
            #    st_echarts(option, height="1000px")
            
        col1, col2 = st.columns(2)
        # Load esmfold data for selected sORF
        sorf_aa_seq = sorf_excel_df[sorf_excel_df['vtx_id']==vtx_id]['aa'].iloc[0]
        structure = esmfold[sorf_aa_seq]['pdb']
        plddt = esmfold[sorf_aa_seq]['plddt']
        # Plot esmfold structure
        st.title('sORF ESMfold')
        view = py3Dmol.view(js='https://3dmol.org/build/3Dmol.js',)
        view.addModel(structure, 'pdb')
        view.setStyle({'cartoon': {'colorscheme': {'prop':'b','gradient': 'roygb','min':50,'max':90}}})
        view.zoomTo()
        with col2:
            components.html(view._make_html(), height = 500,width=500)
        # showmol(view, height=500, width=1200)
        # Plot plDDT
        fig, axes = plt.subplots(2, 1, sharex=True)
        ax_plddt = plot_structure_plddt(plddt, axes[0])
        k = kibby.loc[vtx_id]['conservation']
        ax_kibby = axes[1].plot(k)
        f = protein_features_df[sorf_aa_seq]
        imdf, textdf = features_to_int_arrays(f)
        colorscale = [[0, '#454D59'], [0.25, '#FFFFFF'], [0.5, '#F1C40F'], [0.75, '#336641']]
        # fig = ff.create_annotated_heatmap(z = imdf, x=imdf.columns, y=imdf.index, text=textdf.values)
        plotly_fig = go.Figure(data=go.Heatmap(z=imdf.values, x=imdf.columns, y=imdf.index[:-1],
                                        colorscale='Picnic', text=textdf.values, zmin=0, zmax=3,
                                        showlegend=False, showscale=False))
        col1.pyplot(fig)
        col1.plotly_chart(plotly_fig)
        # Blastp Mouse
        primary_id = sorf_excel_df[sorf_excel_df['vtx_id'] == vtx_id].iloc[0]['primary_id']
        blastp_results_selected_sorf = blastp_mouse_hits[primary_id]
        if len(blastp_results_selected_sorf) == 0:
            long_text = "No alignments with mouse found." #st.write('No alignments with mouse found.')
        else:
            long_text = ""
            for h in blastp_results_selected_sorf:
                hit_text = f"Match IDs: {h['hit_ids']}  \nAlign Stats: Score - {h['score']}, Length - {h['align_len']}  \n"
                long_text+=hit_text
                long_text+= h['alignment'] + '  \n  \n'

        stx.scrollableTextbox(long_text,height = 300, fontFamily='Courier')


def sorf_transcriptome_atlas(sorf_excel_df):
    st.title("sORF Transcriptome Atlas")
    
    data, col_names, row_names, xena_tau_df, xena_vtx_exp_df = load_xena_heatmap()

    with st.container():
        col1, col2, col3 = st.columns(3)
        with col2:
            values = st.slider(
                'Select Tissue Specificity Tau',
                
                0.0, 1.0, (.8, 1.0),
                help='Higher values of [Tau](https://academic.oup.com/bib/article/18/2/205/2562739) indicate tissue specific expression of a given sORF')
        
        option, events, tissue_vtx_ids = plotting.expression_atlas_heatmap_plot(xena_tau_df, data, col_names, row_names, values)
        
        display_cols = ['vtx_id', 'primary_id', 'phase', 'orf_xref', 'protein_xrefs', 'gene_xref', 'transcript_xref', 'source', 'secreted_mean', 'translated_mean', 'isoform_of']
        value = st_echarts(option, height="1000px", events=events)
        if value:
            st.header('Selected sORF')
            st.dataframe(sorf_excel_df[sorf_excel_df['vtx_id'] == value][display_cols])
            
            fig = plotting.expression_vtx_boxplot(value, xena_vtx_exp_df)
            st.plotly_chart(fig, use_container_width=True)

        
        df = sorf_excel_df[sorf_excel_df['vtx_id'].isin(tissue_vtx_ids)][display_cols]
        st.header('Tissue specific sORFs')
        st.dataframe(df)
           
    
def selector(sorf_excel_df):
    st.title("Select sORFs Interactively")
    st.session_state['x_feature'] = st.selectbox("Select X Feature", options = sorf_excel_df.columns, index=11)
    st.session_state['y_feature'] = st.selectbox("Select Y Feature", options = sorf_excel_df.columns, index=12)
    
    fig = px.scatter(sorf_excel_df, x=st.session_state['x_feature'], y=st.session_state['y_feature'])
    selected_points = plotly_events(fig)
    st.write(selected_points)


def genome_browser():
    """
    """
    components.iframe("http://10.65.23.159:8080/velia_collections.html", height=1200, scrolling=True)


def main():
    #st.sidebar.title("Navigation")
    
    st.set_page_config(layout="wide")
    st.markdown(""" <style>iframe[title="streamlit_echarts.st_echarts"]{ height: 1000px !important } """, unsafe_allow_html=True)

    # Define a dictionary with the page names and corresponding functions
    pages = {
        "sORF Table": sorf_table,
        "sORF Transcriptome Atlas": sorf_transcriptome_atlas,
        "sORF Genome Browser": genome_browser,
        "sORF Selector": selector,
    }
    
    sorf_excel_df = load_sorf_excel()
  
    tab1, tab2, tab3, tab4 = st.tabs(list(pages.keys()))

    with tab1:
        sorf_table(sorf_excel_df)

    with tab2:
        sorf_transcriptome_atlas(sorf_excel_df)

    with tab3:
        genome_browser()
        
    with tab4:
        selector(sorf_excel_df)


# Run the app
if __name__ == "__main__":

    main()