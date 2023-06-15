import streamlit as st
import numpy as np
from streamlit_plotly_events import plotly_events
import plotly.express as px
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
from stmol import showmol
import py3Dmol
from scipy.stats import ranksums
from statsmodels.stats.multitest import multipletests
from streamlit_echarts import st_echarts
from scipy.cluster.hierarchy import linkage, leaves_list


from pandas.api.types import (
    is_categorical_dtype,
    is_datetime64_any_dtype,
    is_numeric_dtype,
    is_object_dtype,
)

from plotting import plot_structure_plddt

# Load data
@st.cache_data()
def load_sorf_excel():
    """
    """
    sorf_excel_table = pd.read_excel('../data/interim_phase1to6_secreted_hits_20230330.xlsx')
    return sorf_excel_table


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
def load_xena_tcga_gtex_target():
    xena_metadata = pd.read_table('../data/TcgaTargetGTEX_phenotype.txt', encoding='latin-1', index_col=0)
    xena_expression = pd.read_feather('../data/xena_ucsc_phase1to6.feather')
    xena_expression.index = xena_expression.pop('index')
    xena_metadata = xena_metadata.loc[xena_expression.index]
    vtx_id_to_transcripts = json.load(open('../data/vtx_to_ensembl_ids.json', 'r'))
    return xena_metadata, xena_expression, vtx_id_to_transcripts

@st.cache_data()
def load_xena_heatmap():
    
    xena_metadata_df, xena_exp_df, vtx_id_to_transcripts = load_xena_tcga_gtex_target()
    
    transcript_to_vtx_id = {}

    for k, vals in vtx_id_to_transcripts.items():
        for val in vals:
            transcript_to_vtx_id[val] = k
    
    xena_exp_df = xena_exp_df.T
    xena_exp_df['vtx_id'] = xena_exp_df.apply(lambda x: transcript_to_vtx_id[x.name], axis=1)
    xena_exp_df = xena_exp_df.groupby('vtx_id').aggregate(np.mean).T
    xena_exp_df = xena_exp_df.merge(xena_metadata_df, left_index=True, right_index=True)
    transcript_col_names = [i for i in xena_exp_df.columns if i not in xena_metadata_df.columns]
    xena_exp_df = xena_exp_df[transcript_col_names].groupby(xena_exp_df['_primary_site']).aggregate(np.mean)

    # compute the clusters
    row_clusters = linkage(xena_exp_df.T.values, method='complete', metric='euclidean')
    col_clusters = linkage(xena_exp_df.values, method='complete', metric='euclidean')

    # compute the leaves order
    row_leaves = leaves_list(row_clusters)
    col_leaves = leaves_list(col_clusters)

    # reorder the DataFrame according to the clusters
    plot_df = xena_exp_df.T.iloc[row_leaves, col_leaves]

    col_map = {x: i for i,x in enumerate(plot_df.columns)}
    row_map = {x: i for i,x in enumerate(plot_df.index)}
    data = [(row_map[k[0]], col_map[k[1]], v) for k,v in plot_df.stack().items()]
    
    col_names = list(col_map.keys())
    row_names = list(row_map.keys())

    return data, col_names, row_names

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


@st.cache_data()
def load_tcga_tumor_vs_nat(xena_metadata, xena_expression):
    cancers = pd.read_csv('../data/cancer_types.txt', header=None, names = ['Disease', 'Code'])
    cancers.index = [i.lower() for i in cancers['Disease']]
    tcga_paired_normal = {}
    tcga_paired_normal_index = []
    pairs = pd.read_excel('../data/tissue_pairs.xlsx')
    for ix, row in pairs.iterrows():
        if isinstance(row['NAT (Solid Tissue Normal)'], str):
            tcga_groups = xena_metadata[(xena_metadata['_primary_site'] == row['Tissue (Disease)']) & (xena_metadata['_study'] == 'TCGA')].copy()
            for cancer in tcga_groups['primary disease or tissue'].unique():
                primary_tumor_samples = tcga_groups[(tcga_groups['primary disease or tissue']==cancer) & 
                                                   (tcga_groups['_sample_type']=='Primary Tumor')].index
                normal_adjacent_samples = tcga_groups[(tcga_groups['primary disease or tissue']==cancer) & 
                                                   (tcga_groups['_sample_type']=='Solid Tissue Normal')].index
                if len(primary_tumor_samples)>=8 and (len(normal_adjacent_samples)>=8):
                    tcga_paired_normal[cancer] = {'Tumor': primary_tumor_samples, 'NAT': normal_adjacent_samples}
                    tcga_paired_normal_index += list(primary_tumor_samples)+list(normal_adjacent_samples)
        else:
            continue
    tcga_nat_table = xena_expression.loc[tcga_paired_normal_index].copy()
    tcga_nat_table.insert(0, 'Cancer', '')
    tcga_nat_table.insert(1, 'Condition', '')
    for c, s in tcga_paired_normal.items():
        tcga_nat_table.loc[s['Tumor'], 'Cancer'] = c
        tcga_nat_table.loc[s['NAT'], 'Cancer'] = c
        tcga_nat_table.loc[s['Tumor'], 'Condition'] = 'Cancer'
        tcga_nat_table.loc[s['NAT'], 'Condition'] = 'Normal Adjacent'
    ave_per_transcript_per_cancer = tcga_nat_table.groupby(['Cancer', 'Condition']).median()
    ave_per_transcript_per_cancer = ave_per_transcript_per_cancer.reset_index()
    stats_per_tumor = {}
    for c, s in tcga_paired_normal.items():
        tumor = tcga_nat_table.loc[s['Tumor'], tcga_nat_table.columns[2:]]
        normal = tcga_nat_table.loc[s['NAT'], tcga_nat_table.columns[2:]]
        logfcs = tumor.mean(axis=0) - normal.mean(axis=0)
        test_statistic, p_value = ranksums(tumor, normal, alternative='two-sided')

        stats = pd.DataFrame({'ranksum': test_statistic, 'p_value': p_value, 'logFC': logfcs.values, 
                 'FDR': multipletests(p_value, method='fdr_bh')[1]}, index=tumor.columns)
        stats_per_tumor[c] = stats.copy()
    logfcs = tumor.mean(axis=0) - normal.mean(axis=0)
    return ave_per_transcript_per_cancer



@st.cache
def convert_df(df):
    # IMPORTANT: Cache the conversion to prevent computation on every rerun
    return df.to_csv().encode('utf-8')



def filter_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    """
    Adds a UI on top of a dataframe to let viewers filter columns

    Args:
        df (pd.DataFrame): Original dataframe

    Returns:
        pd.DataFrame: Filtered dataframe
    """

    df = df.copy()

    # Try to convert datetimes into a standard format (datetime, no timezone)
    for col in df.columns:
        if is_object_dtype(df[col]):
            try:
                df[col] = pd.to_datetime(df[col])
            except Exception:
                pass

        if is_datetime64_any_dtype(df[col]):
            df[col] = df[col].dt.tz_localize(None)

    modification_container = st.container()

    with modification_container:
        to_filter_columns = st.multiselect("Filter dataframe on", df.columns)
        for column in to_filter_columns:
            left, right = st.columns((1, 20))
            # Treat columns with < 10 unique values as categorical
            if is_categorical_dtype(df[column]) or df[column].nunique() < 10:
                user_cat_input = right.multiselect(
                    f"Values for {column}",
                    df[column].unique(),
                    default=list(df[column].unique()),
                )
                df = df[df[column].isin(user_cat_input)]
            elif is_numeric_dtype(df[column]):
                _min = float(df[column].min())
                _max = float(df[column].max())
                step = (_max - _min) / 100
                user_num_input = right.slider(
                    f"Values for {column}",
                    min_value=_min,
                    max_value=_max,
                    value=(_min, _max),
                    step=step,
                )
                df = df[df[column].between(*user_num_input)]
            elif is_datetime64_any_dtype(df[column]):
                user_date_input = right.date_input(
                    f"Values for {column}",
                    value=(
                        df[column].min(),
                        df[column].max(),
                    ),
                )
                if len(user_date_input) == 2:
                    user_date_input = tuple(map(pd.to_datetime, user_date_input))
                    start_date, end_date = user_date_input
                    df = df.loc[df[column].between(start_date, end_date)]
            else:
                user_text_input = right.text_input(
                    f"Substring or regex in {column}",
                )
                if user_text_input:
                    df = df[df[column].astype(str).str.contains(user_text_input)]

    return df


def sorf_table(sorf_excel_table):
    st.title("sORF Library")
    st.write("Table contains secreted sORFs.")
    df = filter_dataframe(sorf_excel_table)
    st.dataframe(df)
    st.download_button(
        label="Download data as CSV",
        data=df,
        file_name='large_df.csv',
        mime='text/csv',
    )


def sorf_heatmap():
    st.title("sORF Transcriptome Atlas")
    st.write("Table contains secreted sORFs.")
    
    data, col_names, row_names = load_xena_heatmap()

    option = {
        "tooltip": {},
        "xAxis": {
            "type": "category", 
            "data": row_names, 
            "axisLabel": {
                "fontSize": 10,
                "rotate": -90,
                "width": 100,
            }
            },
        "yAxis": {
            "type": "category", 
            "data": col_names, 
            },
        "visualMap": {
            "min": 0,
            "max": 10,
            "calculable": True,
            "realtime": False,
            "orient": "horizontal",
            "left": "center",
            "top": "0%",
        },
        "series": [
            {
                "name": "Log2(TPM+1)",
                "type": "heatmap",
                "data": data,
                #"label": {"show": True},
                "emphasis": {
                    "itemStyle": {
                       "borderColor": '#333',
                        "borderWidth": 1
                    }
                },
                "progressive": 1000,
                "animation": False,
            }
        ],
    }
    with st.container():
        st_echarts(option, height="1000px")



def ucsc_link(chrom, start, end):
    base_link = "http://genome.ucsc.edu/cgi-bin/hgTracks?"
    genome = "db=hg38"
    position = f"{chrom}:{start}-{end}"
    ucsc_link = f"{base_link}{genome}&position={position}"
    return ucsc_link


# Function per page
def details(sorf_excel_table, xena_expression, xena_metadata, vtx_id_to_transcripts, ave_per_transcript_per_cancer, esmfold, blastp_mouse_hits):
    # Header text
    st.title('sORF Details')
    # Select a single sORF to plot and convert ID to vtx_id
    st.session_state['id_type_selected'] = st.radio("ID Type", options=('vtx_id', 'primary_id', 'genscript_id'), index=0, horizontal=True)
    st.session_state['detail_id_seleceted'] = st.selectbox("Select ORF", options = sorf_excel_table[st.session_state['id_type_selected']])
    if st.session_state['id_type_selected'] == 'vtx_id':
        vtx_id = st.session_state['detail_id_seleceted']
    else:
        vtx_id = sorf_excel_table[sorf_excel_table[st.session_state['id_type_selected']] == st.session_state['detail_id_seleceted']].iloc[0]['vtx_id']
    selected_row = sorf_excel_table[sorf_excel_table[st.session_state['id_type_selected']] == st.session_state['detail_id_seleceted']].iloc[0]
    link = ucsc_link(selected_row['chromosome'], selected_row['start'], selected_row['end'])
    st.write(f"UCSC Browser [link]({link})")

    # Plot transcript expression levels
    selected_transcripts = vtx_id_to_transcripts[vtx_id]
    xena_overlap = xena_expression.columns.intersection(selected_transcripts)
    if len(xena_overlap)>1:
        selected_expression = xena_expression[xena_overlap]
        groups = list(map(lambda x: '-'.join(map(str, x)), xena_metadata[['_primary_site', '_study']].values))
        cmap_fig = sns.clustermap(selected_expression.groupby(groups).median(),
                                  cmap='coolwarm', cbar_kws={'label': 'Log(TPM+0.001)'}, center=1, vmin=-3, vmax=6)
        st.pyplot(cmap_fig)
    elif len(xena_overlap) == 1:
        st.write('Only 1 transcript')
    else:
        st.write('No transcripts in TCGA/GTEx/TARGET found containing this sORF')
    st.write(xena_overlap)
    # Barplot Tumor vs NAT
    fig, ax = plt.subplots(figsize=(9, 3))
    exp_bplot = sns.barplot(data = ave_per_transcript_per_cancer, x='Cancer', ax = ax,
            y=ave_per_transcript_per_cancer[xena_overlap].apply(lambda x: np.log2(np.sum(np.exp2(x)-0.001)+0.001), axis=1),
            hue='Condition', alpha=0.75, palette=["r", "k"])
    exp_bplot.set_xticks(exp_bplot.get_xticks(), exp_bplot.get_xticklabels(), rotation=90)
    exp_bplot.set_ylabel('Log2 (TPM+0.001)')
    st.pyplot(fig)
    # Load esmfold data for selected sORF
    sorf_aa_seq = sorf_excel_table[sorf_excel_table['vtx_id']==vtx_id]['aa'].iloc[0]
    structure = esmfold[sorf_aa_seq]['pdb']
    plddt = esmfold[sorf_aa_seq]['plddt']
    # Plot esmfold structure
    st.title('sORF ESMfold')
    view = py3Dmol.view(js='https://3dmol.org/build/3Dmol.js',)
    view.addModel(structure, 'pdb')
    view.setStyle({'cartoon': {'colorscheme': {'prop':'b','gradient': 'roygb','min':50,'max':90}}})
    view.zoomTo()
    showmol(view, height=500, width=1200)
    # Plot plDDT
    fig, ax = plot_structure_plddt(plddt)
    st.pyplot(fig)
    # Blastp Mouse
    if st.session_state['id_type_selected'] == 'primary_id':
        primary_id = st.session_state['detail_id_seleceted']
    else:
        primary_id = sorf_excel_table[sorf_excel_table[st.session_state['id_type_selected']] == st.session_state['detail_id_seleceted']].iloc[0]['primary_id']
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
                  
    
def selector(sorf_excel_table):
    st.title("Select sORFs Interactively")
    st.session_state['x_feature'] = st.selectbox("Select X Feature", options = sorf_excel_table.columns, index=0)
    st.session_state['y_feature'] = st.selectbox("Select Y Feature", options = sorf_excel_table.columns, index=1)
    
    fig = px.scatter(sorf_excel_table, x=st.session_state['x_feature'], y=st.session_state['y_feature'])
    selected_points = plotly_events(fig)
    st.write(selected_points)


# Define the main function that will run the app
def main():
    #st.sidebar.title("Navigation")
    
    st.set_page_config(layout="wide")
    st.markdown(""" <style>iframe[title="streamlit_echarts.st_echarts"]{ height: 1000px !important } """, unsafe_allow_html=True)

    # Define a dictionary with the page names and corresponding functions
    pages = {
        "sORF Table": sorf_table,
        "sORF Transcriptome Atlas": sorf_heatmap,
        "sORF Details": details,
        "sORF Selector": selector
    }
    
    sorf_excel_table = load_sorf_excel()
    xena_metadata, xena_expression, vtx_id_to_transcripts = load_xena_tcga_gtex_target()
    esmfold = load_esmfold()
    blastp_mouse_hits = load_mouse_blastp_results()
    ave_per_transcript_per_cancer = load_tcga_tumor_vs_nat(xena_metadata, xena_expression)
    
    tab1, tab2, tab3, tab4 = st.tabs(list(pages.keys()))

    with tab1:
        sorf_table(sorf_excel_table)

    with tab2:
        sorf_heatmap()
    
    with tab3:
        details(sorf_excel_table, xena_expression, xena_metadata, vtx_id_to_transcripts, ave_per_transcript_per_cancer, esmfold, blastp_mouse_hits)
        
    with tab4:
        selector(sorf_excel_table)

# Run the app
if __name__ == "__main__":

    main()