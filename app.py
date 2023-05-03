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

from plotting import plot_structure_plddt

# Load data
@st.cache_data()
def load_sorf_excel():
    sorf_excel_table = pd.read_excel('./data/interim_phase1to6_secreted_hits_20230330.xlsx')
    return sorf_excel_table

@st.cache_data()
def load_esmfold():
    esmfold = {}
    with jsonlines.open('./data/phase1to6_secreted_esmfold.json') as fopen:
        for l in fopen.iter():
            esmfold[l['sequence']] = l
    return esmfold

@st.cache_data()
def load_xena_tcga_gtex_target():
    xena_metadata = pd.read_table('./data/TcgaTargetGTEX_phenotype.txt', encoding='latin-1', index_col=0)
    xena_expression = pd.read_feather('./data/xena_ucsc_phase1to6.feather')
    xena_expression.index = xena_expression.pop('index')
    xena_metadata = xena_metadata.loc[xena_expression.index]
    vtx_id_to_transcripts = json.load(open('./data/vtx_to_ensembl_ids.json', 'r'))
    return xena_metadata, xena_expression, vtx_id_to_transcripts

@st.cache_data()
def load_mouse_blastp_results():
    hits_per_query = defaultdict(list)
    with open('./data/phase1to6_secreted_mouse_blastp.json', 'r') as fopen:
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
    cancers = pd.read_csv('./data/cancer_types.txt', header=None, names = ['Disease', 'Code'])
    cancers.index = [i.lower() for i in cancers['Disease']]
    tcga_paired_normal = {}
    tcga_paired_normal_index = []
    pairs = pd.read_excel('./data/tissue_pairs.xlsx')
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

def sorf_table():
    st.title("sORF Explorer")
    st.write("Table contains secreted sORFs.")
    gd = GridOptionsBuilder.from_dataframe(sorf_excel_table)
    gd.configure_pagination(enabled=True)
    gd.configure_selection(use_checkbox=True)
    gridOptions = gd.build()
    ag_grid = AgGrid(sorf_excel_table, gridOptions = gridOptions,
                     # fit_columns_on_grid_load = True,
                     height=500,
                     width = "200%",
                     theme = "streamlit",
                     update_mode = GridUpdateMode.GRID_CHANGED,
                     reload_data = False,
                     editable = False)
def ucsc_link(chrom, start, end):
    base_link = "http://genome.ucsc.edu/cgi-bin/hgTracks?"
    genome = "db=hg38"
    position = f"{chrom}:{start}-{end}"
    ucsc_link = f"{base_link}{genome}&position={position}"
    return ucsc_link

# Function per page
def details():
    # Header text
    st.title('Details for a sORF')
    # Select a single sORF to plot and convert ID to vtx_id
    st.session_state['id_type_selected'] = st.radio("ID Type", options=('vtx_id', 'primary_id', 'genscript_id'), index=1)
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
    view = py3Dmol.view(js='https://3dmol.org/build/3Dmol.js',)
    view.addModel(structure, 'pdb')
    view.setStyle({'cartoon': {'colorscheme': {'prop':'b','gradient': 'roygb','min':50,'max':90}}})
    view.zoomTo()
    showmol(view, height=500, width=800)
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
                  
    
def selector():
    st.title("Select sORFs Interactively")
    st.session_state['x_feature'] = st.selectbox("Select X Feature", options = sorf_excel_table.columns, index=0)
    st.session_state['y_feature'] = st.selectbox("Select Y Feature", options = sorf_excel_table.columns, index=1)
    
    fig = px.scatter(sorf_excel_table, x=st.session_state['x_feature'], y=st.session_state['y_feature'])
    selected_points = plotly_events(fig)
    st.write(selected_points)

# Define the main function that will run the app
def main():
    st.sidebar.title("Navigation")
    selection = st.sidebar.radio("Go to", list(pages.keys()))

    # Call the appropriate function based on user selection
    page = pages[selection]
    page()

# Run the app
if __name__ == "__main__":
    st.set_page_config(layout="wide")
    # Define a dictionary with the page names and corresponding functions
    pages = {
        "sORF Table": sorf_table,
        "sORF Details": details,
        "sORF Selector": selector
    }
    sorf_excel_table = load_sorf_excel()
    xena_metadata, xena_expression, vtx_id_to_transcripts = load_xena_tcga_gtex_target()
    esmfold = load_esmfold()
    blastp_mouse_hits = load_mouse_blastp_results()
    ave_per_transcript_per_cancer = load_tcga_tumor_vs_nat(xena_metadata, xena_expression)
    main()