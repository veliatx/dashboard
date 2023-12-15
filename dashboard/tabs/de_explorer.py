import pandas as pd
import plotly.express as px
import streamlit as st

from dashboard import util
from dashboard.etl import DATA_DIR
from collections import defaultdict
from streamlit_plotly_events import plotly_events

import sqlite3


def dataframe_with_selections(df):
    df_with_selections = df.copy()
    df_with_selections.insert(0, "Select", False)
    edited_df = st.data_editor(
        df_with_selections,
        hide_index=True,
        column_config={"Select": st.column_config.CheckboxColumn(required=True)},
        disabled=df.columns,
    )
    # Filter the dataframe using the temporary column
    selected_rows = edited_df[edited_df["Select"]]
    return selected_rows.drop("Select", axis=1)


def plot_gene_volcano(plot_gene_df):
    """
    """
    fig = px.scatter(plot_gene_df, x='log2FoldChange', y='log10_padj', opacity=0.5,
                     color="Significant", color_discrete_map={True: "blue", False: "red"},
                     size_max=10, template='plotly_white',
                     labels={"log2FoldChange": "Log2(FoldChange)"},
                     hover_data=['transcript_id', 'vtx_id', 'contrast', 'velia_study'], render_mode='webgl')

    fig.update_layout(legend_font=dict(size=18))
    fig.update_yaxes(autorange="reversed")

    return fig


def de_page(sorf_df):

    filter_option = st.selectbox('Pre-filtered sORFs:', ('Ribo-Seq sORFs',
                                                         'Secreted',
                                                         'Secreted & Conserved',
                                                         'Secreted & Conserved & Novel',
                                                         'Translated',
                                                         'Translated & Conserved',
                                                         'All sORFs'), index = 0, key='de_explorer_filter')
    
    
    sorf_df = util.filter_dataframe_preset(sorf_df, filter_option)

    db_address = DATA_DIR.joinpath('autoimmune_expression_atlas_v1.db')
    transcript_to_vtx_map = defaultdict(list)
    for ix, row in sorf_df.iterrows():
        for t in row['transcripts_exact']:
            transcript_to_vtx_map[t].append(ix)

    st.caption(f"{sorf_df.shape[0]} ORFs across {len(transcript_to_vtx_map)} transcripts")

    with sqlite3.connect(db_address) as sqliteConnection:
        available_studies_df = pd.read_sql("SELECT DISTINCT velia_study, contrast FROM transcript_de", sqliteConnection)
    
    filter_option = st.selectbox('Transcripts to Show:', ('All Transcripts', 'sORF Transcripts Only'),
                                  index = 1, key='de_selectbox')
    

    available_studies_df.index = available_studies_df.apply(lambda x: f"{x['velia_study']} -- {x['contrast']}", axis=1)
    selection = st.multiselect(
        'Choose one or more studies', 
         ['All Studies'] + list(available_studies_df.index))
    
    if 'All Studies' in selection:
        selection_df = available_studies_df
    else:
        selection_df = available_studies_df.loc[selection]
    
    #selection = dataframe_with_selections(available_studies_df)
     # Filter the dataframe using the temporary column, then drop the column
    # selection = edited_df[edited_df.Select]
    
    if selection_df.shape[0] > 0:
        st.caption(f"Your selection: {selection_df.shape[0]} studies")
        st.dataframe(selection_df)
        #st.write(selection)
        
        fdr = .5
        tpm_mean = 4
        log2fc = 1
        contrast = ', '.join([f"'{i}'" for i in selection_df['contrast'].values])#.iloc[0]
        velia_study = ', '.join([f"'{i}'" for i in selection_df['velia_study'].values])
        
        with sqlite3.connect(db_address) as sqliteConnection:
            sqliteConnection = sqlite3.connect(db_address)
            query = f"""SELECT *
            FROM transcript_de
            WHERE transcript_de.padj <= {fdr}
            AND transcript_de.contrast IN ({contrast})
            AND transcript_de.velia_study IN ({velia_study})
            AND ABS(transcript_de.log2FoldChange) > {log2fc}
            AND (transcript_de.case_mean >= {tpm_mean} OR transcript_de.control_mean >= {tpm_mean})
            """
            gene_de_df = pd.read_sql(query, sqliteConnection)
            gene_de_df['vtx_id'] = [transcript_to_vtx_map[t] if t in transcript_to_vtx_map else 'None' for t in gene_de_df['transcript_id']]
        
        if filter_option == 'sORF Transcripts Only':
            gene_de_df = gene_de_df[gene_de_df['vtx_id']!='None']
        
        display_cols = [
            'velia_study', 'contrast', 'vtx_id', 'transcript_id', 'baseMean',
            'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj', 'log10_padj',
            'gene_id', 'case_mean', 'control_mean', 
        ]
        
        vtx_cnt = gene_de_df['vtx_id'].astype(str).nunique()
        tx_cnt = gene_de_df['transcript_id'].nunique()

        st.caption(f"{vtx_cnt} unique uPs on {tx_cnt} DE transcripts")
        st.dataframe(gene_de_df[display_cols])

        gene_de_df['Significant'] = gene_de_df['padj'] < 0.01
        volcano_fig = plot_gene_volcano(gene_de_df)
        selected_points = plotly_events(volcano_fig, click_event=True, hover_event=False, select_event=True)
        # st.plotly_chart(volcano_fig)
        
        # st.write(selected_points)
        selected_vtx_ids = []
        for x in selected_points:
            hoverdata = volcano_fig.data[x['curveNumber']]['customdata'][x['pointIndex']][1]
            selected_vtx_ids += hoverdata
            # st.write(hoverdata)
        # st.write(selected_vtx_ids)
        # st.write(len(volcano_fig.data), len(volcano_fig.data[0]['customdata']), len(volcano_fig.data[1]['customdata']))
    
        st.dataframe(sorf_df.loc[selected_vtx_ids].copy().drop('show_details', axis=1))
    # dashboard.tabs.sorf_explorer_table.sorf_details(sorf_df.loc[selected_vtx_ids])
    
    # st.dataframe(availabe_studies)