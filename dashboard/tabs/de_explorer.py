import numpy as np
import pandas as pd
import plotly.express as px
import streamlit as st
import streamlit.components.v1 as components
import streamlit_scrollable_textbox as stx


from streamlit_plotly_events import plotly_events

from streamlit_echarts import st_echarts
from scipy.cluster.hierarchy import linkage, leaves_list

from dashboard import plotting, description
from dashboard.etl import CACHE_DIR, DATA_DIR
from collections import defaultdict

import random
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
                    # labels={"logFC": "Log2(FoldChange)"},
                    hover_data=['transcript_id', 'vtx_id', 'contrast', 'velia_study'], render_mode='webgl')

    fig.update_layout(legend_font=dict(size=18))
    fig.update_yaxes(autorange="reversed")

    return fig

def de_page(sorf_df):
    db_address = DATA_DIR / "autoimmune_expression_atlas_v1.db"
    transcript_to_vtx_map = defaultdict(list)
    for ix, row in sorf_df.iterrows():
        for t in row['transcripts_exact']:
            transcript_to_vtx_map[t].append(ix)
    
    with sqlite3.connect(db_address) as sqliteConnection:
        available_studies = pd.read_sql("SELECT DISTINCT velia_study, contrast FROM transcript_de", sqliteConnection)
    
    filter_option = st.selectbox('Transcripts to Show:', ('All Transcripts',
                                                                    'sORF Transcripts Only'
                                                                    ), index = 1, key='de_selectbox')
    
    selection = dataframe_with_selections(available_studies)
     # Filter the dataframe using the temporary column, then drop the column
    # selection = edited_df[edited_df.Select]
    
    if selection.shape[0]!=0:
        st.write("Your selection:")
        st.write(selection)
        
        fdr = .5
        tpm_mean = 4
        log2fc = 1
        contrast = ', '.join([f"'{i}'" for i in selection['contrast'].values])#.iloc[0]
        velia_study = ', '.join([f"'{i}'" for i in selection['velia_study'].values])
        
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
        st.dataframe(gene_de_df)
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