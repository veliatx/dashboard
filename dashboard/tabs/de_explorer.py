import pandas as pd
import plotly.express as px
import streamlit as st

from dashboard import util
from dashboard.etl import DATA_DIR
from dashboard.data_load import load_autoimmune_contrasts, load_xena_heatmap_data
from dashboard.tabs.expression_heatmap import filter_weak_expression, assign_top_tissues

from collections import defaultdict
from streamlit_plotly_events import plotly_events
import dashboard.tabs.sorf_explorer_table
import dashboard.tabs.riboseq_atlas
import streamlit.components.v1 as components
import smart_open



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
    plot_gene_df.fillna('', inplace=True)
    fig = px.scatter(plot_gene_df, x='log2FoldChange', y='log10_padj', opacity=0.5,
                     color="Significant", color_discrete_map={True: "blue", False: "red"},
                     size_max=20, template='plotly_white',
                     labels={"log2FoldChange": "Log2(FoldChange)"},
                     hover_data=['transcript_id', 'hgnc_name', 'vtx_id', 'contrast', 'velia_study', 'tissues', 'tau'], render_mode='webgl')
    fig.update_traces(marker=dict(size=12))
    fig.update_layout(legend_font=dict(size=18))
    fig.update_yaxes(autorange="reversed")

    return fig


def download_qc_report_from_s3(velia_study):
    file_stem = f"s3://velia-piperuns-dev/expression_atlas/v1/{velia_study}/"
    with smart_open.smart_open(file_stem, 'r') as fopen:
        data = fopen.readlines()
    return data


def de_page(sorf_df):

    filter_option = st.selectbox('Pre-filtered sORFs:', ('Ribo-Seq sORFs',
                                                        'Secreted',
                                                        'Secreted & Novel',
                                                        'Secreted & Conserved',
                                                        'Secreted & Conserved & Novel',
                                                        'Translated',
                                                        'Translated & Conserved',
                                                        'All sORFs'), index = 0, key='sorf_detail_filter2')
    
    sorf_df = util.filter_dataframe_preset(sorf_df, filter_option)
    transcript_to_hgnc = pd.read_csv(DATA_DIR / 'veliadb_v1.transcript2hgnc.csv.gz', index_col=0)
    db_address = DATA_DIR.joinpath('autoimmune_expression_atlas_v1.db')
    
    xena_tau_df, _, _ = load_xena_heatmap_data()
    xena_tau_df = filter_weak_expression(xena_tau_df, sorf_df)
    #xena_tau_df = xena_tau_df[xena_tau_df['tau'].between(.85, 1)].copy()
    xena_tau_df, _ = assign_top_tissues(xena_tau_df)
    
    transcript_to_vtx_map = defaultdict(list)
    for ix, row in sorf_df.iterrows():
        for t in row['transcripts_exact']:
            transcript_to_vtx_map[t].append(ix)

    st.caption(f"{sorf_df.shape[0]} ORFs across {len(transcript_to_vtx_map)} transcripts")

    with sqlite3.connect(db_address) as sqliteConnection:
        available_studies_df = pd.read_sql("SELECT DISTINCT velia_study, contrast FROM transcript_de", sqliteConnection)
    
    study_ex = 'GSE194263'
    contrast_ex = 'SJOGRENS_SYNDROME_vs_CONTROL'
    query_str = f'not (velia_study == @study_ex and contrast == @contrast_ex)'
    available_studies_df = available_studies_df.query(query_str)

    filter_option = st.selectbox('Transcripts to Show:', ('All Transcripts', 'sORF Transcripts Only'),
                                  index = 1, key='de_selectbox')

    available_studies_df.index = available_studies_df.apply(lambda x: f"{x['velia_study']} -- {x['contrast']}", axis=1)
    
    contrast_samples = load_autoimmune_contrasts()
    # contrast_samples = load_autoimmune_contrasts().groupby(['study -- contrast','contrast_side']).agg(','.join)
    # contrast_samples['samples_left'] = contrast_samples.loc[contrast_samples.index.get_level_values('contrast_side') == 'left']['sample_id']
    # contrast_samples['samples_right'] = contrast_samples.loc[contrast_samples.index.get_level_values('contrast_side') == 'right']['sample_id']
    # contrast_samples.reset_index('contrast_side', inplace=True, drop=True)
    # contrast_samples = contrast_samples[['samples_left','samples_right']].groupby('study -- contrast').agg(sum)

    selection = st.multiselect(
        'Choose one or more studies', 
         ['All Studies'] + list(available_studies_df.index))
    
    if 'All Studies' in selection:
        selection_df = available_studies_df
    else:
        selection_df = available_studies_df.loc[selection]
    # for ix, row in selection_df.iterrows():
    #     st.download_button(row['velia_study'], download_qc_report_from_s3(velia_study), file_name=f"{row['velia_study']}.html")
    
    if selection_df.shape[0] > 0:
        st.caption(f"Your selection: {selection_df.shape[0]} studies")

        with sqlite3.connect(db_address) as sqliteConnection:
            sample_meta_df = pd.read_sql(
                                f"""SELECT
                                        sample_id, 
                                        sample_condition_1,
                                        sample_condition_2,
                                        sample_condition_3,
                                        sample_type_1,
                                        sample_type_2
                                    FROM sample_metadata where velia_study in 
                                        ({(','.join(selection_df['velia_study'].map(lambda x: "'"+x+"'").unique()))});
                                """,
                                sqliteConnection,
                            )
            
        selection_df = selection_df.merge(
                                        contrast_samples[['sample_id','contrast_side']], 
                                        left_index=True, 
                                        right_index=True, 
                                        how='left',
                                    )
        selection_df = selection_df.merge(sample_meta_df, on='sample_id', how='left').fillna('')

        # Naming sample_n here just so it gets captured by the below startswith.
        # Casting to str to accomodate the groupby below.

        selection_df.insert(
                        1, 
                        'sample_n',
                        selection_df.groupby(
                            ['velia_study', 'contrast', 'contrast_side']
                            )['sample_id'].transform(len).astype(str),
                    )
        selection_df.drop('sample_id', inplace=True, axis=1)

        selection_df = selection_df.groupby(
                                        ['velia_study', 'contrast', 'contrast_side'],
                                        ).agg(lambda x: ','.join(x.unique()))
        
        for c in selection_df.columns[selection_df.columns.str.startswith('sample')]:
             selection_df[f'left_{c}'] = selection_df.loc[
                                                selection_df.index.get_level_values('contrast_side') == 'left'][c]
             selection_df[f'right_{c}'] = selection_df.loc[
                                                selection_df.index.get_level_values('contrast_side') == 'right'][c]
        
        selection_df.reset_index(inplace=True, drop=False)
        selection_df.drop(
                    selection_df.columns[selection_df.columns.str.startswith('sample')].tolist()+['contrast_side'], 
                    inplace=True, 
                    axis=1,
                )

        selection_df.index = selection_df.apply(lambda x: f"{x['velia_study']} -- {x['contrast']}", axis=1)
        selection_df = selection_df.groupby(
                                    [selection_df.index, 'velia_study', 'contrast']
                                    ).agg(sum).reset_index(['velia_study','contrast'], drop=False)

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
            'gene_id', 'case_mean', 'control_mean', 'tau', 'tissues'
        ]
        

        gene_de_df['vtx_single'] = gene_de_df.apply(lambda x: x.vtx_id[0], axis=1)
        gene_de_df = gene_de_df.merge(xena_tau_df[['tau', 'tissues']],
                                      left_on='vtx_single', 
                                      right_index=True, how='left')

        gene_de_df = util.filter_dataframe_dynamic(gene_de_df, f'gene_de_filter')

        vtx_cnt = gene_de_df['vtx_id'].astype(str).nunique()
        tx_cnt = gene_de_df['transcript_id'].nunique()      

        st.caption(f"{vtx_cnt} unique uPs on {tx_cnt} DE transcripts")
        st.dataframe(gene_de_df[display_cols])

        gene_de_df['Significant'] = gene_de_df['padj'] < 0.01
        gene_de_df['hgnc_name'] = [transcript_to_hgnc.loc[i] if i in transcript_to_hgnc.index else 'na' for i in gene_de_df['transcript_id']]
        
        volcano_fig = plot_gene_volcano(gene_de_df)
        selected_points = plotly_events(volcano_fig, click_event=True, hover_event=False, select_event=True)
        # st.plotly_chart(volcano_fig)
        
        # st.write(selected_points)
        selected_vtx_ids = []
        for x in selected_points:
            hoverdata = volcano_fig.data[x['curveNumber']]['customdata'][x['pointIndex']][2]
            selected_vtx_ids += hoverdata
            # st.write(hoverdata)
        # st.write(selected_vtx_ids)
        # st.write(len(volcano_fig.data), len(volcano_fig.data[0]['customdata']), len(volcano_fig.data[1]['customdata']))
        selected_vtx_ids = list(set(selected_vtx_ids))
        try:
            dashboard.tabs.sorf_explorer_table.sorf_details(sorf_df.loc[selected_vtx_ids].copy())
        except Exception as e:
            print(e)
            st.dataframe(sorf_df.loc[selected_vtx_ids].copy().drop('show_details', axis=1))
    # dashboard.tabs.sorf_explorer_table.sorf_details(sorf_df.loc[selected_vtx_ids])
    
    # st.dataframe(availabe_studies)
