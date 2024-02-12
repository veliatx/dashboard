import numpy as np
import pandas as pd
import streamlit as st

import plotly.graph_objects as go
import plotly.express as px

from dashboard import util
from dashboard.etl import DATA_DIR
from dashboard.data_load import load_autoimmune_metadata, load_xena_heatmap_data, query_s3_normed_counts
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

def plot_de_boxplots(
                    counts_df:pd.DataFrame,
                    stats_df:pd.DataFrame,
                    transcript:str,
                    vtx_id:str,
                ) -> go.Figure:
    """
    """
    # Color the boxes by study. 

    colormap = {c:px.colors.qualitative.G10[i%len(px.colors.qualitative.G10)] \
                for i,c in enumerate(counts_df['velia_study'].unique())}

    counts_df['color'] = counts_df['velia_study'].map(lambda x: colormap[x])
                    
    fig_box = px.box(
            counts_df,
            x='display',
            y=transcript,
            points='all',
            hover_data=['sample_condition_1','sample_condition_2','sample_type_1','sample_type_2'],
            hover_name='sample_id',
            color='color',
            boxmode='overlay',
        )
    fig_box.update_layout(
            xaxis_title='',
            yaxis_title='approximate_tpm',
            title=f'{transcript} -- {vtx_id}',
            title_x=0.5,
            title_xanchor='center',
            showlegend=False,
            )

    # Plot vertical lines separating experiments.

    exps_conds = counts_df['display'].unique()
    exps = [*map(lambda x: x.split('--')[0].strip(), exps_conds)]
    exps_idxs =  np.sort(np.unique(exps, return_index=True)[1])

    for i in exps_idxs[1:]:
        fig_box.add_vline(x=i-0.5, line_color='gray', line_dash='dash')

    # Plot horizontal lines denoting significance.
    
    min_counts = counts_df.loc[:,transcript].min()
    max_counts = counts_df.loc[:,transcript].max()

    padj_cutoff = 0.01
    lfc_cutoff = 1.0

    padj = stats_df.iloc[stats_df.index.get_level_values('stat') == 'padj'][[transcript]]
    lfc = stats_df.iloc[stats_df.index.get_level_values('stat') == 'log2FoldChange'][[transcript]]           
                        
    plot_sig = np.where((padj[transcript] < padj_cutoff).values & (abs(lfc[transcript]) > lfc_cutoff).values)[0]
                    
    for i, y in enumerate(np.linspace(
                    max_counts + max_counts*0.1, 
                    max_counts + ((max_counts-min_counts)*0.4), 
                    plot_sig.shape[0],
                    )):
        left_loc = np.where(exps_conds == padj.iloc[[plot_sig[i]]].index.get_level_values('display_left').values)[0][0]
        right_loc = np.where(exps_conds == padj.iloc[[plot_sig[i]]].index.get_level_values('display_right').values)[0][0]
        fig_box.add_shape(
                        type="line", 
                        x0=min(left_loc,right_loc), 
                        y0=y, 
                        x1=max(left_loc,right_loc), 
                        y1=y, 
                        label=dict(text='*', textposition='middle center', font=dict(size=12, color='black')),
                        ) 
    fig_box.update_yaxes(range=[-5.,max_counts + ((max_counts-min_counts)*(0.1*(plot_sig.shape[0]+1)))])  

    return fig_box

@st.cache_data(ttl='45m')
def build_normed_counts_dfs(
                selected_contrasts_df:pd.DataFrame,
                contrast:pd.DataFrame,
                meta:pd.DataFrame,
                selected_transcripts:list,
                ):
    """
    Merges counts dataframes queried from s3 into one dataframe for plotting. df_stats is df with multi-index set on stat, 
    velia_study, contrast, display_left, and display_right.  
    
    Args:
        selected_contrasts pd.DataFrame
        contrast pd.DataFrame
        meta pd.DataFrame
        selected_transcripts list
    Returns:
        df_counts, df_stats Tuple[pd.DataFrame,pd.DataFrame] 
    """
    dfs = [query_s3_normed_counts(r.velia_study, r.contrast) for i,r in selected_contrasts_df.iterrows()]
    
    df_counts = pd.concat(
                [df.loc[
                     df.index.isin(selected_transcripts+['velia_study','contrast']),
                     ~df.columns.isin(('log2FoldChange','padj'))
                    ] for df in dfs], 
                ignore_index=False,
                axis=1,
                ).fillna(0.0).T

    df_stats = pd.concat(
                [df.loc[
                    df.index.isin(selected_transcripts+['velia_study','contrast']),
                    ['log2FoldChange','padj']
                    ].T for df in dfs],
                ignore_index=False,
                axis=0,
                )
                    
    df_stats.loc['log2FoldChange'] = df_stats.loc['log2FoldChange'].fillna(0.0)
    df_stats.loc['padj'] = df_stats.loc['padj'].fillna(1.0)

    df_stats[['left','right']] = df_stats['contrast'].str.upper().str.split('_VS_', expand=True)
    df_stats['display_left'] = df_stats.apply(lambda x: f'{x.velia_study} -- {x.left}', axis=1) 
    df_stats['display_right'] = df_stats.apply(lambda x: f'{x.velia_study} -- {x.right}', axis=1)
    df_stats.drop(['left','right'],axis=1, inplace=True)
    df_stats.index.name = 'stat'
    df_stats.set_index(['velia_study','contrast','display_left','display_right'], append=True, inplace=True)

    df_counts = df_counts.merge(
            contrast[
                contrast['contrast'].isin(selected_contrasts_df['contrast'])
                ][['sample_id','display','contrast_side']],
            left_index=True,
            right_on='sample_id',
            ).merge(
            meta.fillna(''),
            on='sample_id',
            ).sort_values(['velia_study','contrast','contrast_side'], ascending=False)

    df_counts.drop_duplicates('sample_id', inplace=True)
                    
    return df_counts, df_stats


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
    
    contrast_samples, sample_meta_df = load_autoimmune_metadata()

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
        
        fdr = .5
        tpm_mean = 4
        log2fc = 1
        contrast = ', '.join([f"'{i}'" for i in selection_df['contrast'].values])#.iloc[0]
        velia_study = ', '.join([f"'{i}'" for i in selection_df['velia_study'].values])
        
        with sqlite3.connect(db_address) as sqliteConnection:
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
        
        if gene_de_df.shape[0]>0:
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
                hoverdata = volcano_fig.data[x['curveNumber']]['customdata'][x['pointIndex']]
                if isinstance(hoverdata[2], list):
                    selected_vtx_ids.append((hoverdata[2][0], hoverdata[0]))
                else:
                    selected_vtx_ids.append((hoverdata[2], hoverdata[0]))


            # st.write(selected_vtx_ids)
            # st.write(len(volcano_fig.data), len(volcano_fig.data[0]['customdata']), len(volcano_fig.data[1]['customdata']))
            selected_vtx_ids = list(set(selected_vtx_ids))

            gene_de_sorf_df = gene_de_df.explode('vtx_id')
            gene_de_sorf_df = gene_de_sorf_df[
                                gene_de_sorf_df['vtx_id'].isin([x[0] for x in selected_vtx_ids]) & \
                                gene_de_sorf_df['transcript_id'].isin([x[1] for x in selected_vtx_ids])
                                ]
            gene_de_sorf_df = gene_de_sorf_df[['vtx_id','velia_study','contrast','transcript_id','log2FoldChange','padj','Significant']] \
                                .merge(sorf_df.reset_index(drop=True), on='vtx_id') \
                                .drop('show_details', axis=1) \
                                .set_index('vtx_id')
            
            try:
                # dashboard.tabs.sorf_explorer_table.sorf_details(sorf_df.loc[[x[0] for x in selected_vtx_ids]].copy())
                dashboard.tabs.sorf_explorer_table.sorf_details(gene_de_sorf_df)
            except Exception as e:
                print(e)
                # st.dataframe(sorf_df.loc[[x[0] for x in selected_vtx_ids]].copy().drop('show_details', axis=1))
                st.dataframe(gene_de_sorf_df)
            
            st.title('DE boxplots')
   
            if len(selected_vtx_ids) > 0:

                boxplot_studies = st.multiselect(   
                                        'Select contrasts for boxplot',
                                        list(available_studies_df.index),
                                        default=list(selection_df.index),
                                    )

                selected_df_boxplots = available_studies_df.loc[boxplot_studies]

                if selected_df_boxplots.shape[0] > 0:

                    normed_counts_df, stats_df = build_normed_counts_dfs(
                                                        selected_df_boxplots,
                                                        contrast_samples,
                                                        sample_meta_df,
                                                        [x[1] for x in selected_vtx_ids],
                                                    )
                    
                    selected_vtx_id, selected_transcript_id = st.selectbox(
                                                        'Selected transcripts:', 
                                                        selected_vtx_ids, 
                                                        format_func=lambda x: f'{x[0]} -- {x[1]}',
                                                    )
                    boxplot_fig = plot_de_boxplots(
                                                normed_counts_df,
                                                stats_df,
                                                selected_transcript_id,
                                                selected_vtx_id,
                                            )
                st.plotly_chart(boxplot_fig)
            
        else:
            st.write('No transcripts with these filtering were found.')
    # dashboard.tabs.sorf_explorer_table.sorf_details(sorf_df.loc[selected_vtx_ids])
    
    # st.dataframe(availabe_studies)
