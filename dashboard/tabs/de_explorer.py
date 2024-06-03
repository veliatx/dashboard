import numpy as np
import pandas as pd
import streamlit as st
from typing import List, Tuple
import inspect

import plotly.graph_objects as go
import plotly.express as px

from dashboard import util
from dashboard.etl import DATA_DIR, DB_CONNECTION_STRING, REDSHIFT_CONNECTION_STRING
from dashboard.data_load import load_xena_heatmap_data
from dashboard.tabs.expression_heatmap import filter_weak_expression, assign_top_tissues

from expression_atlas_db import base, queries

from collections import defaultdict
from streamlit_plotly_events import plotly_events
import dashboard.tabs.sorf_explorer_table
import dashboard.tabs.riboseq_atlas
import streamlit.components.v1 as components
import smart_open

# # Auto-wrap all the query functions in st.cache_data. 

# for name, f in inspect.getmembers(queries, inspect.isfunction):
#     setattr(
#         queries, 
#         name, 
#         st.cache_data(f, ttl='1h', hash_funcs={
#             'sqlalchemy.orm.session._Session': lambda _: None,
#             'sqlalchemy.orm.attributes.InstrumentedAttribute': lambda _: None,
#             'builtins.function': lambda _: None,
#             },
#         ),
#     )

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


def plot_gene_volcano(plot_gene_df: pd.DataFrame, w_vtx:bool=True) -> Tuple[go.Figure, List[str]]:
    """
    """
    if w_vtx:
        hover_labels=['transcript_id', 'hgnc_name', 'vtx_id', 'contrast', 'velia_study', 'tissues', 'tau'] 
    else:
        hover_labels=['transcript_id', 'hgnc_name', 'contrast', 'velia_study']

    plot_gene_df.fillna('', inplace=True)
    fig = px.scatter(
        plot_gene_df, 
        x='log2FoldChange', 
        y='log10_padj', 
        opacity=0.5,
        color="Significant", 
        color_discrete_map={True: "blue", False: "red"},
        size_max=20, 
        template='plotly_white',
        labels={"log2FoldChange": "Log2(FoldChange)"},
        hover_data=hover_labels, 
        render_mode='webgl',
    )
    fig.update_traces(marker=dict(size=12))
    fig.update_layout(legend_font=dict(size=18))
    fig.update_yaxes(autorange="reversed")

    return fig, hover_labels

def plot_de_boxplots(
                    counts_df:pd.DataFrame,
                    stats_df:pd.DataFrame,
                    transcript:str,
                    label_id:str,
                ) -> go.Figure:
    """
    """
    # Color the boxes by study. 

    colormap = {c:px.colors.qualitative.G10[i%len(px.colors.qualitative.G10)] \
                for i,c in enumerate(counts_df['velia_study'].unique())}

    counts_df['color'] = counts_df['velia_study'].map(lambda x: colormap[x])

    # Set velia_study_contrast in counts_df and stats_df, used for grouping boxplot and pulling stats.
                    
    counts_df['velia_study_contrast'] = counts_df.apply(
        lambda x: f'{x.velia_study} -- {x.contrast} -- {x.contrast_side}', 
        axis=1,
    )
    stats_df.reset_index(inplace=True, drop=False)
    stats_df['velia_study_contrast_left'] = stats_df.apply(
        lambda x: f'{x.velia_study} -- {x.contrast} -- left',
        axis=1,
    )
    stats_df['velia_study_contrast_right'] = stats_df.apply(
        lambda x: f'{x.velia_study} -- {x.contrast} -- right',
        axis=1,
    )
    stats_df.set_index([
            'stat',
            'velia_study',
            'contrast',
            'display_left',
            'display_right',
        ], 
        inplace=True,
    )
    counts_df.sort_values('velia_study_contrast', inplace=True, ascending=False)
    stats_df.sort_values(['velia_study','contrast'], inplace=True, ascending=False)

    fig_box = px.box(
        counts_df,
        x='velia_study_contrast',
        y=transcript,
        points='all',
        hover_data={
                c:False if c not in \
                ['sample_condition_1','sample_condition_2','sample_type_1','sample_type_2'] \
                else True for c in counts_df.columns
                },
        hover_name='sample_id',
        color='color',
        boxmode='overlay',
        category_orders={'velia_study_contrast': counts_df['velia_study_contrast'].unique()},
    )
                    
    # Plot vertical lines separating experiments.

    exps_conds = counts_df['velia_study_contrast'].unique()
    exps = [*map(lambda x: '--'.join(x.split('--')[:-1]).strip(), exps_conds)]
    exps_idxs = np.sort(np.unique(exps, return_index=True)[1])
                    
    for i in exps_idxs[1:]:
        fig_box.add_vline(x=i-0.5, line_color='gray', line_dash='dash')

    # Plot horizontal lines denoting significance.
    
    max_counts = counts_df.loc[:,transcript].max()

    padj_cutoff = 0.01
    lfc_cutoff = 1.0

    padj = stats_df.iloc[stats_df.index.get_level_values('stat') == 'padj']
    lfc = stats_df.iloc[stats_df.index.get_level_values('stat') == 'log2foldchange']                            
    plot_sig = np.where((padj[transcript] < padj_cutoff).values & (abs(lfc[transcript]) > lfc_cutoff).values)[0]
                    
    for l in plot_sig:
        left_loc = np.where(exps_conds == padj.iloc[l]['velia_study_contrast_left'])[0][0]
        right_loc = np.where(exps_conds == padj.iloc[l]['velia_study_contrast_right'])[0][0]
        fig_box.add_shape(
            type="line", 
            x0=left_loc, 
            y0=max_counts + max_counts*0.1, 
            x1=right_loc, 
            y1=max_counts + max_counts*0.1, 
            label=dict(text='*', textposition='middle center', font=dict(size=12, color='black')),
        ) 
        
    fig_box.update_layout(
        xaxis_title='',
        yaxis_title='approximate_tpm',
        title=f'{transcript} -- {label_id}',
        title_x=0.5,
        title_xanchor='center',
        showlegend=False,
        height=600,
        yaxis_range=[-0.05*(max_counts),max_counts + max_counts*0.25],
    )

    return fig_box

def download_qc_report_from_s3(velia_study):
    file_stem = f"s3://velia-piperuns-dev/expression_atlas/v1/{velia_study}/"
    with smart_open.smart_open(file_stem, 'r') as fopen:
        data = fopen.readlines()
    return data


def de_page(sorf_df):

    Session = base.configure(DB_CONNECTION_STRING)
    SessionRedshift = base.configure(REDSHIFT_CONNECTION_STRING)

    filter_option = st.selectbox('Pre-filtered sORFs:', (
                                                        'Secreted',
                                                        'Secreted & Novel',
                                                        'Secreted & Conserved',
                                                        'Secreted & Conserved & Novel',
                                                        ), index = 0, key='sorf_detail_filter2')
    
    if not st.session_state.get('sorf_df') or filter_option != st.session_state['sorf_df']['filter_option']:
        
        # Clear session_state.
        for k in st.session_state.keys():
            st.session_state.pop(k)

        with st.spinner('Pulling dataframes...'):
            sorf_df = util.filter_dataframe_preset(sorf_df, filter_option)
            transcript_to_hgnc = pd.read_csv(DATA_DIR / 'veliadb_v1.transcript2hgnc.csv.gz', index_col=0)
            
            xena_tau_df, _, _ = load_xena_heatmap_data()
            xena_tau_df = filter_weak_expression(xena_tau_df, sorf_df)
            #xena_tau_df = xena_tau_df[xena_tau_df['tau'].between(.85, 1)].copy()
            xena_tau_df, _ = assign_top_tissues(xena_tau_df)
            
            transcript_to_vtx_map = defaultdict(list)
            for ix, row in sorf_df.iterrows():
                for t in row['transcripts_exact']:
                    transcript_to_vtx_map[t].append(ix)

            st.session_state['sorf_df'] = {
                'sorf_df': sorf_df.copy(),
                'xena_tau_df': xena_tau_df.copy(),
                'transcript_to_vtx_map': transcript_to_vtx_map.copy(),
                'filter_option':filter_option,
                'transcript_to_hgnc': transcript_to_hgnc.copy(),
            }
    else:
        transcript_to_vtx_map = st.session_state['sorf_df']['transcript_to_vtx_map'].copy()
        xena_tau_df = st.session_state['sorf_df']['xena_tau_df'].copy()
        sorf_df = st.session_state['sorf_df']['sorf_df'].copy()
        transcript_to_hgnc = st.session_state['sorf_df']['transcript_to_hgnc'].copy()

    st.caption(f"{sorf_df.shape[0]} ORFs across {len(transcript_to_vtx_map)} transcripts")

    with Session.begin() as session, \
        st.spinner('Loading studies...'):
        available_studies_df = queries.fetch_contrasts(session)
        available_studies_df.rename({'velia_id':'velia_study','contrast_name':'contrast'}, axis=1, inplace=True)
        available_studies_df = available_studies_df[['velia_study','contrast']]
            
    study_ex = 'GSE194263'
    contrast_ex = 'SJOGRENS_SYNDROME_vs_CONTROL'
    query_str = f'not (velia_study == @study_ex and contrast == @contrast_ex)'
    available_studies_df = available_studies_df.query(query_str)

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
    # for ix, row in selection_df.iterrows():
    #     st.download_button(row['velia_study'], download_qc_report_from_s3(velia_study), file_name=f"{row['velia_study']}.html")
    
    if selection_df.shape[0] > 0:
        st.caption(f"Your selection: {selection_df.shape[0]} studies")
          
        if not st.session_state.get('selected_contrasts') or \
            st.session_state['selected_contrasts'][0] != set(selection_df.index.tolist()):

            with Session.begin() as session, st.spinner('Loading contrast metadata...'):
                sample_contrasts_df = queries.build_contrast_metatable(
                    session, 
                    studies=selection_df['velia_study'].tolist(), 
                    contrasts=selection_df['contrast'].tolist(),
                )
                sample_contrasts_df.rename({'velia_id':'velia_study','contrast_name':'contrast'}, axis=1, inplace=True)
                st.session_state['selected_contrasts'] = (set(selection_df.index.tolist()), sample_contrasts_df.copy())
        else:
            sample_contrasts_df = st.session_state['selected_contrasts'][1].copy()
        st.dataframe(sample_contrasts_df)
        
        fdr = .5
        tpm_mean = 4
        log2fc = 1
       
        if not st.session_state.get('differential_expression') or \
            st.session_state['differential_expression'][0] != set(selection_df.index.tolist()):

            with Session.begin() as session, SessionRedshift.begin() as session_redshift, \
                st.spinner('Loading differentially-expressed transcripts...'):
                gene_de_df = queries.query_differentialexpression(
                    session,
                    session_redshift,
                    studies=selection_df['velia_study'].tolist(),
                    contrasts=selection_df['contrast'].tolist(),
                    sequenceregions_type='transcript',
                    log10_padj_threshold=np.log10(fdr),
                    log2_fc_threshold=log2fc,
                    mean_threshold=tpm_mean,
                )
                gene_de_df['vtx_id'] = [transcript_to_vtx_map.get(t,'None') for t in gene_de_df['transcript_id']]
                gene_de_df.rename({
                    'velia_id':'velia_study',
                    'contrast_name':'contrast',
                    'basemean':'baseMean',
                    'log2foldchange':'log2FoldChange',
                    'lfcse':'lfcSE',
                    }, axis=1, inplace=True,
                )
                gene_de_df.reset_index(inplace=True, drop=True)
                st.session_state['differential_expression'] = (set(selection_df.index.tolist()), gene_de_df.copy())
        else:
            gene_de_df = st.session_state['differential_expression'][1].copy()
            
        if filter_option == 'sORF Transcripts Only':
            gene_de_df = gene_de_df[gene_de_df['vtx_id']!='None']
        
        if gene_de_df.shape[0]>0:

            if filter_option == 'sORF Transcripts Only':

                display_cols = [
                    'velia_study', 'contrast', 'vtx_id', 'transcript_id', 'baseMean',
                    'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj', 'log10_padj',
                    'hgnc_name', 'case_mean', 'control_mean', 'tau', 'tissues', 
                    'SNPS', 'MAPPED_TRAIT', 'P-VALUE'
                ]

                gene_de_df['vtx_single'] = gene_de_df.apply(lambda x: x.vtx_id[0], axis=1)
                gene_de_df = gene_de_df.merge(xena_tau_df[['tau', 'tissues']],
                                            left_on='vtx_single', 
                                            right_index=True, how='left')
            
                gene_de_df = gene_de_df.merge(
                    sorf_df[['SNPS', 'MAPPED_TRAIT', 'P-VALUE']],
                    left_on='vtx_single', 
                    right_index=True, 
                    how='left',
                )
            else:
                display_cols = [
                    'velia_study', 'contrast', 'transcript_id', 'baseMean',
                    'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj', 'log10_padj',
                    'hgnc_name', 'case_mean', 'control_mean', 
                ]
            
            gene_de_df = util.filter_dataframe_dynamic(gene_de_df, f'gene_de_filter')

            vtx_cnt = gene_de_df['vtx_id'].astype(str).nunique()
            tx_cnt = gene_de_df['transcript_id'].nunique()      

            gene_de_df['Significant'] = gene_de_df['padj'] < 0.01
            gene_de_df['hgnc_name'] = [transcript_to_hgnc.loc[i]['hgnc_name'] if i in transcript_to_hgnc.index else 'na' for i in gene_de_df['transcript_id']]
            
            if filter_option == 'sORF Transcripts Only':
                st.caption(f"{vtx_cnt} unique uPs on {tx_cnt} DE transcripts")
            else:
                st.caption(f"{tx_cnt} DE transcripts")

            st.dataframe(gene_de_df[display_cols])

            volcano_fig, hover_labels = plot_gene_volcano(gene_de_df, w_vtx=filter_option == 'sORF Transcripts Only')
            selected_points = plotly_events(volcano_fig, click_event=False, hover_event=False, select_event=True)

            selected_pts_ids = []

            # If set for sORF transcripts, filter here to pull vtx_id out of hoverdata or hgnc_name.

            if filter_option == 'sORF Transcripts Only':
                pt0_idx = hover_labels.index('vtx_id') if 'vtx_id' in hover_labels else 2
            else:
                pt0_idx = hover_labels.index('hgnc_name')
            pt1_idx = hover_labels.index('transcript_id')

            for x in selected_points:
                hoverdata = volcano_fig.data[x['curveNumber']]['customdata'][x['pointIndex']]

                if isinstance(hoverdata[pt0_idx], list):
                    selected_pts_ids.append((hoverdata[pt0_idx][0], hoverdata[pt1_idx]))
                else:
                    selected_pts_ids.append((hoverdata[pt0_idx], hoverdata[pt1_idx]))
            
            selected_pts_ids = list(set(selected_pts_ids))

            if len(selected_pts_ids) > 100:
                st.write(f"That's a lot of microproteins! You chose {len(selected_pts_ids)}. Try choosing a few less.")
                selected_pts_ids = []

            if filter_option == 'sORF Transcripts Only':

                gene_de_sorf_df = gene_de_df.explode('vtx_id')
                gene_de_sorf_df = gene_de_sorf_df[
                    gene_de_sorf_df['vtx_id'].isin([x[0] for x in selected_pts_ids]) & \
                    gene_de_sorf_df['transcript_id'].isin([x[1] for x in selected_pts_ids])
                ]
                gene_de_sorf_df = gene_de_sorf_df[['vtx_id','velia_study','contrast','transcript_id','log2FoldChange','padj','Significant']] \
                    .merge(sorf_df.reset_index(drop=True), on='vtx_id') \
                    .drop('show_details', axis=1) \
                    .set_index('vtx_id')
                
                try:
                    dashboard.tabs.sorf_explorer_table.sorf_details(gene_de_sorf_df.copy())
                except Exception as e:
                    print(e)
                    st.dataframe(gene_de_sorf_df)
            else:
                st.title('Selected Transcripts')
                gene_de_df = gene_de_df.loc[:,gene_de_df.columns != 'vtx_id']
                st.dataframe(gene_de_df[gene_de_df['transcript_id'].isin([x[1] for x in selected_pts_ids])])
            
            st.title('DE boxplots')
   
            if len(selected_pts_ids) > 0:

                boxplot_studies = st.multiselect(   
                    'Select contrasts for boxplot',
                    list(available_studies_df.index),
                    default=list(selection_df.index),
                )

                selected_df_boxplots = available_studies_df.loc[boxplot_studies]

                if selected_df_boxplots.shape[0] > 0:

                    if not st.session_state.get('sample_measurement') or \
                        st.session_state['sample_measurement'][0] != set(selected_pts_ids) or \
                        st.session_state['sample_measurement'][1] != set(boxplot_studies):
                        with Session.begin() as session, SessionRedshift.begin() as session_redshift, \
                            st.spinner('Loading expression measurements...'):
                            normed_counts_df = queries.query_samplemeasurement(
                                session,
                                session_redshift,
                                studies=selected_df_boxplots['velia_study'].tolist(),
                                contrasts=selected_df_boxplots['contrast'].tolist(),
                                sequenceregions=[x[1] for x in selected_pts_ids],
                                exact_id_match=False,
                            )
                            normed_counts_df.rename(
                                {'velia_id':'velia_study', 'contrast_name':'contrast', 'srx_id':'sample_id'}, 
                                axis=1, 
                                inplace=True,
                            )
                            normed_counts_df = pd.pivot(
                                normed_counts_df, 
                                index=['velia_study','contrast','contrast_side','sample_id']+\
                                    [c for c in normed_counts_df.columns if c.startswith('sample_condition')], 
                                values='normed_counts_transform', 
                                columns='transcript_id',
                            ).fillna(0.0).reset_index()
                            normed_counts_df['display'] = normed_counts_df.apply(
                            lambda x: f"{x.velia_study} -- " +
                                f"{x.contrast.upper().split('VS')[0 if x.contrast_side == 'left' else 1].strip('_')}", 
                            axis=1,
                            )
                            normed_counts_df = normed_counts_df.sort_values(
                                ['velia_study','contrast','contrast_side'], 
                                ascending=False,
                            )
                            
                            stats_df = queries.query_differentialexpression(
                                session,
                                session_redshift,
                                studies=selected_df_boxplots['velia_study'].tolist(),
                                contrasts=selected_df_boxplots['contrast'].tolist(),
                                sequenceregions=[x[1] for x in selected_pts_ids],
                                log10_padj_threshold=None,
                                log2_fc_threshold=None,
                                mean_threshold=None,
                                exact_id_match=False,
                            )
                            stats_df.rename({'velia_id':'velia_study','contrast_name':'contrast'}, axis=1, inplace=True)
                            stats_df[['left','right']] = stats_df['contrast'].str.upper().str.split('_VS_', expand=True)
                            stats_df['display_left'] = stats_df.apply(lambda x: f'{x.velia_study} -- {x.left}', axis=1)
                            stats_df['display_right'] = stats_df.apply(lambda x: f'{x.velia_study} -- {x.right}', axis=1)
                            stats_df.drop(['left','right'],axis=1, inplace=True)
                            lfc_df = pd.pivot(
                                stats_df, 
                                values='log2foldchange', 
                                columns='transcript_id', 
                                index=['velia_study','contrast','display_left','display_right'],
                            ).fillna(0.0)
                            lfc_df['stat'] = 'log2foldchange'
                            padj_df = pd.pivot(
                                stats_df, 
                                values='padj', 
                                columns='transcript_id', 
                                index=['velia_study','contrast','display_left','display_right'],
                            ).fillna(1.0)
                            padj_df['stat'] = 'padj'
                            stats_df = pd.concat([lfc_df,padj_df])
                            stats_df.set_index(['stat'], append=True, inplace=True)
                            stats_df = stats_df.sort_values(
                                ['velia_study','contrast','display_left','display_right','stat'], 
                                ascending=False,
                            )

                            st.session_state['sample_measurement'] = (
                                set(selected_pts_ids), 
                                set(boxplot_studies), 
                                normed_counts_df.copy(), 
                                stats_df.copy(),
                            )
                    else:
                        _, _, normed_counts_df, stats_df = st.session_state['sample_measurement']
                 
                    selected_label_id, selected_transcript_id = st.selectbox(
                        'Selected transcripts:', 
                        selected_pts_ids, 
                        format_func=lambda x: f'{x[0]} -- {x[1]}',
                    )
                    boxplot_fig = plot_de_boxplots(
                        normed_counts_df,
                        stats_df,
                        selected_transcript_id,
                        selected_label_id,
                    )
                st.plotly_chart(boxplot_fig, use_container_width=True)
            
        else:
            st.write('No transcripts with these filtering were found.')
    # dashboard.tabs.sorf_explorer_table.sorf_details(sorf_df.loc[selected_pts_ids])
    
    # st.dataframe(availabe_studies)
