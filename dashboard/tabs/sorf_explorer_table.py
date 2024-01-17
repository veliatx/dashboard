import py3Dmol

import numpy as np
import pandas as pd
import plotly.express as px
import streamlit as st
import streamlit.components.v1 as components
import streamlit_scrollable_textbox as stx

from streamlit_echarts import st_echarts

from dashboard import plotting, description, data_load, util
from dashboard.etl import DATA_DIR, CACHE_DIR

import sqlite3

def sorf_details(sorf_df):
    st.title('sORF Table')
    
    filter_option = st.selectbox('Pre-filtered sORFs:', ('Ribo-Seq sORFs',
                                                        'Secreted',
                                                        'Secreted & Conserved',
                                                        'Secreted & Conserved & Novel',
                                                        'Translated',
                                                        'Translated & Conserved',
                                                        'All sORFs'), index = 0, key='sorf_detail_filter')
    
    df = util.filter_dataframe_preset(sorf_df, filter_option)

    df = util.filter_dataframe_dynamic(df, f'explorer_filter')
    df.sort_index(inplace=True)

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
        },
        key='data_editor',
        hide_index=True
    )

    st.download_button(
        label="Download as CSV",
        data=util.convert_df(update_df),
        file_name='sorf_selection.csv',
        mime='text/csv',
    )

    st.session_state['data_editor_prev'] = st.session_state['data_editor'].copy()

    if 'curr_vtx_id' in st.session_state.keys():

        # Load data
        xena_metadata, xena_transcript_ids = data_load.load_xena_metadata()
        autoimmune_metadata = data_load.load_autoimmune_atlas()
        esmfold = data_load.load_esmfold(CACHE_DIR.joinpath('protein_data', 'esmfold.jsonlines'))
        blastp_mouse_hits, blastp_data_for_sorf_table = data_load.load_mouse_blastp_results()
        protein_features_df = data_load.load_protein_feature_string_representations()
        phylocsf_dataframe = data_load.load_phylocsf_data()
        xena_overlap = []

        vtx_id = st.session_state['curr_vtx_id']
        selected_row = df[df['vtx_id'] == st.session_state['curr_vtx_id']]

        st.divider()
        st.header('sORF Details')
        st.dataframe(selected_row[['vtx_id', 'screening_phase_id', 'orf_xrefs', 'protein_xrefs', 'gene_xrefs', 'aa', 'transcripts_exact']])

        selected_transcripts = sorf_df.loc[vtx_id, 'transcripts_exact']
        
        xena_overlap = xena_transcript_ids.intersection(set([i.split('.')[0] for i in selected_transcripts]))

        value = None
        with st.expander("Transcription Data", expanded=True):
            col1, col2 = st.columns(2)

            with col1:
                title = f'TCGA/GTEx Transcript Specific Expression - {vtx_id}'
                selected_expression_tcga = pd.read_parquet(CACHE_DIR.joinpath('xena_app.parq'), columns=xena_overlap)
                selected_expression_tcga_ave = selected_expression_tcga.groupby(xena_metadata['dashboard_group']).median()
                echart_option_tcga, events_tcga = plotting.expression_heatmap_plot(title, 
                                                                                   selected_expression_tcga_ave,
                                                                                   list(xena_overlap))
                
                if echart_option_tcga:
                    value_tcga = st_echarts(echart_option_tcga,
                                       height="900px", 
                                       events=events_tcga, 
                                       renderer='svg',
                                       key = 'tcga_echart_heatmap_explorer'
                                       )
                else:
                    st.write('No transcripts in TCGA/GTEx/TARGET found containing this sORF')
                    
                st.title('Autoimmune Expression Atlas')
                title = f'Autoimmune Atlas Transcript Specific Expression - {vtx_id}'
                
                formatted_ids = ', '.join(f"'{id_}'" for id_ in selected_transcripts)
                sql_query = f"SELECT * FROM transcript_tpm WHERE transcript_tpm.transcript_id IN ({formatted_ids});"
                selected_expression_ai = pd.read_sql(sql_query, sqlite3.connect(DATA_DIR.joinpath('autoimmune_expression_atlas_v1.db'))).fillna(0.01)
                selected_expression_ai_ave = selected_expression_ai.pivot_table(index='group',
                                                                        columns='transcript_id',
                                                                        values='tpm', 
                                                                        aggfunc=np.nanmedian).fillna(0.01)#.apply(lambda x: np.log2(x+1))
                sample_sizes = selected_expression_ai['group'].value_counts()
                selected_expression_ai_ave.index = [f"{x} n={sample_sizes[x]}" for x in selected_expression_ai_ave.index]
                echart_option_ai, events_ai = plotting.expression_heatmap_plot(title,
                                                                               selected_expression_ai_ave,
                                                                               median_groups=False)
                if echart_option_ai:
                    value_ai = st_echarts(echart_option_ai,
                                       height="900px",
                                       events=events_ai,
                                       renderer='svg',
                                       key = 'ai_echart_heatmap_explorer'
                                       )
                else:
                    st.write('No transcripts in Velia Autoimmune Atlas found containing this sORF')

            with col2:

                if len(xena_overlap) > 0:
                    if value_tcga:
                        selected_transcript_tcga = value_tcga
                    else:
                        selected_transcript_tcga = list(xena_overlap)[0]

                    chart_title = f'Differential Expression - {selected_transcript_tcga}'
                    de_exact_echarts_options_b = plotting.bar_plot_expression_groups_tcga(selected_transcript_tcga.split('.')[0], 
                                                                                            'TCGA', ['GTEx Mean', 'Cancer Mean'],
                                                                                            chart_title)                                     
                    st_echarts(options=de_exact_echarts_options_b, key='b', height='900px', width = '600px', renderer='svg')
                
                st.title('')
                if len(selected_transcripts) > 0:
                    if value_ai:
                        selected_transcript_ai = value_ai
                    else:
                        selected_transcript_ai = selected_transcripts[0]
                    db_address = DATA_DIR.joinpath('autoimmune_expression_atlas_v1.db')
                    option_ai_de = plotting.bar_plot_expression_group_autoimmune(util.query_de_transcripts(selected_transcript_ai, db_address).fillna(0.01),
                                                                                 f'Autoimmune DE - {selected_transcript_ai}',
                                                                                 db_address)
                    st_echarts(options=option_ai_de, key='c', height='900px', width = '650px', renderer='svg')
                    
            if (len(xena_overlap) > 0) and value_tcga:
                st.write(value_tcga)
                xena_vtx_exp_df = xena_metadata.merge(selected_expression_tcga, left_index=True, right_index=True)
                fig_tcga = plotting.expression_vtx_boxplot(selected_transcript_tcga.split('.')[0], xena_vtx_exp_df)
                st.plotly_chart(fig_tcga, use_container_width=True)
                
            if (len(selected_transcripts) > 0) and value_ai:    
                fig_ai = px.box(data_frame = selected_expression_ai[selected_expression_ai['transcript_id']==selected_transcript_ai].sort_values('group').rename({'tpm': selected_transcript_ai}, axis=1),
                    x='group', points = 'all',
                    y=selected_transcript_ai, height=500,
                    width=800 
                )
                st.plotly_chart(fig_ai, use_container_width=True)


        with st.expander("Protein Structure and Function", expanded=True):
            col3, col4 = st.columns(2)
            try:
                # Load esmfold data for selected sORF
                sorf_aa_seq = sorf_df[sorf_df['vtx_id']==vtx_id]['aa'].iloc[0]
                if sorf_aa_seq[-1] == '*':
                    sorf_aa_seq = sorf_aa_seq[:-1]
                plddt = esmfold[sorf_aa_seq]['plddt']
                structure = esmfold[sorf_aa_seq]['pdb']

                with col3:
                    # Plot plDDT and Phylocsf
                    achart = plotting.plot_sequence_line_plots_altair(vtx_id, sorf_aa_seq, phylocsf_dataframe, esmfold)
                    col3.altair_chart(achart, use_container_width=False)
                            
                with col4:
                    modified_structure_colors = plotting.color_protein_terminal_ends(sorf_aa_seq, structure)
                    view = py3Dmol.view(js='https://3dmol.org/build/3Dmol.js',)
                    view.addModel(modified_structure_colors, 'pdb')
                    view.setStyle({'cartoon': {'colorscheme': {'prop':'b','gradient': 'roygb','min':50,'max':90}}})
                    view.zoomTo()
                    st.header('sORF ESMfold', help="Red - Low Confidence  \nBlue - High Confidence  \nConfidence is based on plDDT score from ESMFold  \nN-term is blue and C-term is red")
                    components.html(view._make_html(), height=500, width=600)
            except:
                print('No structure and function for this AA')

            f = protein_features_df[vtx_id]
            imdf = plotting.format_protein_feature_strings_for_altair_heatmap(f)
            altair_signal_features_fig = plotting.altair_protein_features_plot(imdf)
            st.header('AA Feature Predictions', help=description.amino_acid_features_hover_text)
            st.altair_chart(altair_signal_features_fig, use_container_width=True)

        
        with st.expander("BLASTp results", expanded=True):
            # Blastp Mouse
            blastp_results_selected_sorf = blastp_mouse_hits[vtx_id]
            if len(blastp_results_selected_sorf) == 0:
                long_text = "No alignments with mouse found."
            else:
                long_text = ""
                for h in blastp_results_selected_sorf:
                    hit_text = f"Match IDs: {h['hit_ids']}  \nAlign Stats: Score - {h['score']}, Length - {h['align_len']}  \n"
                    long_text+=hit_text
                    long_text+= h['alignment'] + '  \n  \n'
    
            stx.scrollableTextbox(long_text, height = 300, fontFamily='Courier')
