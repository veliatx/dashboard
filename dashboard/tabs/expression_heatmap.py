import streamlit as st
from streamlit_plotly_events import plotly_events
from streamlit_echarts import st_echarts
from scipy.cluster.hierarchy import linkage, leaves_list

from dashboard import plotting, util
from dashboard.util import filter_dataframe_dynamic, convert_list_string, convert_df
from dashboard.data_load import *

def tcga_page(sorf_df):
    st.title("sORF Transcriptome Atlas (TCGA)")

    df = filter_dataframe_dynamic(sorf_df, 'tcga_transcriptome')

    st.write(f'{df.shape[0]} sORF entries')

    with st.container():
        col1, col2, col3 = st.columns(3)
        with col1:
            xena_tau_df, xena_vtx_sum_df, xena_vtx_exp_df = load_xena_heatmap_data()

        with col2:
            values = st.slider(
                'Select Tissue Specificity Tau',
                0.0, 1.0, (.98, 1.0),
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
        st.write(df.shape)
        st.dataframe(df)

        st.download_button(
            label="Download as CSV",
            data=convert_df(df),
            file_name='sorf_atlas_selection.csv',
            mime='text/csv',
        )


def tissue_specific_page(sorf_df):
    st.title("Tissue specificity")

    filter_option = st.selectbox('Pre-filtered sORFs:', ('Ribo-Seq sORFs',
                                                        'Secreted',
                                                        'Secreted & Novel',
                                                        'Secreted & Conserved',
                                                        'Secreted & Conserved & Novel',
                                                        'Translated',
                                                        'Translated & Conserved',
                                                        'All sORFs'), index = 0, key='sorf_detail_filter2')
    
    df = util.filter_dataframe_preset(sorf_df, filter_option)

    st.write(f'{df.shape[0]} sORF entries')

    with st.container():
        col1, col2 = st.columns(2)
        with col1:
            xena_tau_df, xena_vtx_sum_df, xena_vtx_exp_df = load_xena_heatmap_data()
            expr_df = xena_tau_df[[c for c in xena_tau_df.columns if c not in ['tau']]]
            expr_vtx = set(expr_df[expr_df > .1].dropna(how='all').index)
            expr_vtx = list(expr_vtx.intersection(set(df.index)))
            xena_tau_df = xena_tau_df.groupby(level=0).median()
            xena_tau_df = xena_tau_df.reindex(index=expr_vtx)
            tissues = st.multiselect('Select Tissue(s):', 
                                     ['All Tissues'] + list(xena_tau_df.columns),
                                     key='tissue_specific_selection')
            if 'All Tissues' in tissues:
                tissues = list(xena_tau_df.columns)

        with col2:
            values = st.slider(
                'Select Tissue Specificity Tau',
                0.0, 1.0, (.85, 1.0),
                help='Higher values of [Tau](https://academic.oup.com/bib/article/18/2/205/2562739) indicate tissue specific expression of a given sORF')
            
            if tissues:
                subset_df = xena_tau_df[xena_tau_df['tau'].between(*values)]
                tissue_specific_vtx_ids = []

                top_tissues = []
                for i, row in subset_df.iterrows():
                    tissue_list = list(row[row > np.percentile(row, 99)].index.values)
                    for t in tissue_list:
                        if t in tissues:
                            tissue_specific_vtx_ids.append(i)
                    top_tissues.append(tissue_list)
                subset_df['tissues'] = top_tissues

                st.write(f'{len(tissue_specific_vtx_ids)}')
                st.write(subset_df.shape)

        #st.write(xena_vtx_sum_df.astype(float))

        if tissues and len(tissue_specific_vtx_ids) > 1:
            plot_df = xena_tau_df.drop(columns='tau')
            option, events = plotting.expression_atlas_heatmap_plot(tissue_specific_vtx_ids, plot_df.T)

            display_cols = ['vtx_id', 'screening_phase_id', 'screening_phase', 'orf_xrefs', 'protein_xrefs', 'gene_xrefs', 'transcript_xrefs', 'source']#, 'secreted_mean', 'translated_mean', 'isoform_of']
            value = st_echarts(option, height="1000px", events=events)
            if value:
                st.header('Selected sORF')
                st.dataframe(sorf_df[sorf_df['vtx_id'] == value][display_cols])
                
                fig = plotting.expression_vtx_boxplot(value, xena_vtx_exp_df)
                st.plotly_chart(fig, use_container_width=True)

            df = sorf_df[sorf_df['vtx_id'].isin(tissue_specific_vtx_ids)][display_cols]
            exp_df = subset_df.loc[tissue_specific_vtx_ids][['tissues', 'tau']].copy()
            df = df.merge(exp_df, left_index=True, right_index=True)
            st.header('Tissue specific sORFs')
            st.write(df.shape)
            st.dataframe(df)

            st.download_button(
                label="Download as CSV",
                data=convert_df(df),
                file_name='sorf_atlas_selection.csv',
                mime='text/csv',
            )
        
