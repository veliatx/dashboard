import plotly.express as px
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


def filter_weak_expression(xena_tau_df, filtered_sorf_df):
    """
    """
    expr_df = xena_tau_df[[c for c in xena_tau_df.columns if c not in ['tau']]]
    expr_vtx = set(expr_df[expr_df > .5].dropna(how='all').index)
    expr_vtx = list(expr_vtx.intersection(set(filtered_sorf_df.index)))
    xena_tau_df = xena_tau_df.groupby(level=0).median()
    xena_tau_df = xena_tau_df.reindex(index=expr_vtx)

    return xena_tau_df


def assign_top_tissues(xena_tau_df):
    """
    """

    top_tissues = []
    tissue_cnt_dict = defaultdict(int)
    for i, row in xena_tau_df.iterrows():
        tissue_list = list(row[row > np.percentile(row, 99)].index.values)
        top_tissues.append(tissue_list)
        for t in tissue_list:
            tissue_cnt_dict[t] += 1

    xena_tau_df['tissues'] = top_tissues
    tissue_cnt_df = pd.DataFrame.from_dict(tissue_cnt_dict, 
                                           orient='index',
                                           columns=['count'])
    tissue_cnt_df.index.name = 'tissue'
    tissue_cnt_df.drop('tau', inplace=True)
    tissue_cnt_df.reset_index(inplace=True)
    tissue_cnt_df.sort_values(by='count', inplace=True, ascending=False)
    
    return xena_tau_df, tissue_cnt_df


def tissue_specific_page(sorf_df):
    """
    """
    st.title("Tissue specificity")

    with st.container():
        col1, col2 = st.columns(2)
        with col1:
            filter_option = st.selectbox('Select sORF Collection', ('Ribo-Seq sORFs',
                                                        'Secreted',
                                                        'Secreted & Novel',
                                                        'Secreted & Conserved',
                                                        'Secreted & Conserved & Novel',
                                                        'Translated',
                                                        'Translated & Conserved',
                                                        'All sORFs'), index = 0, key='sorf_detail_filter2')

        with col2:
            values = st.slider(
                'Select Tissue Specificity Tau',
                0.0, 1.0, (.85, 1.0),
                help='Higher values of [Tau](https://academic.oup.com/bib/article/18/2/205/2562739) indicate tissue specific expression of a given sORF')
    
            
        filtered_sorf_df = util.filter_dataframe_preset(sorf_df, filter_option)

        xena_tau_df, xena_vtx_sum_df, xena_vtx_exp_df = load_xena_heatmap_data()
        xena_tau_df = filter_weak_expression(xena_tau_df, filtered_sorf_df)
        xena_tau_df = xena_tau_df[xena_tau_df['tau'].between(*values)].copy()
        xena_tau_df, tissue_cnt_df = assign_top_tissues(xena_tau_df)

        st.write(f'{filtered_sorf_df.shape[0]} total sORFs')
        st.write(f'{xena_tau_df.shape[0]} tissue specific sORFs')
        fig = px.bar(tissue_cnt_df, x='tissue', y='count')
        st.plotly_chart(fig, use_container_width=True)

        selected_tissues = st.multiselect('Select Tissue(s):', 
                                     ['All Tissues'] + list(xena_tau_df.columns),
                                     key='tissue_specific_selection')
        
        if 'All Tissues' in selected_tissues:
            selected_tissues = list(xena_tau_df.columns)
        
        if selected_tissues:
            tissue_specific_vtx_ids = []

            for i, row in xena_tau_df.iterrows():
                for t in row.tissues:
                    if t in selected_tissues:
                        tissue_specific_vtx_ids.append(i)

            st.write(f'{len(tissue_specific_vtx_ids)}')
            st.write(xena_tau_df.shape)

        #st.write(xena_vtx_sum_df.astype(float))

        if selected_tissues and len(tissue_specific_vtx_ids) > 1:
            plot_df = xena_tau_df.drop(columns=['tau', 'tissues'])
            option, events = plotting.expression_atlas_heatmap_plot(tissue_specific_vtx_ids, plot_df.T)

            display_cols = ['vtx_id', 'screening_phase_id', 'screening_phase', 'orf_xrefs', 'protein_xrefs', 'gene_xrefs', 'transcript_xrefs', 'source']#, 'secreted_mean', 'translated_mean', 'isoform_of']
            value = st_echarts(option, height="1000px", events=events)
            if value:
                st.header('Selected sORF')
                st.dataframe(sorf_df[sorf_df['vtx_id'] == value][display_cols])
                
                fig = plotting.expression_vtx_boxplot(value, xena_vtx_exp_df)
                st.plotly_chart(fig, use_container_width=True)

            df = sorf_df[sorf_df['vtx_id'].isin(tissue_specific_vtx_ids)][display_cols]
            exp_df = xena_tau_df.loc[tissue_specific_vtx_ids][['tissues', 'tau']].copy()
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
        
