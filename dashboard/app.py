
import streamlit as st
import streamlit.components.v1 as components

from streamlit_plotly_events import plotly_events

from streamlit_echarts import st_echarts
from scipy.cluster.hierarchy import linkage, leaves_list

from dashboard.util import filter_dataframe, convert_list_string, convert_df
from dashboard.etl.sorf_query import load_jsonlines_table
from dashboard.data_load import *
import dashboard.tabs.riboseq_atlas
import dashboard.tabs.sorf_explorer_table
import dashboard.tabs.expression_heatmap
import dashboard.tabs.sorf_prioritization
import dashboard.tabs.de_explorer

import gzip

def genome_browser():
    """
    """
    components.iframe("http://10.65.25.231:8080/velia_collections.html", height=1200, scrolling=True)


def main():
    
    st.set_page_config(layout="wide")
    st.markdown(""" <style>iframe[title="streamlit_echarts.st_echarts"]{ height: 1000px !important } """, unsafe_allow_html=True)
    with open("style.css") as f:
        st.markdown(f"<style>{f.read()}</style>", unsafe_allow_html=True)
    
    # Define a dictionary with the page names and corresponding functions
    pages = {
        "sORF Details": dashboard.tabs.sorf_explorer_table.sorf_details,
        "sORF Transcriptome Atla (TCGA)": dashboard.tabs.expression_heatmap.tcga_page,
        "sORF Genome Browser": genome_browser,
        "sORF Ribo-seq Atlas": dashboard.tabs.riboseq_atlas.page,
       "DE Explorer": dashboard.tabs.de_explorer.de_page
    }
    sorf_df = load_sorf_df_conformed()
    tab1, tab2, tab3, tab4, tab5 = st.tabs(list(pages.keys()))

    with tab1:
        dashboard.tabs.sorf_explorer_table.sorf_details(sorf_df)

    with tab2:
        dashboard.tabs.expression_heatmap.tcga_page(sorf_df)

    with tab3:
        genome_browser()
        
    with tab4:
        dashboard.tabs.riboseq_atlas.page()
    
    with tab5:
        dashboard.tabs.de_explorer.de_page(sorf_df)
    
    # with tab5:
        # dashboard.tabs.sorf_prioritization.page(sorf_df)
        
if __name__ == "__main__":
    # Config from riboseq atlas scripts
    if "experiments" not in st.session_state:
        st.session_state["experiments"] = []
    if "count" not in st.session_state:
        st.session_state["count"] = 0
    if "job_names" not in st.session_state:
        st.session_state["job_names"] = tabs.riboseq_atlas.get_job_names_on_s3()
        
    main()
