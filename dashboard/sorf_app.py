
import streamlit as st
import streamlit.components.v1 as components

from dashboard.data_load import load_sorf_df_conformed
import dashboard.tabs.riboseq_atlas
import dashboard.tabs.sorf_explorer_table
import dashboard.tabs.expression_heatmap
import dashboard.tabs.sorf_prioritization
import dashboard.tabs.de_explorer
import dashboard.tabs.expression_atlas_summary
import pathlib

import gzip

APP_NAME = 'sorf'

def genome_browser():
    """
    """
    components.iframe("http://10.65.25.231:8080/velia_collections.html", height=1200, scrolling=True)


def main():
    
    st.set_page_config(layout="wide")
    st.markdown(""" <style>iframe[title="streamlit_echarts.st_echarts"]{ height: 1000px !important } """, unsafe_allow_html=True)
    app_path = pathlib.Path(__file__.replace('sorf_app.py', ''))
    with open(app_path.joinpath("style.css")) as f:
        st.markdown(f"<style>{f.read()}</style>", unsafe_allow_html=True)
    
    # Define a dictionary with the page names and corresponding functions
    pages = {
        "sORF Details": dashboard.tabs.sorf_explorer_table.sorf_details,
        #"sORF Transcriptome Atlas (TCGA)": dashboard.tabs.expression_heatmap.tcga_page,
        "sORF Genome Browser": genome_browser,
        #"sORF Ribo-seq Atlas": dashboard.tabs.riboseq_atlas.page,
        "DE Explorer": dashboard.tabs.de_explorer.de_page,
        "Expression Atlas": dashboard.tabs.expression_atlas_summary.expression_atlas_summary,
    }
    sorf_df = load_sorf_df_conformed()
    tab1, tab2, tab3, tab4 = st.tabs(list(pages.keys()))

    with tab1:
        dashboard.tabs.sorf_explorer_table.sorf_details(sorf_df)

    with tab2:
        genome_browser()
        
    with tab3:
       dashboard.tabs.de_explorer.de_page(sorf_df)

    with tab4:
        dashboard.tabs.expression_atlas_summary.expression_atlas_summary()

    
    
if __name__ == "__main__":
    # Config from riboseq atlas scripts
    if "experiments" not in st.session_state:
        st.session_state["experiments"] = []
    if "count" not in st.session_state:
        st.session_state["count"] = 0
    if "job_names" not in st.session_state:
        st.session_state["job_names"] = dashboard.tabs.riboseq_atlas.get_job_names_on_s3()
        
    main()
