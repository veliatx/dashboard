"""
Main Streamlit application for the sORF differential expression dashboard.

This module sets up the main application interface with multiple tabs for exploring
differential expression data, genome browser visualization, and Ribo-seq data.
"""

import streamlit as st
import streamlit.components.v1 as components

from dashboard.data_load import load_sorf_df_conformed
import dashboard.tabs.riboseq_atlas
import dashboard.tabs.sorf_explorer_table
import dashboard.tabs.expression_heatmap
import dashboard.tabs.sorf_prioritization
import dashboard.tabs.de_explorer

APP_NAME = 'de'


def genome_browser():
    """Display an embedded genome browser iframe."""
    components.iframe(
        "http://10.65.25.231:8080/velia_collections.html",
        height=1200,
        scrolling=True
    )


def main():
    """Initialize and run the main Streamlit application."""
    st.set_page_config(layout="wide")
    
    # Set custom styling
    st.markdown(
        """ <style>iframe[title="streamlit_echarts.st_echarts"]{ height: 1000px !important } """,
        unsafe_allow_html=True
    )
    with open("style.css") as f:
        st.markdown(f"<style>{f.read()}</style>", unsafe_allow_html=True)
    
    # Define available pages/tabs
    pages = {
        "DE Explorer": dashboard.tabs.de_explorer.de_page,
        "sORF Genome Browser": genome_browser,
        "sORF Ribo-seq Atlas": dashboard.tabs.riboseq_atlas.page,
    }
    
    # Load data and create tabs
    sorf_df = load_sorf_df_conformed()
    tab1, tab2, tab3 = st.tabs(list(pages.keys()))
    
    # Populate tabs with content
    with tab1:
        dashboard.tabs.de_explorer.de_page(sorf_df)
        
    with tab2:
        genome_browser()
        
    with tab3:
        dashboard.tabs.riboseq_atlas.page()


if __name__ == "__main__":
    # Initialize session state variables for Ribo-seq atlas
    if "experiments" not in st.session_state:
        st.session_state["experiments"] = []
    if "count" not in st.session_state:
        st.session_state["count"] = 0
    if "job_names" not in st.session_state:
        st.session_state["job_names"] = dashboard.tabs.riboseq_atlas.get_job_names_on_s3()
        
    main()
