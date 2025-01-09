"""
Main application file for the sORF dashboard.

This module contains the core functionality for the sORF web application built with Streamlit.
It defines the main page layout, navigation, and integrates various dashboard components.
"""

import streamlit as st
import streamlit.components.v1 as components
import pathlib

from dashboard.data_load import load_sorf_df_conformed
import dashboard.tabs.riboseq_atlas
import dashboard.tabs.sorf_explorer_table
import dashboard.tabs.expression_heatmap
import dashboard.tabs.sorf_prioritization
import dashboard.tabs.de_explorer
import dashboard.tabs.expression_atlas_summary

APP_NAME = 'sorf'


def genome_browser():
    """
    Renders an iframe containing the genome browser visualization.
    
    The genome browser is loaded from a local server and displayed in a scrollable iframe.
    """
    components.iframe("http://10.65.25.231:8080/velia_collections.html", height=1200, scrolling=True)


def main():
    """
    Main function that sets up and runs the Streamlit dashboard application.
    
    Configures the page layout, loads styles, initializes navigation tabs,
    and renders the appropriate content for each tab.
    """
    st.set_page_config(layout="wide")
    st.markdown(
        """ <style>iframe[title="streamlit_echarts.st_echarts"]{ height: 1000px !important } """,
        unsafe_allow_html=True
    )
    
    app_path = pathlib.Path(__file__.replace('sorf_app.py', ''))
    with open(app_path.joinpath("style.css")) as f:
        st.markdown(f"<style>{f.read()}</style>", unsafe_allow_html=True)
    
    # Define a dictionary with the page names and corresponding functions
    pages = {
        "sORF Details": dashboard.tabs.sorf_explorer_table.sorf_details,
        "sORF Genome Browser": genome_browser,
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
    # Initialize session state variables for riboseq atlas functionality
    if "experiments" not in st.session_state:
        st.session_state["experiments"] = []
    if "count" not in st.session_state:
        st.session_state["count"] = 0
    if "job_names" not in st.session_state:
        st.session_state["job_names"] = dashboard.tabs.riboseq_atlas.get_job_names_on_s3()
        
    main()
