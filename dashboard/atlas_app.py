"""
Main Streamlit application for the sORF tissue specificity atlas.

This module sets up the main application interface for exploring sORF tissue
expression patterns and specificity across different tissues and cell types.
"""

import pathlib
import streamlit as st
import streamlit.components.v1 as components

from dashboard.data_load import load_sorf_df_conformed
import dashboard.tabs.expression_heatmap

APP_NAME = 'atlas'


def main():
    """Initialize and run the main Streamlit application."""
    # Configure page layout
    st.set_page_config(layout="wide")
    
    # Set custom styling for charts
    st.markdown(
        """ <style>iframe[title="streamlit_echarts.st_echarts"]{ height: 1000px !important } """,
        unsafe_allow_html=True
    )
    
    # Load custom CSS
    app_path = pathlib.Path(__file__.replace('atlas_app.py', ''))
    with open(app_path.joinpath("style.css")) as f:
        st.markdown(f"<style>{f.read()}</style>", unsafe_allow_html=True)
    
    # Define available pages
    pages = {
        "sORF Tissue Specificity": dashboard.tabs.expression_heatmap.tissue_specific_page,
    }
    
    # Load data and display tissue specificity page
    sorf_df = load_sorf_df_conformed()
    dashboard.tabs.expression_heatmap.tissue_specific_page(sorf_df)


if __name__ == "__main__":
    main()
