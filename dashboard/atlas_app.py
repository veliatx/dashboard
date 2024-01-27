
import streamlit as st
import streamlit.components.v1 as components

from dashboard.data_load import load_sorf_df_conformed
import dashboard.tabs.riboseq_atlas
import dashboard.tabs.sorf_explorer_table
import dashboard.tabs.expression_heatmap
import dashboard.tabs.sorf_prioritization
import dashboard.tabs.de_explorer
import pathlib

import gzip

APP_NAME = 'atlas'


def main():
    
    st.set_page_config(layout="wide")
    st.markdown(""" <style>iframe[title="streamlit_echarts.st_echarts"]{ height: 1000px !important } """, unsafe_allow_html=True)
    app_path = pathlib.Path(__file__.replace('atlas_app.py', ''))
    with open(app_path.joinpath("style.css")) as f:
        st.markdown(f"<style>{f.read()}</style>", unsafe_allow_html=True)
    
    pages = {
        "sORF Tissue Specificity": dashboard.tabs.expression_heatmap.tissue_specific_page,
    }
    sorf_df = load_sorf_df_conformed()
    #tab1, tab2, tab3, tab4, tab5 = st.tabs(list(pages.keys()))

    dashboard.tabs.expression_heatmap.tissue_specific_page(sorf_df)
    
    
if __name__ == "__main__":
        
    main()
