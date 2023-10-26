import json
import jsonlines
import os
import pickle
import py3Dmol

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px
import streamlit as st
import streamlit.components.v1 as components
import streamlit_scrollable_textbox as stx

from collections import defaultdict
from streamlit_plotly_events import plotly_events
from tqdm import tqdm

from streamlit_echarts import st_echarts
from scipy.cluster.hierarchy import linkage, leaves_list

from dashboard import plotting, description
from dashboard.util import filter_dataframe, convert_list_string, convert_df
from dashboard.etl.sorf_query import load_jsonlines_table
from dashboard.data_load import *
import dashboard.tabs.riboseq_atlas
import dashboard.tabs.sorf_explorer_table
import dashboard.tabs.expression_heatmap

from veliadb import base
from veliadb.base import Orf, Protein, ProteinXref
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
        # "sORF Transcriptome Atlas (Autoimmune)": dashboard.tabs.expression_heatmap.autoimmune_page,
        "sORF Genome Browser": genome_browser,
        "sORF Ribo-seq Atlas": dashboard.tabs.riboseq_atlas.page
    }
    sorf_df = load_sorf_df_conformed()
    tcga_data = load_xena_tcga_gtex_target(sorf_df)
    autoimmune_data = load_autoimmune_atlas(sorf_df)
    tab1, tab2, tab3, tab4 = st.tabs(list(pages.keys()))

    with tab1:
        dashboard.tabs.sorf_explorer_table.sorf_details(sorf_df, tcga_data, autoimmune_data)

    with tab2:
        dashboard.tabs.expression_heatmap.tcga_page(sorf_df, tcga_data)
    # with tab3:
        # dashboard.tabs.expression_heatmap.autoimmune_page(sorf_df, autoimmune_data)
    with tab3:
        genome_browser()
        
    with tab4:
        dashboard.tabs.riboseq_atlas.page()
        
if __name__ == "__main__":
    # Config from riboseq atlas scripts
    if "experiments" not in st.session_state:
        st.session_state["experiments"] = []
    if "count" not in st.session_state:
        st.session_state["count"] = 0
    if "job_names" not in st.session_state:
        st.session_state["job_names"] = tabs.riboseq_atlas.get_job_names_on_s3()
        
    main()
