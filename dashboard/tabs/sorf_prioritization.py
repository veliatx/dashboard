import os
import streamlit as st
import datetime
import pandas as pd
from dashboard.etl import CACHE_DIR

# Function to save the edited DataFrame
# def save_edited_df(edited_df):

    
def page(sorf_df):
    os.makedirs(CACHE_DIR / 'sorf_priority_tables', exist_ok=True)
    existing_versions = os.listdir(CACHE_DIR / 'sorf_priority_tables')
    existing_versions = sorted([i for i in existing_versions if i.endswith('_table.parq')])
    df = pd.read_parquet(CACHE_DIR / 'sorf_priority_tables' / existing_versions[-1])

    # Function to apply coloring based on value
    def color_code(value):
        if value == 'High':
            color = 'red'
        elif value == 'Medium':
            color = 'orange'
        elif value == 'Low':
            color = 'green'
        else:
            color = 'black'
        return f'color: {color};'

    # Apply coloring to DataFrame
    styled_df = df.style.applymap(lambda x: color_code(x) if isinstance(x, str) else '')

    # Session state to store the DataFrame
    if 'editable_df' not in st.session_state:
        st.session_state['editable_df'] = df


        # Display the styled DataFrame
    edited_df = st.data_editor(styled_df, num_rows='dynamic')
        
        # if submitted:
    st.session_state['editable_df'] = edited_df
    edited_df.to_parquet(CACHE_DIR / 'sorf_priority_tables' / f'{datetime.datetime.now().strftime("%Y%m%d-%H%M%S")}_table.parq')
        

            
            