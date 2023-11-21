from pandas.api.types import (
    is_categorical_dtype,
    is_datetime64_any_dtype,
    is_numeric_dtype,
    is_object_dtype,
    is_list_like
)

import pandas as pd
import streamlit as st
import sqlite3

def strip_ensembl_versions(transcript_ids):
    return [i.split('.')[0] if i.startswith('ENST') else i for i in transcript_ids]

def query_de_transcripts(transcript_id, 
                      db_address, 
                      log10padj_threshold = -2, 
                      minimum_expression = 2):
    # if isinstance(transcripts, str):
    #     transcripts = [transcripts]
    con = sqlite3.connect(db_address)
    query = """SELECT *
    FROM transcript_de
    WHERE transcript_de.transcript_id = '{0}'
    AND transcript_de.log10_padj <= {1}
    AND (transcript_de.case_mean >= {2} OR transcript_de.control_mean >= {2})
    """.format(transcript_id, log10padj_threshold, minimum_expression)
    return pd.read_sql(query, con)

def query_transcript_tpms(transcript_id_list, 
                          db_address):
    con = sqlite3.connect(db_address)
    query = """SELECT * FROM transcript_tpm
                WHERE transcript_tpm.transcript_id IN ({0});
                """.format(', '.join(transcript_id_list))
    return pd.read_sql(query, con)
    

@st.cache_data()
def convert_df(df):
    return df.to_csv().encode('utf-8')


def filter_dataframe(df: pd.DataFrame, key='details') -> pd.DataFrame:
    """
    Adds a UI on top of a dataframe to let viewers filter columns

    Args:
        df (pd.DataFrame): Original dataframe

    Returns:
        pd.DataFrame: Filtered dataframe
    """

    df = df.copy()

    # Try to convert datetimes into a standard format (datetime, no timezone)
    for col in df.columns:
        if is_object_dtype(df[col]):
            try:
                df[col] = pd.to_datetime(df[col], format="%d/%m/%Y")
            except Exception:
                pass

        if is_datetime64_any_dtype(df[col]):
            df[col] = df[col].dt.tz_localize(None)


    modification_container = st.container()

    with modification_container:
        to_filter_columns = st.multiselect("Filter dataframe on", df.columns)
        for column in to_filter_columns:
            left, right = st.columns((1, 20))
            # Treat columns with < 10 unique values as categorical
            if is_categorical_dtype(df[column]) or df[column].astype(str).nunique() < 10:
                user_cat_input = right.multiselect(
                    f"Values for {column}",
                    df[column].unique(),
                    default=list(df[column].unique()),
                    key=column
                )
                df = df[df[column].isin(user_cat_input)]
            elif is_numeric_dtype(df[column]):
                _min = float(df[column].min())
                _max = float(df[column].max())
                step = (_max - _min) / 100
                user_num_input = right.slider(
                    f"Values for {column}",
                    min_value=_min,
                    max_value=_max,
                    value=(_min, _max),
                    step=step,
                )
                df = df[df[column].between(*user_num_input)]
            elif is_datetime64_any_dtype(df[column]):
                user_date_input = right.date_input(
                    f"Values for {column}",
                    value=(
                        df[column].min(),
                        df[column].max(),
                    ),
                )
                if len(user_date_input) == 2:
                    user_date_input = tuple(map(pd.to_datetime, user_date_input))
                    start_date, end_date = user_date_input
                    df = df.loc[df[column].between(start_date, end_date)]
            
            #elif is_list_like(df[col]):
            #    df[col] = df[col].astype(str)
            else:
                user_text_input = right.text_input(
                    f"Substring or regex in {column}",
                )
                if user_text_input:
                    df = df[df[column].astype(str).str.contains(user_text_input)]

    return df


def ucsc_link(chrom, start, end):
    base_link = "http://genome.ucsc.edu/cgi-bin/hgTracks?"
    genome = "db=hg38"
    position = f"{chrom}:{start}-{end}"
    ucsc_link = f"{base_link}{genome}&position={position}"
    return ucsc_link


def convert_list_string(x):
    l = list(map(float, x.strip('][').split(',')))
    return l