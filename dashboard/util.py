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


def strip_ensembl_versions(transcript_ids: list) -> list:
    """Remove version numbers from Ensembl transcript IDs.
    
    Args:
        transcript_ids: List of transcript IDs
        
    Returns:
        List of transcript IDs with version numbers removed for Ensembl IDs
    """
    return [i.split('.')[0] if i.startswith('ENST') else i for i in transcript_ids]


def query_de_transcripts(transcript_id: str,
                        db_address: str,
                        log10padj_threshold: float = -2,
                        minimum_expression: float = 2) -> pd.DataFrame:
    """Query differential expression data for a transcript.
    
    Args:
        transcript_id: Transcript identifier
        db_address: Path to SQLite database
        log10padj_threshold: Maximum log10 adjusted p-value threshold
        minimum_expression: Minimum expression threshold in TPM
        
    Returns:
        DataFrame containing differential expression results
    """
    con = sqlite3.connect(db_address)
    query = f"""SELECT *
    FROM transcript_de
    WHERE transcript_de.transcript_id = '{transcript_id}'
    AND transcript_de.log10_padj <= {log10padj_threshold}
    AND (transcript_de.case_mean >= {minimum_expression} OR transcript_de.control_mean >= {minimum_expression})
    """
    return pd.read_sql(query, con)


def query_transcript_tpms(transcript_id_list: list,
                        db_address: str) -> pd.DataFrame:
    """Query TPM values for a list of transcripts.
    
    Args:
        transcript_id_list: List of transcript identifiers
        db_address: Path to SQLite database
        
    Returns:
        DataFrame containing TPM values for requested transcripts
    """
    con = sqlite3.connect(db_address)
    query = f"""SELECT * FROM transcript_tpm
                WHERE transcript_tpm.transcript_id IN ({', '.join(transcript_id_list)});
             """
    return pd.read_sql(query, con)


@st.cache_data()
def convert_df(df: pd.DataFrame) -> bytes:
    """Convert DataFrame to CSV bytes.
    
    Args:
        df: Input DataFrame
        
    Returns:
        UTF-8 encoded CSV bytes
    """
    return df.to_csv().encode('utf-8')


def filter_dataframe_dynamic(df: pd.DataFrame, key: str = 'details') -> pd.DataFrame:
    """Add interactive filtering UI on top of a DataFrame.

    Args:
        df: Original DataFrame
        key: Unique key for Streamlit widget state

    Returns:
        Filtered DataFrame based on user selections
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
            else:
                user_text_input = right.text_input(
                    f"Substring or regex in {column}",
                )
                if user_text_input:
                    df = df[df[column].astype(str).str.contains(user_text_input)]

    return df


def ucsc_link(chrom: str, start: int, end: int) -> str:
    """Generate UCSC Genome Browser link for a genomic region.
    
    Args:
        chrom: Chromosome name
        start: Start position
        end: End position
        
    Returns:
        URL to UCSC Genome Browser showing specified region
    """
    base_link = "http://genome.ucsc.edu/cgi-bin/hgTracks?"
    genome = "db=hg38"
    position = f"{chrom}:{start}-{end}"
    ucsc_link = f"{base_link}{genome}&position={position}"
    return ucsc_link


def convert_list_string(x: str) -> list:
    """Convert string representation of list to actual list of floats.
    
    Args:
        x: String representation of list (e.g. '[1.0, 2.0, 3.0]')
        
    Returns:
        List of float values
    """
    return list(map(float, x.strip('][').split(',')))


def filter_dataframe_preset(sorf_df: pd.DataFrame, filter_option: str) -> pd.DataFrame:
    """Filter sORF DataFrame based on preset categories.

    Args:
        sorf_df: Original DataFrame containing sORF data
        filter_option: Key specifying preset filtering strategy

    Returns:
        Filtered DataFrame based on selected preset category
    """
    view_cols = list(sorf_df.columns)
    df = sorf_df[view_cols].copy()
    
    signal_cols = ['SignalP 4.1_cut', 'SignalP 5b_cut', 'SignalP 6slow_cut', 'Deepsig_cut']
    conservation_cols = ['tblastn_align_identity', 'blastp_align_identity', 'nonsig_blastp_align_identity', 'nonsig_tblastn_align_identity']
    isoform_cols = ['swissprot_isoform', 'ensembl_isoform', 'refseq_isoform']
    conservation_threshold = 50
    exist_on_transcript = df['transcripts_exact'].apply(len).astype('bool')
    measured_secreted_or_predicted_secreted = df['secreted_hibit'] | (df[signal_cols] > -1).any(axis=1)
    is_not_isoform = df[isoform_cols].apply(lambda x: [not i=='None' for i in x]).max(axis=1)==0
    
    if filter_option == 'Secreted & Transmembrane':
        df = df[df['Ribo-Seq sORF']]
    elif filter_option == 'Secreted':
        df = df[(df['Ribo-Seq sORF']) & measured_secreted_or_predicted_secreted]
    elif filter_option == 'Secreted & Novel':
        df = df[(df['Ribo-Seq sORF']) & measured_secreted_or_predicted_secreted & is_not_isoform]
    elif filter_option == 'Secreted & Conserved':
        df = df[(df['Ribo-Seq sORF']) & 
                measured_secreted_or_predicted_secreted & 
                (df[conservation_cols] > conservation_threshold).any(axis=1)]
    elif filter_option == 'Secreted & Conserved & Novel':
        df = df[(df['Ribo-Seq sORF']) & 
                measured_secreted_or_predicted_secreted & 
                is_not_isoform & 
                (df[conservation_cols] > conservation_threshold).any(axis=1)]
    elif filter_option == 'Transmembrane':
        df = df[(df['Ribo-Seq sORF']) & 
                ~measured_secreted_or_predicted_secreted & 
                (df['DeepTMHMM_prediction'])]
    elif filter_option == 'Transmembrane & Conserved':
        df = df[(df['Ribo-Seq sORF']) & 
                ~measured_secreted_or_predicted_secreted & 
                (df['DeepTMHMM_prediction']) & 
                (df[conservation_cols] > conservation_threshold).any(axis=1)]
    elif filter_option == 'Transmembrane & Conserved & Novel':
        df = df[(df['Ribo-Seq sORF']) & 
                ~measured_secreted_or_predicted_secreted & 
                (df['DeepTMHMM_prediction']) & 
                is_not_isoform & 
                (df[conservation_cols] > conservation_threshold).any(axis=1)]
    elif filter_option == 'All sORFs':
        pass

    return df
