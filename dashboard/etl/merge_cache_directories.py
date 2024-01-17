import pandas as pd
import os
import json
import click
import sqlite3
import pathlib
from dashboard.data_load import load_esmfold
import pickle

def copy_table(source_cursor, dest_cursor, table_name):
    # Copy table structure
    source_cursor.execute(f"SELECT sql FROM sqlite_master WHERE type='table' AND name='{table_name}'")
    create_table_sql = source_cursor.fetchone()[0]
    dest_cursor.execute(create_table_sql)

    # Copy data
    source_cursor.execute(f"SELECT * FROM {table_name}")
    rows = source_cursor.fetchall()
    dest_cursor.executemany(f"INSERT INTO {table_name} VALUES ({','.join(['?']*len(rows[0]))})", rows)

def merge_blastp_json(path1, path2):
    hits_p1, table_dict_p1 = pickle.load(open(path1, 'rb'))
    hits_p2, table_dict_p2 = pickle.load(open(path2, 'rb'))
    for k, v in hits_p2.items():
        if k not in hits_p1:
            hits_p1[k] = v
    for k, v in table_dict_p2.items():
        if k not in table_dict_p1:
            table_dict_p1[k] = v
    return hits_p1, table_dict_p1
# @click.command()
# @click.argument('EXISTING_CACHE_DIR', type=click.Path(exists=True, resolve_path=True))
# @click.argument('CACHE_DIR_TO_MERGE', type=click.Path(exists=True, resolve_path=True))
# @click.argument('MERGE_OUTPUT_FOLDER', type=click.Path(resolve_path=True))

def merge_cache_directories(EXISTING_CACHE_DIR, CACHE_DIR_TO_MERGE, MERGE_OUTPUT_FOLDER):
    EXISTING_CACHE_DIR = pathlib.Path(EXISTING_CACHE_DIR)
    CACHE_DIR_TO_MERGE = pathlib.Path(CACHE_DIR_TO_MERGE)
    MERGE_OUTPUT_FOLDER = pathlib.Path(MERGE_OUTPUT_FOLDER)
    MERGE_OUTPUT_FOLDER.joinpath('protein_data').mkdir(parents=True, exist_ok=True)
    pd.read_parquet(EXISTING_CACHE_DIR.joinpath('xena_metadata.parq')).to_parquet(MERGE_OUTPUT_FOLDER.joinpath('xena_metadata.parq'))
    #sorf dataframe
    sorf_df_existing = pd.read_parquet(EXISTING_CACHE_DIR.joinpath('sorf_df.parq'))
    sorf_df_merging = pd.read_parquet(CACHE_DIR_TO_MERGE.joinpath('sorf_df.parq'))
    pd.concat((sorf_df_existing, sorf_df_merging), join='outer').to_parquet(MERGE_OUTPUT_FOLDER.joinpath('sorf_df.parq'))

    #blastp json
    query_hits, json_table = merge_blastp_json(
        EXISTING_CACHE_DIR.joinpath('protein_data', 'blastp.results.pkl'),
        CACHE_DIR_TO_MERGE.joinpath('protein_data', 'blastp.results.pkl')
    )
    with open(MERGE_OUTPUT_FOLDER.joinpath('protein_data', 'blastp.results.pkl'), 'wb') as fwrite:
        pickle.dump((query_hits, json_table), fwrite)
    # Protein Features
    protein_scores_existing = pd.read_csv(os.path.join(EXISTING_CACHE_DIR, 'protein_data', 'sequence_features_scores.csv'), index_col=0)
    protein_strings_existing = pd.read_csv(os.path.join(EXISTING_CACHE_DIR, 'protein_data', 'sequence_features_strings.csv'), index_col=0)

    protein_scores_merging = pd.read_csv(os.path.join(CACHE_DIR_TO_MERGE, 'protein_data', 'sequence_features_scores.csv'), index_col=0)
    protein_strings_merging = pd.read_csv(os.path.join(CACHE_DIR_TO_MERGE, 'protein_data', 'sequence_features_strings.csv'), index_col=0)

    pd.concat((protein_scores_existing, protein_scores_merging), join='outer').to_csv(MERGE_OUTPUT_FOLDER.joinpath('protein_data', 'sequence_features_scores.csv'))
    pd.concat((protein_strings_existing, protein_strings_merging), join='outer').to_csv(MERGE_OUTPUT_FOLDER.joinpath('protein_data', 'sequence_features_strings.csv'))

    # ESMFold
    esmfold_existing = load_esmfold(EXISTING_CACHE_DIR.joinpath('protein_data', 'esmfold.jsonlines'))
    esmfold_merging = load_esmfold(CACHE_DIR_TO_MERGE.joinpath('protein_data', 'esmfold.jsonlines'))
    for k, v in esmfold_merging.items():
        if k not in esmfold_existing:
            esmfold_existing[k] = v
    with open(MERGE_OUTPUT_FOLDER.joinpath('protein_data', 'esmfold.jsonlines'), 'w') as fwrite:
        for v in esmfold_existing.values():
            fwrite.write(json.dumps(v))
            fwrite.write('\n')
            

    # Protein search tools
    isoform_check_existing = pd.read_parquet(os.path.join(EXISTING_CACHE_DIR, 'protein_data', 'isoforms_search.parq'))
    isoform_check_merging = pd.read_parquet(os.path.join(CACHE_DIR_TO_MERGE, 'protein_data', 'isoforms_search.parq'))

    pd.concat((isoform_check_existing, isoform_check_merging), join='outer').to_parquet(MERGE_OUTPUT_FOLDER.joinpath('protein_data', 'isoforms_search.parq'))

    nonsignal_blasts_existing = pd.read_parquet(os.path.join(EXISTING_CACHE_DIR, 'protein_data', 'nonsignal_seq_blast_tblastn.parq'))
    nonsignal_blasts_merging = pd.read_parquet(os.path.join(CACHE_DIR_TO_MERGE, 'protein_data', 'nonsignal_seq_blast_tblastn.parq'))

    pd.concat((nonsignal_blasts_existing, nonsignal_blasts_merging), join='outer').to_parquet(MERGE_OUTPUT_FOLDER.joinpath('protein_data', 'nonsignal_seq_blast_tblastn.parq'))

    # Mouse homology
    mouse_blastp_existing = pd.read_parquet(os.path.join(EXISTING_CACHE_DIR, 'protein_data', 'mouse_blastp.parq'))
    mouse_blastp_merging = pd.read_parquet(os.path.join(CACHE_DIR_TO_MERGE, 'protein_data', 'mouse_blastp.parq'))

    pd.concat((mouse_blastp_existing, mouse_blastp_merging), join='outer').to_parquet(MERGE_OUTPUT_FOLDER.joinpath('protein_data', 'mouse_blastp.parq'))

    mouse_tblastn_existing = pd.read_parquet(os.path.join(EXISTING_CACHE_DIR, 'protein_data', 'mouse_tblastn.parq'))
    mouse_tblastn_merging = pd.read_parquet(os.path.join(CACHE_DIR_TO_MERGE, 'protein_data', 'mouse_tblastn.parq'))

    pd.concat((mouse_tblastn_existing, mouse_tblastn_merging), join='outer').to_parquet(MERGE_OUTPUT_FOLDER.joinpath('protein_data', 'mouse_tblastn.parq'))



    # TCGA
    xena_exact_heatmap_data_existing = pickle.load(open(EXISTING_CACHE_DIR.joinpath('xena_exact_heatmap.pkl'), 'rb'))
    xena_exact_heatmap_data_merging = pickle.load(open(CACHE_DIR_TO_MERGE.joinpath('xena_exact_heatmap.pkl'), 'rb'))
    xena_exact_heatmap_data_merged = [None, None, None]
    xena_exact_heatmap_data_merged[0] = pd.concat((xena_exact_heatmap_data_existing[0], xena_exact_heatmap_data_merging[0]), axis=0, join='outer')
    xena_exact_heatmap_data_merged[1] = pd.concat((xena_exact_heatmap_data_existing[1], xena_exact_heatmap_data_merging[1]), axis=1, join='outer')
    xena_exact_heatmap_data_merged[2] = pd.concat((xena_exact_heatmap_data_existing[2], xena_exact_heatmap_data_merging[2]), axis=1, join='outer')
    pickle.dump(xena_exact_heatmap_data_merged, open(MERGE_OUTPUT_FOLDER.joinpath('xena_exact_heatmap.pkl'), 'wb'))

    xena_parq_merged = pd.concat((pd.read_parquet(EXISTING_CACHE_DIR.joinpath('xena_app.parq')), pd.read_parquet(CACHE_DIR_TO_MERGE.joinpath('xena_app.parq'))), axis=1, join='outer')
    xena_parq_merged.to_parquet(MERGE_OUTPUT_FOLDER.joinpath('xena_app.parq'))

    # Xena.db, (xena_de_tables_dict.pkl, CODE_de.parq) -> these actually shouldn't be needed since they're merged into xena.db
    # should be no merging autoimmune

    # Connect to the databases
    conn1 = sqlite3.connect(EXISTING_CACHE_DIR.joinpath('xena.db'))
    conn2 = sqlite3.connect(CACHE_DIR_TO_MERGE.joinpath('xena.db'))
    output_db = sqlite3.connect(MERGE_OUTPUT_FOLDER.joinpath('xena.db'))
    conn2 = sqlite3.connect('database2.db')
    # Create a cursor for the first database
    cursor1 = conn1.cursor()
    # Attach the second database to the first database connection
    cursor1.execute(f"ATTACH DATABASE '{CACHE_DIR_TO_MERGE.joinpath('xena.db')}' AS db2")
    merged_cursor = output_db.cursor()

    # Check for the existence of the table in the first database (optional)
    cursor1.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='transcript_de'")
    # Merge the data from the second database into the first
    # This assumes that the table structures are identical and there are no conflicting primary keys
    cursor1.execute("INSERT INTO transcript_de SELECT * FROM db2.transcript_de")
    cursor1.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='sample_metadata_de'")
    # Merge the data from the second database into the first
    # This assumes that the table structures are identical and there are no conflicting primary keys
    cursor1.execute("INSERT INTO sample_metadata_de SELECT * FROM db2.sample_metadata_de")

    conn1.commit()

    merged_cursor.execute(f"ATTACH DATABASE '{EXISTING_CACHE_DIR.joinpath('xena.db')}' AS merged_db")
    cursor1.execute("SELECT name FROM sqlite_master WHERE type='table'")
    tables = cursor1.fetchall()
    
    # Copy each table from the merged database to the new database
    for table in tables:
        table_name = table[0]
        copy_table(cursor1, merged_cursor, table_name)

    output_db.commit()

    # Commit the changes

    # Close the connections
    conn1.close()
    conn2.close()
    output_db.close()

    