"""
Utility script for merging two cache directories containing protein sequence data and analysis results.

This script merges data from two cache directories containing protein sequence analysis results,
including BLAST results, protein features, ESMFold predictions, and expression data. The merged
results are written to a new output directory.

The script handles merging of:
- Parquet files containing sORF data and metadata
- BLAST results in JSON/pickle format  
- Protein feature predictions (scores and string representations)
- ESMFold structure predictions
- Protein search tool results
- Mouse homology data
- TCGA/Xena expression data
- SQLite databases containing differential expression results
"""

import json
import os
import pathlib
import pickle
import sqlite3
from typing import Dict, List, Tuple, Any

import click
import pandas as pd

from dashboard.data_load import load_esmfold


def copy_table(source_cursor: sqlite3.Cursor,
               dest_cursor: sqlite3.Cursor,
               table_name: str) -> None:
    """Copy a table from one SQLite database to another.

    Args:
        source_cursor: Cursor for source database
        dest_cursor: Cursor for destination database
        table_name: Name of table to copy
    """
    # Copy table structure
    source_cursor.execute(
        f"SELECT sql FROM sqlite_master WHERE type='table' AND name='{table_name}'"
    )
    create_table_sql = source_cursor.fetchone()[0]
    dest_cursor.execute(create_table_sql)

    # Copy data
    source_cursor.execute(f"SELECT * FROM {table_name}")
    rows = source_cursor.fetchall()
    dest_cursor.executemany(
        f"INSERT INTO {table_name} VALUES ({','.join(['?']*len(rows[0]))})",
        rows
    )


def merge_blastp_json(path1: pathlib.Path,
                     path2: pathlib.Path) -> Tuple[Dict, Dict]:
    """Merge BLAST results from two pickle files.

    Args:
        path1: Path to first BLAST results pickle file
        path2: Path to second BLAST results pickle file

    Returns:
        Tuple containing:
        - Dictionary of merged BLAST hits
        - Dictionary of merged BLAST results tables
    """
    hits_p1, table_dict_p1 = pickle.load(open(path1, 'rb'))
    hits_p2, table_dict_p2 = pickle.load(open(path2, 'rb'))

    # Merge hits dictionaries
    for k, v in hits_p2.items():
        if k not in hits_p1:
            hits_p1[k] = v

    # Merge table dictionaries
    for k, v in table_dict_p2.items():
        if k not in table_dict_p1:
            table_dict_p1[k] = v

    return hits_p1, table_dict_p1


@click.command()
@click.argument('existing_cache_dir',
                type=click.Path(exists=True, resolve_path=True))
@click.argument('cache_dir_to_merge',
                type=click.Path(exists=True, resolve_path=True))
@click.argument('merge_output_folder',
                type=click.Path(resolve_path=True))
def merge_cache_directories(existing_cache_dir: str,
                          cache_dir_to_merge: str,
                          merge_output_folder: str) -> None:
    """Merge two cache directories containing protein sequence analysis results.

    Args:
        existing_cache_dir: Path to existing cache directory
        cache_dir_to_merge: Path to cache directory to merge in
        merge_output_folder: Path to output merged results
    """
    existing_cache_dir = pathlib.Path(existing_cache_dir)
    cache_dir_to_merge = pathlib.Path(cache_dir_to_merge)
    merge_output_folder = pathlib.Path(merge_output_folder)

    # Create output directory structure
    merge_output_folder.joinpath('protein_data').mkdir(parents=True,
                                                     exist_ok=True)

    # Copy metadata
    pd.read_parquet(existing_cache_dir.joinpath('xena_metadata.parq')).to_parquet(
        merge_output_folder.joinpath('xena_metadata.parq')
    )

    # Merge sORF dataframes
    sorf_df_existing = pd.read_parquet(
        existing_cache_dir.joinpath('sorf_df.parq')
    )
    sorf_df_merging = pd.read_parquet(
        cache_dir_to_merge.joinpath('sorf_df.parq')
    )
    pd.concat(
        (sorf_df_existing,
         sorf_df_merging.drop(
             sorf_df_existing.index.intersection(sorf_df_merging.index)
         )),
        join='outer'
    ).to_parquet(merge_output_folder.joinpath('sorf_df.parq'))

    # Merge BLAST results
    query_hits, json_table = merge_blastp_json(
        existing_cache_dir.joinpath('protein_data', 'blastp.results.pkl'),
        cache_dir_to_merge.joinpath('protein_data', 'blastp.results.pkl')
    )
    with open(merge_output_folder.joinpath('protein_data',
                                         'blastp.results.pkl'), 'wb') as fwrite:
        pickle.dump((query_hits, json_table), fwrite)

    # Merge protein features
    protein_scores_existing = pd.read_csv(
        os.path.join(existing_cache_dir, 'protein_data',
                    'sequence_features_scores.csv'),
        index_col=0
    )
    protein_strings_existing = pd.read_csv(
        os.path.join(existing_cache_dir, 'protein_data',
                    'sequence_features_strings.csv'),
        index_col=0
    )

    protein_scores_merging = pd.read_csv(
        os.path.join(cache_dir_to_merge, 'protein_data',
                    'sequence_features_scores.csv'),
        index_col=0
    )
    protein_strings_merging = pd.read_csv(
        os.path.join(cache_dir_to_merge, 'protein_data',
                    'sequence_features_strings.csv'),
        index_col=0
    )

    pd.concat(
        (protein_scores_existing,
         protein_scores_merging.drop(
             protein_scores_existing.index.intersection(
                 protein_scores_merging.index)
         )),
        join='outer'
    ).to_csv(merge_output_folder.joinpath('protein_data',
                                        'sequence_features_scores.csv'))

    pd.concat(
        (protein_strings_existing,
         protein_strings_merging.drop(
             protein_strings_existing.index.intersection(
                 protein_strings_merging.index)
         )),
        join='outer'
    ).to_csv(merge_output_folder.joinpath('protein_data',
                                        'sequence_features_strings.csv'))

    # Merge ESMFold results
    esmfold_existing = load_esmfold(
        existing_cache_dir.joinpath('protein_data', 'esmfold.jsonlines')
    )
    esmfold_merging = load_esmfold(
        cache_dir_to_merge.joinpath('protein_data', 'esmfold.jsonlines')
    )
    for k, v in esmfold_merging.items():
        if k not in esmfold_existing:
            esmfold_existing[k] = v
    with open(merge_output_folder.joinpath('protein_data',
                                         'esmfold.jsonlines'), 'w') as fwrite:
        for v in esmfold_existing.values():
            fwrite.write(json.dumps(v))
            fwrite.write('\n')

    # Merge protein search tool results
    isoform_check_existing = pd.read_parquet(
        os.path.join(existing_cache_dir, 'protein_data', 'isoforms_search.parq')
    )
    isoform_check_merging = pd.read_parquet(
        os.path.join(cache_dir_to_merge, 'protein_data', 'isoforms_search.parq')
    )

    pd.concat(
        (isoform_check_existing,
         isoform_check_merging.drop(
             isoform_check_existing.index.intersection(
                 isoform_check_merging.index)
         )),
        join='outer'
    ).to_parquet(merge_output_folder.joinpath('protein_data',
                                            'isoforms_search.parq'))

    nonsignal_blasts_existing = pd.read_parquet(
        os.path.join(existing_cache_dir, 'protein_data',
                    'nonsignal_seq_blast_tblastn.parq')
    )
    nonsignal_blasts_merging = pd.read_parquet(
        os.path.join(cache_dir_to_merge, 'protein_data',
                    'nonsignal_seq_blast_tblastn.parq')
    )

    pd.concat(
        (nonsignal_blasts_existing,
         nonsignal_blasts_merging.drop(
             nonsignal_blasts_existing.index.intersection(
                 nonsignal_blasts_merging.index)
         )),
        join='outer'
    ).to_parquet(merge_output_folder.joinpath(
        'protein_data', 'nonsignal_seq_blast_tblastn.parq'))

    # Merge mouse homology data
    mouse_blastp_existing = pd.read_parquet(
        os.path.join(existing_cache_dir, 'protein_data', 'mouse_blastp.parq')
    )
    mouse_blastp_merging = pd.read_parquet(
        os.path.join(cache_dir_to_merge, 'protein_data', 'mouse_blastp.parq')
    )

    pd.concat(
        (mouse_blastp_existing,
         mouse_blastp_merging.drop(
             mouse_blastp_existing.index.intersection(
                 mouse_blastp_merging.index),
             axis=0
         )),
        join='outer'
    ).to_parquet(merge_output_folder.joinpath('protein_data',
                                            'mouse_blastp.parq'))

    mouse_tblastn_existing = pd.read_parquet(
        os.path.join(existing_cache_dir, 'protein_data', 'mouse_tblastn.parq')
    )
    mouse_tblastn_merging = pd.read_parquet(
        os.path.join(cache_dir_to_merge, 'protein_data', 'mouse_tblastn.parq')
    )

    pd.concat(
        (mouse_tblastn_existing,
         mouse_tblastn_merging.drop(
             mouse_tblastn_existing.index.intersection(
                 mouse_tblastn_merging.index),
             axis=0
         )),
        join='outer'
    ).to_parquet(merge_output_folder.joinpath('protein_data',
                                            'mouse_tblastn.parq'))

    # Merge TCGA data
    xena_exact_heatmap_data_existing = pickle.load(
        open(existing_cache_dir.joinpath('xena_exact_heatmap.pkl'), 'rb')
    )
    xena_exact_heatmap_data_merging = pickle.load(
        open(cache_dir_to_merge.joinpath('xena_exact_heatmap.pkl'), 'rb')
    )
    xena_exact_heatmap_data_merged = [None, None, None]
    xena_exact_heatmap_data_merged[0] = pd.concat(
        (xena_exact_heatmap_data_existing[0],
         xena_exact_heatmap_data_merging[0]),
        axis=0,
        join='outer'
    )
    xena_exact_heatmap_data_merged[1] = pd.concat(
        (xena_exact_heatmap_data_existing[1],
         xena_exact_heatmap_data_merging[1]),
        axis=1,
        join='outer'
    )
    xena_exact_heatmap_data_merged[2] = pd.concat(
        (xena_exact_heatmap_data_existing[2],
         xena_exact_heatmap_data_merging[2]),
        axis=1,
        join='outer'
    )
    pickle.dump(xena_exact_heatmap_data_merged,
                open(merge_output_folder.joinpath('xena_exact_heatmap.pkl'),
                     'wb'))

    xena_parq_merged = pd.concat(
        (pd.read_parquet(existing_cache_dir.joinpath('xena_app.parq')),
         pd.read_parquet(cache_dir_to_merge.joinpath('xena_app.parq'))),
        axis=1,
        join='outer'
    )
    xena_parq_merged = xena_parq_merged.T.drop_duplicates().T
    xena_parq_merged.to_parquet(merge_output_folder.joinpath('xena_app.parq'))

    # Merge SQLite databases
    conn1 = sqlite3.connect(existing_cache_dir.joinpath('xena.db'))
    conn2 = sqlite3.connect(cache_dir_to_merge.joinpath('xena.db'))
    output_db = sqlite3.connect(merge_output_folder.joinpath('xena.db'))

    cursor1 = conn1.cursor()
    cursor1.execute(
        f"ATTACH DATABASE '{cache_dir_to_merge.joinpath('xena.db')}' AS db2"
    )
    merged_cursor = output_db.cursor()

    # Merge transcript differential expression data
    cmd = """INSERT INTO transcript_de 
            SELECT * FROM db2.transcript_de 
            WHERE NOT EXISTS (
                SELECT 1 
                FROM transcript_de 
                WHERE transcript_de.transcript_id = db2.transcript_de.transcript_id
            )"""
    cursor1.execute(cmd)
    conn1.commit()

    merged_cursor.execute(
        f"ATTACH DATABASE '{existing_cache_dir.joinpath('xena.db')}' AS merged_db"
    )
    cursor1.execute("SELECT name FROM sqlite_master WHERE type='table'")
    tables = cursor1.fetchall()

    # Copy tables to new database
    for table in tables:
        table_name = table[0]
        copy_table(cursor1, merged_cursor, table_name)

    output_db.commit()

    # Close database connections
    conn1.close()
    conn2.close()
    output_db.close()