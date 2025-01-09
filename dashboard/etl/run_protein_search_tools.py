"""
Script to run protein sequence search tools and process results.

This script performs sequence similarity searches using MMSeqs2 and BLAST against various databases:
- SwissProt, Ensembl, and RefSeq protein databases for isoform detection
- Mouse protein and RNA databases for homology detection
- Searches with and without signal peptides

Results are saved as parquet files containing alignment statistics and matches.
"""

import gzip
import pathlib
import shlex
import subprocess
from typing import List, Optional

import click
import pandas as pd

from dashboard import plotting, description
from dashboard.etl.sorf_query import load_jsonlines_table


@click.command()
@click.argument('input_vtx', type=click.Path(exists=True, resolve_path=True))
@click.argument('output_prefix', type=click.Path(resolve_path=True))
@click.option("--number_threads", help="Number of threads to use for parallel computing.",
              default=8, type=int)
def run_protein_search_tools(input_vtx: str,
                           output_prefix: str,
                           number_threads: int) -> None:
    """Run protein sequence search tools and process results.

    Args:
        input_vtx: Path to input FASTA file containing protein sequences
        output_prefix: Path prefix for output files
        number_threads: Number of threads to use for parallel processing

    Raises:
        ValueError: If required reference files are missing
    """
    input_vtx = pathlib.Path(input_vtx)
    output_prefix = pathlib.Path(output_prefix)

    # Define reference file paths
    REFSEQ_FASTA_FILE = output_prefix.joinpath('GRCh38_latest_protein.faa')
    ENSEMBL_FASTA_FILE = output_prefix.joinpath('swissprot_proteins.fa')
    SWISSPROT_FASTA_FILE = output_prefix.joinpath('ensembl_proteins.fa')
    VTX_WITHOUT_SIGNAL_SEQ = output_prefix.joinpath('vtx_aa_seq_signal_sequence_removed.fa')

    # Check for missing reference files
    any_missing: List[pathlib.Path] = []
    for f in [REFSEQ_FASTA_FILE, ENSEMBL_FASTA_FILE,
              SWISSPROT_FASTA_FILE, VTX_WITHOUT_SIGNAL_SEQ]:
        if not f.exists():
            any_missing.append(f)
    if any_missing:
        raise ValueError(
            "Reference files are missing. Please check that update_cache.py has completed, "
            f"or create these manually in the output_prefix directory: {any_missing}"
        )

    # Run SwissProt isoform search
    swissprot_isoform_df = _run_mmseqs_search(
        input_vtx, SWISSPROT_FASTA_FILE, output_prefix
    )
    swissprot_isoform_df.rename(columns={'target': 'swissprot_isoform'}, inplace=True)

    # Run Ensembl isoform search
    ensembl_isoform_df = _run_mmseqs_search(
        input_vtx, ENSEMBL_FASTA_FILE, output_prefix
    )
    ensembl_isoform_df.rename(columns={'target': 'ensembl_isoform'}, inplace=True)

    # Run RefSeq isoform search
    refseq_isoform_df = _run_mmseqs_search(
        input_vtx, REFSEQ_FASTA_FILE, output_prefix
    )
    refseq_isoform_df.rename(columns={'target': 'refseq_isoform'}, inplace=True)

    # Run BLAST searches against mouse databases
    _run_mouse_blast_searches(input_vtx, output_prefix, number_threads)

    # Merge and save isoform results
    isoform_df = pd.DataFrame([
        swissprot_isoform_df['swissprot_isoform'],
        ensembl_isoform_df['ensembl_isoform'],
        refseq_isoform_df['refseq_isoform']
    ]).T
    isoform_df.to_parquet(pathlib.Path(output_prefix).joinpath('isoforms_search.parq'))

    # Run BLAST searches for sequences with signal peptides removed
    _run_signal_peptide_blast_searches(
        VTX_WITHOUT_SIGNAL_SEQ, output_prefix, number_threads
    )


def _run_mmseqs_search(input_file: pathlib.Path,
                      target_file: pathlib.Path,
                      output_prefix: pathlib.Path) -> pd.DataFrame:
    """Run MMSeqs2 search and process results.

    Args:
        input_file: Path to query sequences
        target_file: Path to target database
        output_prefix: Output directory path

    Returns:
        DataFrame containing filtered and processed search results
    """
    query_db = f'/root/{input_file.name}'
    target_db = f'/root/{target_file.name}'
    output_file = f'/root/{input_file.stem}_{target_file.stem}_alignments.m8'
    options = ('--format-output query,target,evalue,qstart,qend,qlen,'
              'qcov,gapopen,pident,fident,alnlen,raw,bits')

    base_cmd = f"docker run --rm -v {output_prefix}:/root soedinglab/mmseqs2 mmseqs easy-search"
    full_cmd = f'{base_cmd} {options} {query_db} {target_db} {output_file} /root/tmp'
    subprocess.run(shlex.split(full_cmd))

    # Process results
    isoform_df = pd.read_csv(
        output_prefix.joinpath(f'{input_file.stem}_{target_file.stem}_alignments.m8'),
        sep='\t',
        names=options.split()[1].split(',')
    )
    isoform_df = isoform_df[isoform_df['pident'] >= 100]
    isoform_df.drop_duplicates(inplace=True)
    return isoform_df.groupby('query').aggregate(list)


def _run_mouse_blast_searches(input_file: pathlib.Path,
                            output_prefix: pathlib.Path,
                            number_threads: int) -> None:
    """Run BLAST searches against mouse protein and RNA databases.

    Args:
        input_file: Path to query sequences
        output_prefix: Output directory path
        number_threads: Number of threads for parallel processing
    """
    blast_db_path = pathlib.Path('/efs/databases/blast')
    
    # Run protein BLAST
    _run_blast_search(
        input_file=input_file,
        output_prefix=output_prefix,
        blast_db_path=blast_db_path,
        db_name='mouse.protein.genbank.faa',
        output_name='mouse_blastp.parq',
        program='blastp',
        number_threads=number_threads
    )

    # Run translated BLAST
    _run_blast_search(
        input_file=input_file,
        output_prefix=output_prefix,
        blast_db_path=blast_db_path,
        db_name='mouse.rna.fna',
        output_name='mouse_tblastn.parq',
        program='tblastn',
        number_threads=number_threads
    )

    # Run fast BLAST for JSON output
    cmd = (f'docker run --rm -v {output_prefix}:/data -v /efs/databases/blast:/db '
           f'ncbi/blast blastp -task blastp-fast -outfmt 15 -db /db/mouse.protein.faa '
           f'-query /data/{input_file.name} -max_target_seqs 20 -num_threads 1 '
           f'-out /data/blastp.results.json')
    subprocess.run(shlex.split(cmd))


def _run_signal_peptide_blast_searches(input_file: pathlib.Path,
                                     output_prefix: pathlib.Path,
                                     number_threads: int) -> None:
    """Run BLAST searches for sequences with signal peptides removed.

    Args:
        input_file: Path to sequences without signal peptides
        output_prefix: Output directory path
        number_threads: Number of threads for parallel processing
    """
    blast_db_path = pathlib.Path('/efs/databases/blast')

    # Run protein BLAST
    blastp_df = _run_blast_search(
        input_file=input_file,
        output_prefix=output_prefix,
        blast_db_path=blast_db_path,
        db_name='mouse.protein.genbank.faa',
        output_name=None,
        program='blastp',
        number_threads=number_threads,
        prefix='nonsig_blastp'
    )

    # Run translated BLAST
    tblastn_df = _run_blast_search(
        input_file=input_file,
        output_prefix=output_prefix,
        blast_db_path=blast_db_path,
        db_name='mouse.rna.fna',
        output_name=None,
        program='tblastn',
        number_threads=number_threads,
        prefix='nonsig_tblastn'
    )

    # Merge and save results
    output_df = blastp_df.merge(
        tblastn_df, left_index=True, right_index=True, how='outer'
    )
    output_df.to_parquet(output_prefix.joinpath('nonsignal_seq_blast_tblastn.parq'))


def _run_blast_search(input_file: pathlib.Path,
                     output_prefix: pathlib.Path,
                     blast_db_path: pathlib.Path,
                     db_name: str,
                     output_name: Optional[str],
                     program: str,
                     number_threads: int,
                     prefix: str = '') -> pd.DataFrame:
    """Run a BLAST search with specified parameters.

    Args:
        input_file: Path to query sequences
        output_prefix: Output directory path
        blast_db_path: Path to BLAST databases
        db_name: Name of target database
        output_name: Name for output file (optional)
        program: BLAST program to use
        number_threads: Number of threads for parallel processing
        prefix: Prefix for column names

    Returns:
        DataFrame containing filtered and processed search results
    """
    output_fmt = '6 qaccver saccver stitle bitscore qcovs length pident gaps evalue'
    options = f'-outfmt "{output_fmt}" -num_threads {number_threads}'
    query = f'-query /blast/data/{input_file.name}'
    output = f'-out /blast/data/{input_file.stem}.{program}.out'
    
    base_cmd = (f'docker run --rm -v {blast_db_path}:/blast/blastdb '
                f'-v {output_prefix}:/blast/data ncbi/blast')
    full_cmd = f'{base_cmd} {program} {options} -db {db_name} {query} {output}'
    subprocess.run(shlex.split(full_cmd))

    # Process results
    header = [
        'vtx_id', f'{prefix}_hit_id', f'{prefix}_description',
        f'{prefix}_score', f'{prefix}_query_coverage',
        f'{prefix}_align_length', f'{prefix}_align_identity',
        f'{prefix}_gaps', f'{prefix}_evalue'
    ]
    
    df = pd.read_csv(
        output_prefix.joinpath(f'{input_file.stem}.{program}.out'),
        sep='\t',
        names=header
    )
    df = df.sort_values(
        by=f'{prefix}_score', ascending=False
    ).groupby('vtx_id').first()

    if output_name:
        df.to_parquet(output_prefix.joinpath(output_name))
    
    return df
