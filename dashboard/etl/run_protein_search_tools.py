from dashboard import plotting, description
from dashboard.etl.sorf_query import load_jsonlines_table

import gzip
import shlex

import pathlib
import subprocess
import pandas as pd
import click

@click.command()
@click.argument('input_vtx', type=click.Path(exists=True, resolve_path=True))
@click.argument('output_prefix', type=click.Path(resolve_path=True))
@click.option("--number_threads", help="Number of threads to use for parallel computing.", default=8, type=int)

# input_vtx = #pathlib.Path('/home/ec2-user/repos/dashboard/cache/protein_data/protein_tools_input.fasta')
# CACHE_DIR = pathlib.Path('/home/ec2-user/repos/dashboard/cache')
# output_prefix = pathlib.Path('/home/ec2-user/repos/dashboard/cache/protein_data')

def run_protein_search_tools(input_vtx, output_prefix, number_threads):
    input_vtx = pathlib.Path(input_vtx)
    output_prefix = pathlib.Path(output_prefix)
    
    REFSEQ_FASTA_FILE = output_prefix.joinpath('GRCh38_latest_protein.faa')
    ENSEMBL_FASTA_FILE = output_prefix.joinpath('swissprot_proteins.fa')
    SWISSPROT_FASTA_FILE = output_prefix.joinpath('ensembl_proteins.fa')
    VTX_WITHOUT_SIGNAL_SEQ = output_prefix.joinpath('vtx_aa_seq_signal_sequence_removed.fa')
    any_missing = []
    for f in [REFSEQ_FASTA_FILE, ENSEMBL_FASTA_FILE, SWISSPROT_FASTA_FILE, VTX_WITHOUT_SIGNAL_SEQ]:
        if not f.exists():
            any_missing.append(f)
    if len(any_missing) > 0:
        raise ValueError(f"Reference files are missing the following files don't exist. \
            please check that update_cache.py has completed, or create these manually in the output_prefix directory.: {any_missing}")
    # Swissprot Isoform Check
    fasta_file = SWISSPROT_FASTA_FILE
    query_db = f'/root/{input_vtx.name}'
    target_db = f'/root/{fasta_file.name}'
    output_file = f'/root/{input_vtx.stem}_{fasta_file.stem}_alignments.m8'
    options = '--format-output query,target,evalue,qstart,qend,qlen,qcov,gapopen,pident,fident,alnlen,raw,bits'

    base_cmd = f"docker run --rm -v {output_prefix}:/root soedinglab/mmseqs2 mmseqs easy-search" 
    full_cmd = f'{base_cmd} {options} {query_db} {target_db} {output_file} /root/tmp'
    subprocess.run(shlex.split(full_cmd))

    isoform_df = pd.read_csv(output_prefix.joinpath(f'{input_vtx.stem}_{fasta_file.stem}_alignments.m8'), sep='\t', names=options.split()[1].split(','))
    isoform_df = isoform_df[(isoform_df['pident'] >= 100)]
    isoform_df.drop_duplicates(inplace=True)
    swissprot_isoform_df = isoform_df.groupby('query').aggregate(list)
    swissprot_isoform_df.rename(columns={'target': 'swissprot_isoform'}, inplace=True)

    # Ensembl Isoforms
    fasta_file = ENSEMBL_FASTA_FILE
    query_db = f'/root/{input_vtx.name}'
    target_db = f'/root/{fasta_file.name}'
    output_file = f'/root/{input_vtx.stem}_{fasta_file.stem}_alignments.m8'
    options = '--format-output query,target,evalue,qstart,qend,qlen,qcov,gapopen,pident,fident,alnlen,raw,bits'

    base_cmd = f"docker run --rm -v {output_prefix}:/root soedinglab/mmseqs2 mmseqs easy-search" 
    full_cmd = f'{base_cmd} {options} {query_db} {target_db} {output_file} /root/tmp'
    subprocess.run(shlex.split(full_cmd))
    
    isoform_df = pd.read_csv(output_prefix.joinpath(f'{input_vtx.stem}_{fasta_file.stem}_alignments.m8'), sep='\t', names=options.split()[1].split(','))
    isoform_df = isoform_df[(isoform_df['pident'] >= 100)]
    isoform_df.drop_duplicates(inplace=True)
    ensembl_isoform_df = isoform_df.groupby('query').aggregate(list)
    ensembl_isoform_df.rename(columns={'target': 'ensembl_isoform'}, inplace=True)

    # Refseq Isoforms
    fasta_file = REFSEQ_FASTA_FILE
    query_db = f'/root/{input_vtx.name}'
    target_db = f'/root/{fasta_file.name}'
    output_file = f'/root/{input_vtx.stem}_{fasta_file.stem}_alignments.m8'
    options = '--format-output query,target,evalue,qstart,qend,qlen,qcov,gapopen,pident,fident,alnlen,raw,bits'

    base_cmd = f"docker run --rm -v {output_prefix}:/root soedinglab/mmseqs2 mmseqs easy-search" 
    full_cmd = f'{base_cmd} {options} {query_db} {target_db} {output_file} /root/tmp'
    subprocess.run(shlex.split(full_cmd))
    
    isoform_df = pd.read_csv(output_prefix.joinpath(f'{input_vtx.stem}_{fasta_file.stem}_alignments.m8'), sep='\t', names=options.split()[1].split(','))
    isoform_df = isoform_df[(isoform_df['pident'] >= 100)]
    isoform_df.drop_duplicates(inplace=True)

    refseq_isoform_df = isoform_df.groupby('query').aggregate(list)
    refseq_isoform_df.rename(columns={'target': 'refseq_isoform'}, inplace=True)

#     swissprot_isoform_df.to_csv(pathlib.Path(output_prefix) / 'swissprot_isoform.csv')
#     ensembl_isoform_df.to_csv(pathlib.Path(output_prefix) / 'ensembl_isoform.csv')
#     refseq_isoform_df.to_csv(pathlib.Path(output_prefix) / 'refseq_isoform.csv')
    
    # Blast Tools
    blast_db_path = pathlib.Path('/efs/databases/blast')
    blast_db = '-db mouse.protein.genbank.faa'
    output_fmt = '6 qaccver saccver stitle bitscore qcovs length pident gaps evalue'
    options = f'-outfmt "{output_fmt}" -num_threads {number_threads}'
    query = f'-query /blast/data/{input_vtx.name}'
    output = f'-out /blast/data/{input_vtx.stem}.blastp.out'
    base_cmd = f'docker run --rm -v {blast_db_path}:/blast/blastdb -v {output_prefix}:/blast/data ncbi/blast'
    full_cmd = f'{base_cmd} blastp {options} {blast_db} {query} {output}'
    subprocess.run(shlex.split(full_cmd))
    
    header = ['vtx_id', 'blastp_hit_id', 'blastp_description', 'blastp_score',
          'blastp_query_coverage', 'blastp_align_length', 'blastp_align_identity', 
          'blastp_gaps', 'blastp_evalue']

    blastp_df = pd.read_csv(output_prefix.joinpath(f'{input_vtx.stem}.blastp.out'), sep='\t', names=header)
    bdf = blastp_df.sort_values(by='blastp_score', ascending=False).groupby('vtx_id').first()
    bdf.to_parquet(pathlib.Path(output_prefix).joinpath('mouse_blastp.parq'))
    
    blast_db_path = pathlib.Path('/efs/databases/blast')
    blast_db = '-db mouse.rna.fna'
    output_fmt = '6 qaccver saccver stitle score qcovs length pident gaps evalue'
    options = f'-outfmt "{output_fmt}" -num_threads {number_threads}'
    query = f'-query /blast/data/{input_vtx.name}'
    output = f'-out /blast/data/{input_vtx.stem}.tblastn.out'
    
    base_cmd = f'docker run --rm -v {blast_db_path}:/blast/blastdb -v {output_prefix}:/blast/data ncbi/blast'
    full_cmd = f'{base_cmd} tblastn {options} {blast_db} {query} {output}'
    subprocess.run(shlex.split(full_cmd))
    
    header = ['vtx_id', 'tblastn_hit_id', 'tblastn_description', 'tblastn_score',
          'tblastn_query_coverage', 'tblastn_align_length', 'tblastn_align_identity', 
          'tblastn_gaps', 'tblastn_evalue']

    tblastn_df = pd.read_csv(output_prefix.joinpath(f'{input_vtx.stem}.tblastn.out'), sep='\t', names=header)
    tdf = tblastn_df.sort_values(by='tblastn_score', ascending=False).groupby('vtx_id').first()
    tdf.to_parquet(pathlib.Path(output_prefix).joinpath('mouse_tblastn.parq'))
    
    
    subprocess.run(shlex.split(f'docker run --rm -v {output_prefix}:/data -v /efs/databases/blast:/db ncbi/blast blastp -task blastp-fast -outfmt 15 -db /db/mouse.protein.faa -query /data/{input_vtx.name} -max_target_seqs 20 -num_threads 1 -out /data/blastp.results.json'))
    
    # Merge isoform results
    isoform_df = pd.DataFrame([swissprot_isoform_df['swissprot_isoform'], ensembl_isoform_df['ensembl_isoform'], refseq_isoform_df['refseq_isoform']]).T
    # isoform_df.replace(pd.NA, 'None', inplace=True)
    isoform_df.to_parquet(pathlib.Path(output_prefix).joinpath('isoforms_search.parq'))
    
    # Run blast against sequneces with signal peptide removed
    #Blast Tools
    blast_db_path = pathlib.Path('/efs/databases/blast')
    blast_db = '-db mouse.protein.genbank.faa'
    output_fmt = '6 qaccver saccver stitle bitscore qcovs length pident gaps evalue'
    options = f'-outfmt "{output_fmt}" -num_threads {number_threads}'
    query = f'-query /blast/data/{VTX_WITHOUT_SIGNAL_SEQ.name}'
    output = f'-out /blast/data/{VTX_WITHOUT_SIGNAL_SEQ.stem}.blastp.out'
    base_cmd = f'docker run --rm -v {blast_db_path}:/blast/blastdb -v {output_prefix}:/blast/data ncbi/blast'
    full_cmd = f'{base_cmd} blastp {options} {blast_db} {query} {output}'
    subprocess.run(shlex.split(full_cmd))
    
    header = [
        'vtx_id', 'nonsig_blastp_refseq_id', 'nonsig_blastp_hit_id', 'nonsig_blastp_description', 'nonsig_blastp_score',
        'nonsig_blastp_align_length', 'nonsig_blastp_align_identity', 'nonsig_blastp_gaps', 'nonsig_blastp_evalue']
    
    non_sig_blastp_df = pd.read_csv(output_prefix.joinpath(f'{VTX_WITHOUT_SIGNAL_SEQ.stem}.blastp.out'), sep='\t', names=header)
    
    header = [
        'vtx_id', 'nonsig_tblastn_refseq_id', 'nonsig_tblastn_hit_id', 'nonsig_tblastn_description', 'nonsig_tblastn_score',
        'nonsig_tblastn_align_length', 'nonsig_tblastn_align_identity', 'nonsig_tblastn_gaps', 'nonsig_tblastn_evalue']
    output = f'-out /blast/data/{VTX_WITHOUT_SIGNAL_SEQ.stem}.tblastn.out'
    blast_db = '-db mouse.rna.fna'
    base_cmd = f'docker run --rm -v {blast_db_path}:/blast/blastdb -v {output_prefix}:/blast/data ncbi/blast'
    full_cmd = f'{base_cmd} tblastn {options} {blast_db} {query} {output}'
    subprocess.run(shlex.split(full_cmd))
    tblastn_df = pd.read_csv(output_prefix.joinpath(f'{VTX_WITHOUT_SIGNAL_SEQ.stem}.tblastn.out'), sep='\t', names=header)
    tdf = tblastn_df.sort_values(by='nonsig_tblastn_score', ascending=False).groupby('vtx_id').first()
    bdf = non_sig_blastp_df.sort_values('nonsig_blastp_score', ascending=False).groupby('vtx_id').first()
    output_nonsig_df = bdf.merge(tdf, left_index=True, right_index=True, how='outer')
    output_nonsig_df.to_parquet(output_prefix.joinpath('nonsignal_seq_blast_tblastn.parq'))

