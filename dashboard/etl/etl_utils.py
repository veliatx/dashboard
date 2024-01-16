import pandas as pd
from dashboard.etl import CACHE_DIR, DATA_DIR
import pathlib

from veliadb import base
from veliadb.base import Protein, ProteinXref, Orf, Dataset
import smart_open

def merge_sorf_df_blast(sorf_df, blastp_table):
    sorf_df = sorf_df.copy()
    sorf_df = sorf_df.merge(pd.DataFrame(blastp_table).T, left_index=True, right_index=True, how='left')
    protein_scores = pd.read_csv(CACHE_DIR.joinpath('protein_data', 'sequence_features_scores.csv'), index_col=0)

    sorf_df = sorf_df.merge(protein_scores[['Deepsig', 'SignalP 6slow', 'SignalP 5b', 'SignalP 4.1']],
                    left_index=True, right_index=True, how='left')
    protein_strings = pd.read_csv(CACHE_DIR.joinpath('protein_data', 'sequence_features_strings.csv'), index_col=0)
    protein_cutsite = protein_strings.apply(lambda x: x.str.find('SO')+1).replace(0, -1).drop('Sequence', axis=1)
    sorf_df = sorf_df.merge(protein_cutsite,
                    left_index=True, right_index=True, how='left', suffixes=('_score', '_cut'))

    id_data = pd.read_csv('s3://velia-data-dev/VDC_001_screening_collections/all_phases/interim_phase1to7_non-sigp_20230723.csv')

    id_data.index = id_data['vtx_id']
    sorf_df = sorf_df.merge(id_data[['trans1',
            'trans2', 'trans3', 'sec1', 'sec2', 'sec3', 'translated_mean',
            'secreted_mean', 'translated', 'swissprot_isoform', 'ensembl_isoform',
            'refseq_isoform', 'phylocsf_58m_avg', 'phylocsf_58m_max', 'phylocsf_58m_min',
            'phylocsf_vals']], left_index=True, right_index=True, how='left')

    with open(DATA_DIR.joinpath('all_secreted_phase1to7.txt'), 'r') as f:
        secreted_ids = [i.strip() for i in f.readlines()]

    sorf_df.insert(int(np.where(sorf_df.columns=='translated')[0][0]), 'secreted',
                   [True if i in secreted_ids else False for i in sorf_df.index])
    sorf_df['transcripts_exact'] = [tuple(i) for i in sorf_df['transcripts_exact']]
    return sorf_df

def fasta_write_veliadb_protein_sequences():
    session = base.Session()
    # Additional Protein Features
    swissprot_query = \
    session.query(Protein)\
           .join(ProteinXref, ProteinXref.protein_id == Protein.id)\
           .join(Dataset, Dataset.id == ProteinXref.xref_dataset_id)\
           .filter(Dataset.name == 'swissprot')\
           .distinct(ProteinXref.protein_id)

    fasta_file = CACHE_DIR.joinpath('protein_data', 'swissprot_proteins.fa')

    with open(fasta_file, 'w') as outfile:
        for protein in swissprot_query.all():
            outfile.write(f'>{protein.uniprot_id}\n{protein.aa_seq}\n')

    fasta_file = CACHE_DIR.joinpath('protein_data', 'ensembl_proteins.fa')

    with open(fasta_file, 'w') as outfile:
        for protein in session.query(Protein).filter(Protein.ensembl_protein_id.ilike('ENSP%')).all():
            outfile.write(f'>{protein.ensembl_protein_id}\n{protein.aa_seq}\n')
    session.close()
    
    with smart_open.smart_open("s3://velia-analyses-dev/VAP_20230327_phase1to6_secreted/data/GRCh38_latest_protein.faa") as fopen:
        with open(CACHE_DIR.joinpath('protein_data', 'GRCh38_latest_protein.faa'), 'wb') as fwrite:
            for line in fopen.readlines():
                fwrite.write(line)