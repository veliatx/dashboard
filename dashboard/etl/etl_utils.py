import numpy as np
import pandas as pd
from dashboard.etl import DATA_DIR
import pathlib
from collections import defaultdict
import json

from orfdb import base
from orfdb.base import Protein, ProteinXref, Orf, Dataset
import smart_open

def merge_sorf_df_blast(sorf_df, blastp_table, protein_data_path):
    sorf_df = sorf_df.copy()
    sorf_df = sorf_df.merge(pd.DataFrame(blastp_table).T, left_index=True, right_index=True, how='left')
    protein_scores = pd.read_csv(protein_data_path.joinpath('sequence_features_scores.csv'), index_col=0)

    sorf_df = sorf_df.merge(protein_scores[['Deepsig', 'SignalP 6slow', 'SignalP 5b', 'SignalP 4.1']],
                    left_index=True, right_index=True, how='left')
    protein_strings = pd.read_csv(protein_data_path.joinpath('sequence_features_strings.csv'), index_col=0)
    protein_strings.drop(['DeepTMHMM_prediction', 'DeepTMHMM_length'], axis=1, inplace=True)
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

def parse_blastp_json(file_path):
    hits_per_query = defaultdict(list)
    sorf_table_data = {}
    with open(file_path, 'r') as fopen:
        blastp = json.load(fopen)
        if 'BlastOutput2' in blastp:
            blastp = blastp['BlastOutput2']
    for entry in blastp:
        entry = entry['report']['results']
        q = entry['search']['query_title']
        hits = entry['search']['hits']
        if len(hits) == 0:
            # print('No alignments found with mouse')
            pass
        else:
            for h in hits:
                ids = []
                for item in h['description']:
                    ids.append(item['accession'])
                alignment = h['hsps']
                alignment = alignment[0]
                align_str = '  \n'.join([h['description'][0]['title'], alignment['qseq'], alignment['midline'], alignment['hseq']])
                alignment['hit_ids'] = ';'.join(ids)
                alignment['alignment'] = align_str
                hits_per_query[q].append(alignment)
                if isinstance(alignment, dict):
                    best_hit = alignment
                else:
                    best_hit = pd.DataFrame(alignment).sort_values('score', ascending=False).iloc[0]
                best_hit_description = [h for h in hits if h['num'] == best_hit['num']][0]['description'][0]
                sorf_table_data[q] = {'blastp_score': best_hit['score'],
                 'blastp_query_coverage': best_hit['align_len']/len(best_hit['qseq']),
                 'blastp_align_length': best_hit['align_len'],
                 'blastp_gaps': best_hit['gaps'],
                 'blastp_align_identity': best_hit['identity']/best_hit['align_len'],
                'blastp_subject': best_hit_description['id'],
                'blastp_hit_description': best_hit_description['title']
                }
    return hits_per_query, sorf_table_data

def fasta_write_veliadb_protein_sequences(protein_data_path):
    session = base.Session()
    # Additional Protein Features
    swissprot_query = \
    session.query(Protein)\
           .join(ProteinXref, ProteinXref.protein_id == Protein.id)\
           .join(Dataset, Dataset.id == ProteinXref.xref_dataset_id)\
           .filter(Dataset.name == 'swissprot')\
           .distinct(ProteinXref.protein_id)

    fasta_file = protein_data_path.joinpath( 'swissprot_proteins.fa')

    with open(fasta_file, 'w') as outfile:
        for protein in swissprot_query.all():
            outfile.write(f'>{protein.uniprot_id}\n{protein.aa_seq}\n')

    fasta_file = protein_data_path.joinpath('ensembl_proteins.fa')

    with open(fasta_file, 'w') as outfile:
        for protein in session.query(Protein).filter(Protein.ensembl_protein_id.ilike('ENSP%')).all():
            outfile.write(f'>{protein.ensembl_protein_id}\n{protein.aa_seq}\n')
    session.close()
    
    with smart_open.smart_open("s3://velia-analyses-dev/VAP_20230327_phase1to6_secreted/data/GRCh38_latest_protein.faa") as fopen:
        with open(protein_data_path.joinpath('GRCh38_latest_protein.faa'), 'wb') as fwrite:
            for line in fopen.readlines():
                fwrite.write(line)
                
def write_nonsignal_aa_sequences(protein_data_path):
    feature_df = pd.read_csv(protein_data_path.joinpath('sequence_features_strings.csv'), index_col=0)
    signal_cols = ['SignalP 6slow', 'SignalP 4.1', 'SignalP 5b', 'Deepsig']
    nonsignal_seqs = []

    for i, row in feature_df.iterrows():
        nonsignal_aa = ''
        for col in signal_cols:
            if row[col].startswith('S'):
                signal_len = row[col].count('S')
                nonsignal_aa = row['Sequence'][signal_len:].strip('*')
                break
        nonsignal_seqs.append(nonsignal_aa)

    feature_df['nonsignal_seqs'] = nonsignal_seqs

    feature_df['DeepTMHMM_prediction'] = feature_df.apply(lambda x: 'M' in x['DeepTMHMM'] or 'B' in x['DeepTMHMM'], axis=1)
    feature_df['DeepTMHMM_length'] = feature_df.apply(lambda x: x['DeepTMHMM'].count('M'), axis=1)
    feature_df.to_csv(protein_data_path.joinpath('sequence_features_strings.csv'))
    with open(protein_data_path.joinpath('vtx_aa_seq_signal_sequence_removed.fa'), 'w') as fwrite:
        for vtx, seq in feature_df['nonsignal_seqs'].items():
            if len(seq)>0:
                fwrite.write(f">{vtx}\n{seq}\n")
