import boto3
import smart_open
import pandas as pd
import numpy as np
from tqdm import tqdm

# def load_s3_transcript_DE_stats

def load_CCLE_data_from_s3():
    exp = pd.read_csv('s3://velia-analyses-dev/VAP_20230324_pull_ccle/data/OmicsExpressionProteinCodingGenesTPMLogp1_22Q4.csv')


def load_xena_transcripts_with_metadata_from_s3(transcripts_to_load = None):
    """
    Wrapper to load data from s3.

    Potentially this will be moved to redshift, but for now this handles loading Xena data.

    Loads as TPM + 0.001
    
    """
    xena_transcripts = pd.read_csv('s3://velia-data-dev/VDC_004_annotation/tcga/20230329_tcga_data/xena_transcript_names.csv')
    if transcripts_to_load is None:
        selected_transcripts = [str(t[0]) for t in xena_transcripts.values]
    else:
        overlapping_xena_sorf_transcripts = set(xena_transcripts['Xena Transcripts']).intersection(set(transcripts_to_load))
        selected_transcripts = [str(t) for t in overlapping_xena_sorf_transcripts]

    xena = pd.read_feather('s3://velia-data-dev/VDC_004_annotation/tcga/20230329_tcga_data/xena.feather', columns = ['index'] + selected_transcripts)
    xena.index = xena.pop('index')
    metadata = pd.read_table('s3://velia-data-dev/VDC_004_annotation/tcga/20230329_tcga_data/Xena/TcgaTargetGTEX_phenotype.txt',
                             encoding = "ISO-8859-1", index_col=0)
    metadata = metadata[~metadata['_primary_site'].isna()]
    ordr = xena.index.intersection(metadata.index)
    xena = xena.loc[ordr].copy()
    metadata = metadata.loc[ordr].copy()

    tissue_pairs = pd.read_excel('s3://velia-data-dev/VDC_004_annotation/tcga/20230329_tcga_data/Xena/GTEx_TCGA_comparisons.xlsx')
    tissue_pairs['GTEx Tissue Type'] = tissue_pairs['GTEx Tissue Type'].str.replace('whole ', '')
    tissue_pairs = tissue_pairs[tissue_pairs['GTEx Tissue Type']!='--']

    xena = np.exp2(xena)
    return metadata.merge(xena, left_index=True, right_index=True), metadata, tissue_pairs


def read_tcga_de_from_s3(bucket, object_prefix, output_dir = None,
                         tcga_cancer_codes = None):
    session = boto3.Session()
    s3 = session.resource('s3')
    my_bucket = s3.Bucket(bucket)
    
    cancer_de_results = {}
    for f in my_bucket.objects.filter(Prefix=object_prefix):
        fname = f.key
        if fname.endswith('.csv'):
            cancer_de_results[fname.split('/')[-1].split('_')[0]] = fname
    if tcga_cancer_codes is None:
        tcga_cancer_codes = cancer_de_results.keys()
    tables = {}
    for c in tqdm(tcga_cancer_codes):
        if c in cancer_de_results.keys():
            table = pd.read_csv(f"s3://{bucket}/{cancer_de_results[c]}", index_col=0)
            if not output_dir is None:
                table.to_parquet(f"{output_dir}/{c}_de.parq")
            tables[c] = table    
        else:
            print(f"{c} not found in results")
    return tables


def create_comparison_groups_xena_tcga_vs_normal(xena, tissue_pairs):
    groups = {}
    for ix, row in tissue_pairs.iterrows():
        # Normal sample indexes for specific tissues
        
        normal = xena.index[(xena['_primary_site'].str.lower() == row['GTEx Tissue Type'].lower())
                                & (xena['_sample_type'] == 'Normal Tissue')]
        # Cancer sample indexes
        if row['TCGA Cancer Type'] == 'LAML':
            cancer = xena.index[(xena['detailed_category'].str.lower() == row['Description'].lower())
                                    & (xena['_sample_type'] == 'Primary Blood Derived Cancer - Peripheral Blood')]
        else:
            cancer = xena.index[(xena['detailed_category'].str.lower() == row['Description'].lower())
                                    & (xena['_sample_type'] == 'Primary Tumor')]
        groups[row['TCGA Cancer Type']] = {'normal_indices': normal, 'cancer_indices': cancer,
                                       'GTEx Tissue': row['GTEx Tissue Type'], 'TCGA Cancer': row['Description']}
    return groups


def load_ccle_from_s3(path_to_ccle):
    pass


from sqlalchemy import create_engine, Column, Integer, String, Float
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker
import boto3

Base = declarative_base()



class TranscriptDE(Base):
    __tablename__ = 'transcript_differential_expression'
    id = Column(Integer, autoincrement = True, primary_key = True)
    transcript_id = Column(String, primary_key=False)
    reference_mean = Column(Float)
    group_mean = Column(Float)
    log2fc = Column(Float)
    log10padj = Column(Float)
    log10pval = Column(Float)
    contrast = Column(String)
    study = Column(String)
    
    
def load_autoimmune_de_tables_from_s3(bucket = 'velia-piperuns-dev', root_path = 'expression_atlas/v1/'):
    s3 = boto3.client('s3') 
    # List the objects in the specified bucket and prefix, without descending into subdirectories
    response = s3.list_objects_v2(
        Bucket=bucket,
        Prefix=root_path, 
        Delimiter = '/',
        MaxKeys = 21474837)

    # Pattern to match
    pattern = '*_transcript_*.csv'

    de_results = {}
    for p in tqdm(response['CommonPrefixes']):
        p = p['Prefix']
        files = s3.list_objects_v2(
                        Bucket='velia-piperuns-dev',
                        Prefix=f"{p}de_results/{p.split('/')[-2]}_transcript_", 
                        MaxKeys = 21474837)
        # if len(files) == 1:
        for f in files.get('Contents', []):
            pth = f"s3://velia-piperuns-dev/{f['Key']}"
            de_results[pth.split('/')[-1]] = pd.read_csv(pth, index_col=0)
        
    return de_results
    
def create_de_database(db_address):
    engine = create_engine(db_address)
    Base.metadata.create_all(engine)

    session = sessionmaker(bind=engine)()
    transcripts_to_add = []
    atlas_sample_expression_tables = []
    for n, df in load_autoimmune_de_tables_from_s3().items():
        contrast = n.replace('.csv', '').replace('_transcript_', '-')
        study = contrast.split('-')[0]
        contrast = contrast.replace(study+'-', '')
        for row in tqdm(df.itertuples()):
            t = TranscriptDE(transcript_id = row.Index, 
                        reference_mean = row._10, 
                        group_mean = row._9,
                        log2fc = row.log2FoldChange, 
                        log10padj = np.log10(row.padj), 
                        log10pval = np.log10(row.pvalue),
                        contrast = contrast,
                        study = study
                        )
            session.add(t)
        sample_value_column_names = df.columns[10:]
        sample_names = [i.split('_')[-1] for i in sample_value_column_names]
        group = ['CONTROL' if '_CONTROL' in i else 'DISEASE' for i in sample_value_column_names]
        multi_index = pd.MultiIndex.from_arrays([[study]*len(sample_value_column_names), sample_names, [contrast]*len(sample_value_column_names), group],
                                                names=('study', 'sample', 'contrast', 'group'))
        sample_expression = pd.DataFrame(df[sample_value_column_names].values, columns = multi_index, index = df.index)
        atlas_sample_expression_tables.append(sample_expression)
    atlas_sample_expression_tables = pd.concat(atlas_sample_expression_tables, axis=1, join='outer')
    atlas_sample_expression_tables.to_sql('sample_expression',
                                          con=engine,
                                          index=False,
                                          if_exists='replace')

    
    session.commit()
    return session, atlas_sample_expression_tables