import pathlib
module_path = __file__.replace('dashboard/etl/__init__.py', '')
CACHE_DIR = pathlib.Path(module_path) / 'cache_20240612_v_1_13'
DATA_DIR = pathlib.Path(module_path) / 'data_20240610_v_1_12'
NOTEBOOK_DATA_DIR = pathlib.Path(module_path) / 'scripts' / 'notebooks'
HMMER_S3_LOC = 's3://velia-data-dev/VDC_004_annotation/hmmer/latest/'
TPM_DESEQ2_FACTOR = 80
DB_CONNECTION_STRING = 'postgresql://postgres:veliadata@vla-iac-dev01-rds-postgres-dbinstance-etx7xhyugv9t.csszhbmgzzwc.us-west-2.rds.amazonaws.com/expression_atlas_db'
REDSHIFT_CONNECTION_STRING = 'redshift+psycopg2://VeliaData1:VeliaData1@expressionatlasdev.328315166908.us-west-2.redshift-serverless.amazonaws.com:5439/expression_atlas_db'

