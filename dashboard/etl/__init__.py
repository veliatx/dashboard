import pathlib
module_path = __file__.replace('dashboard/etl/__init__.py', '')
CACHE_DIR = pathlib.Path(module_path) / 'cache'
DATA_DIR = pathlib.Path(module_path) / 'data'
NOTEBOOK_DATA_DIR = pathlib.Path(module_path) / 'scripts' / 'notebooks'
TPM_DESEQ2_FACTOR = 80
DB_CONNECTION_STRING = 'postgresql://postgres:veliadata@ip-10-65-16-75.us-west-2.compute.internal/expression_atlas_db'
REDSHIFT_CONNECTION_STRING = 'redshift+psycopg2://VeliaData1:VeliaData1@redshiftdev01.328315166908.us-west-2.redshift-serverless.amazonaws.com:5439/redshiftdev01'

