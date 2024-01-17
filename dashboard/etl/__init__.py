import pathlib
module_path = __file__.replace('dashboard/etl/__init__.py', '')
CACHE_DIR = pathlib.Path(module_path) / 'cache_merged'
DATA_DIR = pathlib.Path(module_path) / 'data'
NOTEBOOK_DATA_DIR = pathlib.Path(module_path) / 'scripts' / 'notebooks'
TPM_DESEQ2_FACTOR = 80

PROTEIN_TOOLS_PATH = '/home/ubuntu/repos/protein_tools/dashboard_etl.py'
