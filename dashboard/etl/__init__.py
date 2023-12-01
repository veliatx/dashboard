import pathlib
module_path = __file__.replace('dashboard/etl/__init__.py', '')
CACHE_DIR = pathlib.Path(module_path) / 'cache'
DATA_DIR = pathlib.Path(module_path) / 'data'
TPM_DESEQ2_FACTOR = 80

PROTEIN_TOOLS_PATH = '/home/ubuntu/repos/protein_tools/dashboard_etl.py'