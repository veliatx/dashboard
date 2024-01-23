
import streamlit as st
import streamlit.components.v1 as components

import pathlib
import urllib
import pandas as pd
import subprocess
from typing import Union
import sys
import numpy as np
import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt
import plotly.express as px

from pydeseq2.dds import DeseqDataSet

APP_NAME = 'geo_explorer'
DATA_PATH = '/home/ubuntu/repos/expression_atlas/geo_explorer'
RESULTS_FOLDER = 'geo_qc'

def download_counts(GEO_ID:str) -> Union[pd.DataFrame,None]:
    """Downloads counts from GEO dataBeta.
    
    Args:
        GEO_ID (str) Geo Dataset ID.
        DATA_PATH (str) Location of project.
        RESULTS_FOLDER (str) Name of folder to nest inside the project. 
    Returns:
        counts_df (Union[pd.DataFrame,None]) A counts dataframe, else None if download error.
    """
    global DATA_PATH, RESULTS_FOLDER
    data_path = pathlib.Path(DATA_PATH, GEO_ID, RESULTS_FOLDER, 'data')
    data_path.mkdir(parents=True, exist_ok=True)

    download_path = data_path.joinpath(f'{GEO_ID}_raw_counts_GRCh38.p13_NCBI.tsv.gz')
    download_url = 'https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc='+ \
                        f'{GEO_ID}&format=file&file={GEO_ID}_raw_counts_GRCh38.p13_NCBI.tsv.gz'
    print(download_path)
    if not download_path.exists():
        print("File doesn't exist. downloading...")
        try:
            path, message = urllib.request.urlretrieve(
                    download_url,
                    download_path,
                )
        except urllib.error.URLError as e:
            print('Check URL: %s' % download_url)
            print(e)
            return e
        except urllib.error.HTTPError as e:
            print('%s experiment probably does not exist geo' % GEO_ID)
            print(e)
            return e
        except Exception as e:
            print(e)
    else:
        print('File exists')
    counts_df = pd.read_csv(download_path, sep='\t')
    counts_df.set_index('GeneID', inplace=True)
    print(counts_df.shape)
    return counts_df

def download_metadata(
                    GEO_ID:str, 
                    SRA_ID:str,
                    counts_columns_order:pd.Series,
                    edirect_root:str='/usr/local/bin/edirect'
                    ) -> Union[pd.DataFrame,None]:
    """Download metadata associated with SRA project from NCBI.
    
    Args:
        GEO_ID (str) Geo Dataset ID.
        SRA_ID (str) SRP Identifier.
        counts_columns_order (pd.Series) Ordering of columns or sample names in counts dataframe.
        DATA_PATH (str) Location of project.
        RESULTS_FOLDER (str) Name of folder to nest inside the project. 
    Returns:
        metadata_df (Union[pd.DataFrame,None]) A metadata dataframe, else None if download error.
    """
    global DATA_PATH, RESULTS_FOLDER
    metadata_path = pathlib.Path(DATA_PATH, GEO_ID, RESULTS_FOLDER, 'data', f"{GEO_ID}_metadata.tsv")
    
    if not metadata_path.exists():
        cmd = f'{edirect_root}/esearch -db sra -query {SRA_ID} | '+ \
                        f'{edirect_root}/efetch -format xml | '+ \
                        f'{edirect_root}/xtract '+ \
                        '-pattern SAMPLE -element IDENTIFIERS/EXTERNAL_ID SAMPLE_ATTRIBUTE/VALUE SAMPLE_TITLE'
        print(cmd)
        result = subprocess.run(
                        cmd,
                        shell=True, 
                        capture_output=True,
                    )
        if len(result.stderr.decode()) != 0:
            return 
        
        with open(metadata_path, 'w') as f_out:
            f_out.write(result.stdout.decode())
    
    if metadata_path.stat().st_size == 0:
        print('Object downloaded is 0 bytes.')
        return
    
    metadata_df = pd.read_csv(metadata_path, sep='\t', header=None)
    metadata_df.rename({0: 'sample_name', 1: 'geo_id'}, inplace=True, axis=1)
    metadata_df = metadata_df[
                [c for c in metadata_df.columns if metadata_df[c].nunique() > 1 or metadata_df[c].nunique() == metadata_df.shape[0]]
                ]
    metadata_df.set_index('geo_id', inplace=True)
    metadata_df.drop('sample_name', inplace=True, axis=1)
    metadata_df = metadata_df.loc[counts_columns_order]
    metadata_df.columns = metadata_df.columns.map(lambda x: f'X{x}')
    metadata_df = metadata_df.fillna(lambda x: 'NA' if isinstance(x,str) else 0.)
    if metadata_df.shape[1] < 1:
        print('Did not find relevant metdata for: %s' % SRP_ID)
        return 
    return metadata_df

@st.cache_data()
def build_adata_process_counts(
                            SRP_ID:str,
                            GEO_ID:str,
                            counts_df:pd.DataFrame, 
                            metadata_df:pd.DataFrame, 
                            GENE_SUM_FILTER:int=10,
                            highly_variable_ngenes:int=1500) -> None:
    """Create adata, process and filter counts, plot PCA, etc.
    
    Args:
        SRP_ID (str)
        GEO_ID (str)
        counts_df (pd.DataFrame)
        metadata_df (pd.DataFrame)
        GENE_SUM_FILTER (int)
        highly_variable_ngenes (int)
    """
    
    adata = ad.AnnData(
                X=counts_df.T.copy(),
                obs=metadata_df.copy(),
                var=pd.DataFrame(None, index=counts_df.index.copy()),
                )
        
    # Filter adata on number counts per row. 
        
    adata = adata[:,
                ((adata.X > GENE_SUM_FILTER).sum(axis=0) > max(adata.shape[0] / 5, 3)) & 
                (adata.X.mean(axis=0) > GENE_SUM_FILTER)
                ]
    
    # Build dds for running vst calculations. 
        
    dds = DeseqDataSet(adata=adata.copy(), design_factors=adata.obs.columns.tolist())
    dds.vst()
        
    # Calculate highly variable genes. 
        
    dds.layers['log1p'] = sc.pp.log1p(dds.X.copy())
    dds.X = dds.layers['vst_counts'].copy()
    sc.pp.highly_variable_genes(dds, n_top_genes=highly_variable_ngenes, layer='log1p')
                                
    sc.pp.scale(dds)
    sc.pp.pca(dds, use_highly_variable=True)
    return adata, dds

def plot_pca(dds, pca_x, pca_y, color_column = 'X2'):
    # Step 2: Extract PCA coordinates
    pca_coordinates = dds.obsm["X_pca"][:, np.array([pca_x-1, pca_y-1])]  # Assuming you want the first two PCs
    pca_df = pd.DataFrame(pca_coordinates, columns=[f'PC{pca_x}', f'PC{pca_y}'], index=dds.obs_names)
    pca_df = pca_df.merge(dds.obs, left_index=True, right_index=True)

    # Step 3: Plot using Plotly Express
    fig = px.scatter(pca_df, x=f'PC{pca_x}', y=f'PC{pca_y}', color=color_column)  # Color by cluster if added
    return fig

def main():
    st.title("GEO Data Explorer")
    with st.container() as study_container:
        geo = st.text_input("Enter GEO Accession", value="")
        sra = st.text_input("Enter SRA Accession", value="")
    with st.container() as pca_container:
        if (geo == "") or (sra == ""):
            pass
        else:
            counts_df = download_counts(geo)
            metadata_df = download_metadata(geo, sra, counts_df.columns)
            adata, dds = build_adata_process_counts(sra, geo, counts_df, metadata_df)
            pca_items = [i for i in range(1, 45)]
            pca_x_selected = st.selectbox("PCA X", options=pca_items, index=0)
            remaining_pca_options = [i for i in pca_items if i != pca_x_selected]
            pca_y_selected = st.selectbox("PCA Y", options=remaining_pca_options, index=0)
            metadata_column = st.selectbox("Metadata Column", options = metadata_df.columns, index=0)
            fig = plot_pca(dds, pca_x_selected, pca_y_selected, color_column = metadata_column)
            st.plotly_chart(fig)

        
if __name__ == "__main__":
    
    main()