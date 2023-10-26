"""
Dashboard has two components of input information.

1. data
    a. contains information that either doesn't change, can't be automated, or updates 
        on a different cadence than cache files below
        i. ribo-seq results
        ii. expression atlas
        iii. synonym mappings etc and screening data 
2. cache
    a. contains automatically generated files that update on a cadence according to 
        new VTX generation. (e.g. everytime we screen a new collection)
    
"""
from dashboard.etl import transcript_features
import pandas as pd
import sqlalchemy
import shutil

def update_autoimmune_expression_atlas(transcript_subset = None):
    session, transcript_atlas_sample_expression = transcript_features.create_de_database('sqlite:///example.db')
    transcript_atlas_sample_expression.to_parquet('../../data/autoimmune_de_samples.parq')
    shutil.move('example.db', '../../data/autoimmune_de.db')
    
if __name__ == '__main__':
    update_autoimmune_expression_atlas()