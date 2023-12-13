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
import boto3


def download_file_from_s3(bucket_name, s3_object_key, local_file_path):
    # Create an S3 client
    s3 = boto3.client('s3')

    # Download the file
    try:
        s3.download_file(bucket_name, s3_object_key, local_file_path)
        print(f"File downloaded successfully to {local_file_path}")
    except Exception as e:
        print(f"Error downloading the file: {e}")


def update_autoimmune_expression_atlas(transcript_subset = None):
    db_path = '/home/ec2-user/dashboard/repos/data/autoimmune_de.db'
    # Specify your bucket name and object key
    bucket_name = "velia-piperuns-dev"
    object_key = "expression_atlas/v1/autoimmune_expression_atlas_v1.db"
    local_file_path = "/home/ec2-user/repos/dashboard/data/autoimmune_expression_atlas_v1.db"
    download_file_from_s3(bucket_name, object_key, local_file_path)

    
if __name__ == '__main__':
    update_autoimmune_expression_atlas()