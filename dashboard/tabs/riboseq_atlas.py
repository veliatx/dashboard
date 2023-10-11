"""This program visulizes Ribo-seq data through streamlit."""

import os
import streamlit as st
import pandas as pd
import boto3
import json
from smart_open import open
import streamlit.components.v1 as components


def get_job_names_on_s3():
    """Retrieves top-level folder names from the velia-piperuns-dev bucket starting with VPR_orfcalling."""

    folder_names = []
    client = boto3.client('s3')
    paginator = client.get_paginator('list_objects')
    result = paginator.paginate(Bucket='velia-piperuns-dev', Delimiter='/')
    for prefix in result.search('CommonPrefixes'):
        folder_name = prefix.get('Prefix')
        if folder_name.endswith("/"):
            folder_name = folder_name[:-1]

        if folder_name.startswith("VPR_orfcalling"):
            folder_names.append(folder_name)
    return folder_names


def get_description(experiment):
    """Returns a short description for the given ribo-seq pipeline run."""

    with open(f"s3://velia-piperuns-dev/{experiment}/{experiment}.json") as json_file:
        parameters = json.load(json_file)
    return parameters["note"]


def shorten_experiment_names(name):
    """Returns a reading-friendly experiment name."""

    return "-".join(name.split("_")[3:])


def select_experiments(state):
    """Adds user-input areas and collects inputs from users."""

    experiment_list = st.multiselect(
        'Experiments to display',
        st.session_state["job_names"],
        placeholder="Choose one or more experiments",
        format_func=shorten_experiment_names)
    
    if experiment_list != st.session_state["experiments"]:
        st.session_state["experiments"] = experiment_list


def get_coverge(experiment_name):
    """Returns coverage data for the given pipeline run."""

    bed_filename = os.path.join("./data", "%s_extended.final.coverage.bed" %experiment_name)
    riboseq_data = {}
    for line in open(f"s3://velia-piperuns-dev/{experiment_name}"
                     f"/output/orfrater_results/{experiment_name}_extended.final.coverage.bed"):
        elements = line.strip().split()
        candidate_id = elements[3]
        n_reads = int(elements[12])
        n_covered_bases = int(elements[13])
        sorf_length = int(elements[14])
        coverage = float(elements[15])
        orfrater = int(elements[16]) if elements[16] != "False" else -1
        
        if candidate_id in riboseq_data:
            print(line)
            print(riboseq_data[candidate_id][0])
        else:
            riboseq_data[candidate_id] = [
                line.strip(), n_reads, orfrater, n_covered_bases, sorf_length, coverage]

    # fix the inaccurate n_reads
    for line in open(f"s3://velia-piperuns-dev/{experiment_name}"
                     f"/output/orfrater_results/{experiment_name}_all.final.coverage.bed"):
        elements = line.strip().split()
        candidate_id = elements[3]
        n_reads = int(elements[12])
        orfrater = int(elements[16]) if elements[16] != "False" else -1
        riboseq_data[candidate_id][1] = n_reads
        if orfrater != riboseq_data[candidate_id][2]:
            exit("incorrect")
    return riboseq_data


def get_all_coverage(experiment_names):
    """Returns a dataframe table containing ribo-seq coverage infomation for the given experiments."""

    # load coverage data
    coverage = {}
    for experiment in experiment_names:
        for sorf_id, riboseq_data in get_coverge(experiment).items():
            if sorf_id not in coverage:
                coverage[sorf_id] = {}
            coverage[sorf_id][experiment] = riboseq_data

    # create dataframe
    # list df indices
    index_list = []
    for sorf_id in coverage:
        index_list.append(sorf_id)

    # list df columns
    column_list = []
    column_names = ["#mapped reads", "orfrater score", "#covered bases", "sorf length", "%coverage",]
    for experiment_idx, experiment in enumerate(experiment_names):
        for column_name in column_names:
            column_list.append(column_name+"-"+str(experiment_idx))

    # load data into a pandas dataframe
    coverage_df = []
    for sorf_id in index_list:
        current_sorf = []
        for experiment in experiment_names:
            for element_idx, element in enumerate(coverage[sorf_id][experiment]):
                if element_idx == 0:
                    continue
                else:
                    current_sorf.append(element)
        coverage_df.append(current_sorf)

    # color dataframe cells based on cell values
    coverage_df = pd.DataFrame(coverage_df, index=index_list, columns=column_list).style
    for experiment_idx in range(len(experiment_names)):
        coverage_df = coverage_df.background_gradient(
            axis=0, subset=[f"%coverage-{experiment_idx}"], vmin=0, vmax=1, cmap='YlOrRd')
        coverage_df = coverage_df.background_gradient(
            axis=0, subset=[f"#mapped reads-{experiment_idx}"], vmin=0, vmax=1000, cmap='Blues')

    return coverage_df

def page():
    st.title("sORF Ribo-seq Atlas")

    select_experiments(st.session_state)
    if st.session_state["experiments"]:
        for experiment_idx, experiment in enumerate(st.session_state["experiments"]):
            st.divider()
            st.write(f"Group {experiment_idx}: {shorten_experiment_names(experiment)}")
            st.write(f"\t{get_description(experiment)}")
        st.divider()
        st.dataframe(get_all_coverage(st.session_state["experiments"]))