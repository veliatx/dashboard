"""This program visulizes Ribo-seq data through streamlit."""

import os
import streamlit as st
import pandas as pd
import boto3
import json
from smart_open import open
import streamlit.components.v1 as components

session = boto3.Session()

def get_job_names_on_s3():
    """Retrieves top-level folder names from the velia-piperuns-dev bucket starting with VPR_orfcalling."""

    folder_names = []
    with open("s3://velia-piperuns-dev/summary/experiments_with_orfrater_results.csv", transport_params={'session': session}) as experiments_csv:
        for index, row in pd.read_csv(experiments_csv, index_col=0).iterrows():
            folder_names.append(row.iloc[0])
    return folder_names


def get_description(experiment):
    """Returns a short description for the given ribo-seq pipeline run."""

    with open(f"s3://velia-piperuns-dev/{experiment}/{experiment}.json", transport_params={'session': session}) as json_file:
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

    # load n_reads, n_covered_bases, sorf_length, coverage, and orfrater score
    riboseq_data = {}
    try:
        for line in open(f"s3://velia-piperuns-dev/{experiment_name}"
                         f"/output/orfrater_results/{experiment_name}_extended.final.coverage.bed", transport_params={'session': session}):
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
    except:
        return riboseq_data
            

    # read total
    total_reads = int(open(f"s3://velia-piperuns-dev/{experiment_name}/output/aligned/{experiment_name}_count.txt", transport_params={'session': session}).readline())

    # update and normalize(RPKM) n_reads
    for line in open(f"s3://velia-piperuns-dev/{experiment_name}"
                     f"/output/orfrater_results/{experiment_name}_all.final.coverage.bed", transport_params={'session': session}):
        elements = line.strip().split()
        candidate_id = elements[3]
        n_reads = int(elements[12])
        sorf_length = int(elements[14])
        rpkm = n_reads*(10**9)/(float(sorf_length*total_reads))
        riboseq_data[candidate_id][1] = rpkm

    return riboseq_data


def get_average_coverage():
    """Returns a dataframe table containing ribo-seq coverage avergaed across each cell-type."""

    with open("s3://velia-piperuns-dev/summary/avg_rpkm_all_experiments.csv", transport_params={'session': session}) as experiments_csv:
        coverage_df = pd.read_csv(experiments_csv, index_col=0)
    return coverage_df.round(0)


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
    column_names = ["RPKM", "orfrater score", "#covered bases", "sorf length", "%coverage",]
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
            axis=0, subset=[f"RPKM-{experiment_idx}"], vmin=0, vmax=500, cmap='Blues')

    return coverage_df

def page():
    st.title("sORF Ribo-seq Atlas")
    st.title("Average RPKM")
    st.dataframe(get_average_coverage())

    st.title("ORF-Rater scores")
    select_experiments(st.session_state)
    if st.session_state["experiments"]:
        for experiment_idx, experiment in enumerate(st.session_state["experiments"]):
            st.divider()
            st.write(f"Group {experiment_idx}: {shorten_experiment_names(experiment)}")
            st.write(f"\t{get_description(experiment)}")
        st.divider()
        st.dataframe(get_all_coverage(st.session_state["experiments"]))
