"""
ETL script for processing protein sequences through various prediction tools.

This script runs multiple protein feature prediction tools including:
- Phobius for transmembrane topology and signal peptide prediction
- DeepTMHMM for transmembrane topology prediction  
- SignalP 6.0, 5.0b, and 4.1 for signal peptide prediction
- DeepSig for signal peptide prediction
- ESMFold for protein structure prediction
- BLAST for sequence similarity search

The results are parsed and saved in standardized formats.
"""

import os
import re
import json
import pathlib
import argparse
import subprocess
import shlex
import pandas as pd
from Bio import SeqIO


def phobius_format_conversion(lines: list, aa_seq: str) -> str:
    """Convert Phobius output to string representation.
    
    Args:
        lines: List of Phobius output lines
        aa_seq: Amino acid sequence
        
    Returns:
        String representation of predictions with:
        S = signal peptide
        I = cytoplasmic
        O = non-cytoplasmic  
        M = transmembrane
    """
    results = ['' for _ in range(len(aa_seq))]
    for line in lines:
        elements = [i for i in line.split(' ') if i]
        if elements[1] == 'SIGNAL':
            for ix in range(int(elements[2])-1, int(elements[3])):
                if results[ix] != '':
                    print('WARNING: Overlapping predictions')
                results[ix] = 'S'
        if elements[1] == 'DOMAIN':
            for ix in range(int(elements[2])-1, int(elements[3])):
                if elements[-1] == 'CYTOPLASMIC.' and elements[-2] != 'NON':
                    if results[ix] != '':
                        print('WARNING: Overlapping predictions')
                    results[ix] = 'I'
                if elements[-1] == 'CYTOPLASMIC.' and elements[-2] == 'NON':
                    if results[ix] != '':
                        print('WARNING: Overlapping predictions')
                    results[ix] = 'O'
        if elements[1] == 'TRANSMEM':
            for ix in range(int(elements[2])-1, int(elements[3])):
                if results[ix] != '':
                    print('WARNING: Overlapping predictions')
                results[ix] = 'M'
    return ''.join(results)


def parse_phobius(data: str, id_to_sequence: dict) -> dict:
    """Parse Phobius output into dictionary of predictions.
    
    Args:
        data: Raw Phobius output text
        id_to_sequence: Dictionary mapping sequence IDs to sequences
        
    Returns:
        Dictionary of predictions for each sequence
    """
    groups = {}
    current_group = []

    lines = [i for i in data.split('\n') if i]

    for line in lines:
        line = line.strip()
        if line.startswith('ID') and current_group:
            groups[id_to_sequence[current_group[0].replace('ID   ', '')]] = current_group[1:-1]
            current_group = [line]
        elif line == '\\':
            continue
        else:
            current_group.append(line)

    if current_group:
        groups[id_to_sequence[current_group[0].replace('ID   ', '')]] = current_group[1:-1]
    
    return {k: {'string_rep': phobius_format_conversion(v, k)} for k, v in groups.items()}


def parse_deeptmhmm(input_file: str) -> dict:
    """Parse DeepTMHMM output file.
    
    Args:
        input_file: Path to DeepTMHMM output file
        
    Returns:
        Dictionary of predictions for each sequence
    """
    with open(input_file, 'r') as fh:
        lines = [i.strip() for i in fh.readlines()]
    entry_indexes = [ix for ix, l in enumerate(lines) if l.startswith('>')]
    results = {}
    for ix in entry_indexes:
        predicted_class = lines[ix].split('|')[-1].strip()
        seq = lines[ix+1]
        prediction = lines[ix+2]
        results[seq] = {'string_rep': prediction, 'predicted_class': predicted_class}
    return results


def parse_signalp6(input_file: str, id_to_sequence: dict) -> dict:
    """Parse SignalP 6.0 output file.
    
    Args:
        input_file: Path to SignalP 6.0 output file
        id_to_sequence: Dictionary mapping sequence IDs to sequences
        
    Returns:
        Dictionary of predictions for each sequence
    """
    pattern = r"Cleavage site between pos. (\d+)"
    with open(input_file, 'r') as fh:
        lines = json.load(fh)
    results = {}
    for ix, r in lines['SEQUENCES'].items():
        match = re.search(pattern, r['CS_pos'])
        seq_length = len(id_to_sequence[r['Name']])
        if match:
            cut_site = int(match.group(1))
            string_representation = ['O'] * seq_length
            for i in range(0, cut_site):
                if i < seq_length:
                    string_representation[i] = 'S'
            string_representation = ''.join(string_representation)
        else:
            cut_site = -1
            string_representation = 'C' * seq_length
        
        results[id_to_sequence[ix]] = {
            'score': r['Likelihood'][1],
            'cut_site': cut_site,
            'string_rep': string_representation
        }
    return results


def parse_signalp5(input_file: str, id_to_sequence: dict) -> dict:
    """Parse SignalP 5.0b output file.
    
    Args:
        input_file: Path to SignalP 5.0b output file
        id_to_sequence: Dictionary mapping sequence IDs to sequences
        
    Returns:
        Dictionary of predictions for each sequence
    """
    lines = pd.read_table(
        input_file,
        skiprows=2,
        names=['ID', 'Prediction', 'score', 'OTHER', 'cut_site']
    )
    lines.fillna(-1, inplace=True)
    lines['ID'] = [str(i) for i in lines['ID']]
    results = {}
    pattern = r"CS pos: (\d+)"
    for _, row in lines.iterrows():
        seq_length = len(id_to_sequence[row['ID']])
        match = re.search(pattern, str(row['cut_site']))
        if match:
            cs = int(match.group(1))
            string_representation = ['O'] * seq_length
            for i in range(0, cs):
                if i < seq_length:
                    string_representation[i] = 'S'
            string_representation = ''.join(string_representation)
        else:
            cs = -1
            string_representation = 'C' * seq_length
        results[id_to_sequence[row['ID']]] = {
            'score': row['score'],
            'cut_site': cs,
            'string_rep': string_representation
        }
    return results


def parse_signalp41(input_file: str, id_to_sequence: dict) -> dict:
    """Parse SignalP 4.1 output file.
    
    Args:
        input_file: Raw SignalP 4.1 output text
        id_to_sequence: Dictionary mapping sequence IDs to sequences
        
    Returns:
        Dictionary of predictions for each sequence
    """
    pattern_cut = r"Cleavage site between pos. (\d+)"
    pattern_dscore = r"D=([\d.]+)"
    pattern_name = r"Name=(\S+)"
    results = {}
    for line in input_file.split('\n'):
        if not line.startswith('Name='):
            continue
        seq_id = re.search(pattern_name, line).group(1)
        seq_length = len(id_to_sequence[seq_id])
        cut_match = re.search(pattern_cut, line)
        if cut_match:
            cut_site = int(cut_match.group(1))
            string_representation = ['O'] * seq_length
            for i in range(0, cut_site):
                if i < seq_length:
                    string_representation[i] = 'S'
            string_representation = ''.join(string_representation)
        else:
            cut_site = -1
            string_representation = 'C' * seq_length
        dscore = float(re.search(pattern_dscore, line).group(1))
        
        results[id_to_sequence[seq_id]] = {
            'score': dscore,
            'cut_site': cut_site,
            'string_rep': string_representation
        }
    return results


def parse_deepsig(input_file: str, id_to_sequence: dict) -> dict:
    """Parse DeepSig output file.
    
    Args:
        input_file: Path to DeepSig output file
        id_to_sequence: Dictionary mapping sequence IDs to sequences
        
    Returns:
        Dictionary of predictions for each sequence
    """
    with open(input_file, 'r') as fh:
        lines = json.load(fh)
    results = {}
    for line in lines:
        seq_length = len(id_to_sequence[line['accession']])
        features = line['features']
        if (not features) or not any(i['type'] == 'SIGNAL' for i in features):
            string_representation = 'C' * seq_length
            results[line['sequence']['sequence']] = {
                'score': 0,
                'cut_site': -1,
                'string_rep': string_representation
            }
        else:
            feature = [i for i in features if i['type'] == 'SIGNAL'][0]
            string_representation = ['O'] * seq_length
            for i in range(0, int(feature['end'])):
                if i < seq_length:
                    string_representation[i] = 'S'
            string_representation = ''.join(string_representation)
            results[line['sequence']['sequence']] = {
                'score': feature['score'],
                'cut_site': feature['end'],
                'string_rep': string_representation
            }
    return results


def run_protein_feature_tools(input_fasta_file: str, output_prefix: str, number_threads: int = 8):
    """Run protein feature prediction tools and parse results.
    
    Args:
        input_fasta_file: Path to input FASTA file
        output_prefix: Path prefix for output files
        number_threads: Number of threads to use for parallel processing
    """
    os.makedirs(output_prefix, exist_ok=True)
    input_filename = pathlib.Path(input_fasta_file).parts[-1]
    seqs = {k: str(v.seq) for k, v in SeqIO.to_dict(SeqIO.parse(input_fasta_file, 'fasta')).items()}
    
    # Run Phobius
    phobius_exec = (
        f"docker run -v {output_prefix}:/data "
        "328315166908.dkr.ecr.us-west-2.amazonaws.com/secretions_tools:latest "
        f"phobius -long /data/{input_filename}"
    )
    phobius_data = subprocess.check_output(shlex.split(phobius_exec)).decode()
    with open(f'{output_prefix}/phobius.results.txt', 'w') as fopen:
        fopen.write(phobius_data)
    
    # Run DeepTMHMM
    os.system(f"BIOLIB_DOCKER_RUNTIME=nvidia biolib run --local DTU/DeepTMHMM --fasta {input_fasta_file}")
    os.system(f"mv biolib_results {output_prefix}/biolib_results")

    # Run SignalP 6.0 slow
    cmd = (
        "docker run --gpus all --rm "
        "-v /efs/models/signalp6:/usr/local/lib/python3.8/dist-packages/signalp/model_weights "
        f"-v {output_prefix}:/home/work streptomyces/signalp "
        f"signalp6 --fastafile /home/work/{input_filename} "
        "--output_dir /home/work --mode slow --bsize 200 --format txt --organism eukarya"
    )
    subprocess.run(shlex.split(cmd))
    subprocess.run(shlex.split(f'rm -f {output_prefix}/output_*plot.txt'))
    
    # Run SignalP 5.0b
    cmd = (
        f"docker run --gpus all --rm -v {output_prefix}:/home/work streptomyces/signalp "
        f"signalp -fasta /home/work/{input_filename} -format short -org euk -plot png "
        "-prefix /home/work/results"
    )
    subprocess.run(shlex.split(cmd))
    
    # Run SignalP 4.1
    cmd = (
        f"docker run --gpus all -v {output_prefix}:/data "
        "328315166908.dkr.ecr.us-west-2.amazonaws.com/secretions_tools:latest "
        f"signalp -f long -t euk /data/{input_filename}"
    )
    signalp41_data = subprocess.check_output(shlex.split(cmd)).decode()
    with open(f'{output_prefix}/signalp41.results.txt', 'w') as fopen:
        fopen.write(signalp41_data)
    
    # Run DeepSig
    cmd = (
        f"docker run --gpus all --rm -v {output_prefix}:/data bolognabiocomp/deepsig "
        f"-f /data/{input_filename} -o /data/deepsig.results -k euk -m json"
    )
    subprocess.run(shlex.split(cmd))

    # Parse results
    with open(f'{output_prefix}/phobius.results.txt', 'r') as fopen:
        phobius_data = ''.join(fopen.readlines())
    
    # Collect results from all tools
    phobius_results = parse_phobius(phobius_data, seqs)
    deep_tmhmm_results = parse_deeptmhmm(f'{output_prefix}/biolib_results/predicted_topologies.3line')
    signalp_6_results = parse_signalp6(f'{output_prefix}/output.json', seqs)
    signalp_5_results = parse_signalp5(f'{output_prefix}/results_summary.signalp5', seqs)
    with open(f'{output_prefix}/signalp41.results.txt', 'r') as fopen:
        signalp41_data = ''.join(fopen.readlines())
    signalp_41_results = parse_signalp41(signalp41_data, seqs)
    deepsig_results = parse_deepsig(f'{output_prefix}/deepsig.results', seqs)

    # Convert results to DataFrames
    deepsig_results = pd.DataFrame(deepsig_results).T
    signalp_6_results = pd.DataFrame(signalp_6_results).T
    signalp_5_results = pd.DataFrame(signalp_5_results).T
    signalp_41_results = pd.DataFrame(signalp_41_results).T
    phobius_results = pd.DataFrame(phobius_results).T
    deep_tmhmm_results = pd.DataFrame(deep_tmhmm_results).T

    # Combine string representations
    string_representations = pd.DataFrame()
    string_representations['Deepsig'] = deepsig_results['string_rep']
    string_representations['SignalP 6slow'] = signalp_6_results['string_rep']
    string_representations['SignalP 5b'] = signalp_5_results['string_rep']
    string_representations['SignalP 4.1'] = signalp_41_results['string_rep']
    string_representations['Phobius'] = phobius_results['string_rep']
    string_representations['DeepTMHMM'] = deep_tmhmm_results['string_rep']
    
    # Add sequence IDs
    vtx_ids = []
    df_ordr = []
    for k, v in seqs.items():
        vtx_ids.append(k)
        df_ordr.append(v)
    string_representations = string_representations.loc[df_ordr].copy()
    string_representations.index = vtx_ids
    string_representations['Sequence'] = [seqs[i] for i in string_representations.index]

    # Save string representations
    string_representations.to_csv(f'{output_prefix}/sequence_features_strings.csv')
    
    # Combine scores
    scores = pd.DataFrame()
    scores['Deepsig'] = deepsig_results['score']
    scores['SignalP 6slow'] = signalp_6_results['score']
    scores['SignalP 5b'] = signalp_5_results['score']
    scores['SignalP 4.1'] = signalp_41_results['score']
    scores['Phobius'] = -1
    
    # Add sequence IDs to scores
    scores = scores.loc[df_ordr].copy()
    scores.index = vtx_ids
    scores['Sequence'] = [seqs[i] for i in scores.index]
    
    # Save scores
    scores.to_csv(f'{output_prefix}/sequence_features_scores.csv')
        
    # Run BLAST search
    output_fmt = '6 qaccver saccver stitle score qcovs length pident gaps evalue'
    cmd = (
        f"docker run --rm -v {output_prefix}:/data -v /efs/databases/blast:/db "
        f"ncbi/blast tblastn -outfmt '{output_fmt}' -db /db/mouse.rna.fna "
        f"-query /data/{input_filename} -num_threads 32 -out /data/tblastn.results.csv"
    )
    subprocess.run(shlex.split(cmd))
    
    # Run ESMFold
    cmd = (
        f'docker run --gpus all -it '
        f'-v {os.path.dirname(os.path.abspath(__file__))}:/opt/openfold '
        f'-v {output_prefix}:/data '
        '328315166908.dkr.ecr.us-west-2.amazonaws.com/esmfold:latest '
        f'python /opt/openfold/run_batch_fasta.py '
        f'/data/{input_filename} /data/esmfold.jsonlines'
    )
    subprocess.run(shlex.split(cmd))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run protein feature prediction tools')
    parser.add_argument('--input', '-i', type=str, help='Input FASTA file path')
    parser.add_argument('--output_prefix', '-o', type=str, help='Output directory path')
    parser.add_argument('--verbose', '-v', action='store_true', help='Enable verbose mode')
    
    args, unknown = parser.parse_known_args()
    
    run_protein_feature_tools(args.input, args.output_prefix, number_threads=8)