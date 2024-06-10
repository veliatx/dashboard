import os
import pathlib
import pandas as pd
import json
import argparse
import biolib
import re
import json
from Bio import SeqIO
import subprocess
import shlex


working_directory = os.path.abspath(os.curdir)

parser = argparse.ArgumentParser(description='')

# Add arguments
parser.add_argument('--input', '-i', type=str, help='Input file path')
parser.add_argument('--output_prefix', '-o', type=str, help='Output file path')
parser.add_argument('--verbose', '-v', action='store_true', help='Enable verbose mode')

def phobius_format_conversion(lines, aa_seq):
    results = ['' for i in range(len(aa_seq))]
    for l in lines:
        elements = [i for i in l.split(' ') if not i == '']
        if elements[1] == 'SIGNAL':
            for ix in range(int(elements[2])-1, int(elements[3])):
                if results[ix] != '':
                    print('WARNING')
                results[ix] = 'S'
        if elements[1] == 'DOMAIN':
            for ix in range(int(elements[2])-1, int(elements[3])):
                if elements[-1] == 'CYTOPLASMIC.' and elements[-2] != 'NON':
                    if results[ix] != '':
                        print('WARNING')
                    results[ix] = 'I'
                if elements[-1] == 'CYTOPLASMIC.' and elements[-2] == 'NON':
                    if results[ix] != '':
                        print('WARNING')
                    results[ix] = 'O'
        if elements[1] == 'TRANSMEM':
            for ix in range(int(elements[2])-1, int(elements[3])):
                if results[ix] != '':
                    print('WARNING')
                results[ix] = 'M'
    return ''.join(results)

def parse_phobius(data, id_to_sequence):
    groups = {}
    current_group = []

    lines = data.split('\n')
    lines = [i for i in lines if len(i)>0]

    for line in lines:
        line = line.strip()
        if line.startswith('ID') and len(current_group)>0:
            if current_group:
                # groups.append(current_group)
                groups[id_to_sequence[current_group[0].replace('ID   ', '')]] = current_group[1:-1]
                current_group = [line]
        elif line == '\\':
            continue
        else:
            current_group.append(line)

    if current_group:
        groups[id_to_sequence[current_group[0].replace('ID   ', '')]] = current_group[1:-1]
    
    return {k: {'string_rep': phobius_format_conversion(v, k)} for k, v in groups.items()}

def parse_deeptmhmm(input_file):
    # lines = input_file.get_data().decode().split('\n')
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
                         
def parse_signalp6(input_file, id_to_sequence):
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
                if i<seq_length:
                    string_representation[i] = 'S'
            string_representation = ''.join(string_representation)
        else:
            cut_site = -1
            string_representation = 'C' * seq_length  
        
        results[id_to_sequence[ix]] = {'score': r['Likelihood'][1], 'cut_site': cut_site, 'string_rep': string_representation}
    return results

def parse_signalp5(input_file, id_to_sequence):
    lines = pd.read_table(input_file, skiprows=2, names=['ID', 'Prediction', 'score', 'OTHER', 'cut_site'])
    lines.fillna(-1, inplace=True)
    lines['ID'] = [str(i) for i in lines['ID']]
    results = {}
    pattern = r"CS pos: (\d+)"
    for ix, row in lines.iterrows():
        seq_length = len(id_to_sequence[row['ID']])
        match = re.search(pattern, str(row['cut_site']))
        if match:
            cs = int(match.group(1))
            string_representation = ['O'] * seq_length
            for i in range(0, cs):
                if i<seq_length:
                    string_representation[i] = 'S'
            string_representation = ''.join(string_representation)
        else:
            cs = -1
            string_representation = 'C' * seq_length
        results[id_to_sequence[row['ID']]] =  {'score': row['score'], 'cut_site': cs, 'string_rep': string_representation}
    return results

def parse_signalp41(input_file, id_to_sequence):
    pattern_cut = r"Cleavage site between pos. (\d+)"
    pattern_dscore = r"D=([\d.]+)"
    pattern_name = r"Name=(\S+)"
    results = {}
    for line in input_file.split('\n'):
        if not line.startswith('Name='):
            continue
        id = re.search(pattern_name, line).group(1)
        seq_length = len(id_to_sequence[id])
        cut_match = re.search(pattern_cut, line)
        if cut_match:
            cut_site = int(cut_match.group(1))
            string_representation = ['O'] * seq_length
            for i in range(0, cut_site):
                if i<seq_length:
                    string_representation[i] = 'S'
            string_representation = ''.join(string_representation)
        else:
            cut_site = -1
            string_representation = 'C' * seq_length
        dmatch = float(re.search(pattern_dscore, line).group(1))
        
        results[id_to_sequence[id]] =  {'score': dmatch, 'cut_site': cut_site, 'string_rep': string_representation}
    return results

def parse_deepsig(input_file, id_to_sequence):
    with open(input_file, 'r') as fh:
        lines = json.load(fh)
    results = {}
    for l in lines:
        seq_length = len(id_to_sequence[l['accession']])
        features = l['features']
        if (len(features) == 0) or not any([i['type'] == 'SIGNAL' for i in features]):
            string_representation = 'C' * seq_length
            results[l['sequence']['sequence']] = {'score': 0, 'cut_site': -1, 'string_rep': string_representation}
        else:
            f = [i for i in features if i['type']=='SIGNAL'][0]
            string_representation = ['O'] * seq_length
            for i in range(0, int(f['end'])):
                if i<seq_length:
                    string_representation[i] = 'S'
            string_representation = ''.join(string_representation)
            results[l['sequence']['sequence']] = {'score': f['score'], 'cut_site': f['end'], 'string_rep': string_representation}
    return results

def run_protein_feature_tools(INPUT_FASTA_FILE, OUTPUT_PREFIX, number_threads=8):
    ## Parse the arguments

    os.makedirs(OUTPUT_PREFIX, exist_ok=True)
    input_filename = pathlib.Path(INPUT_FASTA_FILE).parts[-1]
    seqs = {k: str(v.seq) for k, v in SeqIO.to_dict(SeqIO.parse(INPUT_FASTA_FILE, 'fasta')).items()}
    
    
    ## Phobius
    phobius_exec = f"docker run -v {OUTPUT_PREFIX}:/data 328315166908.dkr.ecr.us-west-2.amazonaws.com/secretions_tools:latest phobius -long /data/{input_filename}"
    phobius_data = subprocess.check_output(shlex.split(phobius_exec)).decode()
    with open(f'{OUTPUT_PREFIX}/phobius.results.txt', 'w') as fopen:
        fopen.write(phobius_data)
    
    ## DeepTMHMM
    os.system(f"BIOLIB_DOCKER_RUNTIME=nvidia biolib run --local DTU/DeepTMHMM --fasta {INPUT_FASTA_FILE}")
    os.system(f"mv biolib_results {OUTPUT_PREFIX}/biolib_results")

    ## SignalP6 slow
    source_weight_dir = "/efs/models/signalp6"
    dest_weight_dir = "/home/jupyter-rob:/usr/local/lib/python3.8/dist-packages/signalp/model_weights"
    cmd = f"docker run --gpus all --rm -v /efs/models/signalp6:/usr/local/lib/python3.8/dist-packages/signalp/model_weights -v {OUTPUT_PREFIX}:/home/work streptomyces/signalp signalp6 --fastafile /home/work/{input_filename} --output_dir /home/work --mode slow --bsize 200 --format txt --organism eukarya"
    subprocess.run(shlex.split(cmd))
    subprocess.run(shlex.split(f'rm -f {OUTPUT_PREFIX}/output_*plot.txt'))
    
    ## SignalP 5.0b
    cmd = f"docker run --gpus all --rm -v {OUTPUT_PREFIX}:/home/work streptomyces/signalp signalp -fasta /home/work/{input_filename} -format short -org euk -plot png -prefix /home/work/results"
    subprocess.run(shlex.split(cmd))
    
    ## SignalP 4.1
    cmd = f"docker run --gpus all -v {OUTPUT_PREFIX}:/data 328315166908.dkr.ecr.us-west-2.amazonaws.com/secretions_tools:latest signalp -f long -t euk /data/{input_filename}"
    signalp41_data = subprocess.check_output(shlex.split(cmd)).decode()
    with open(f'{OUTPUT_PREFIX}/signalp41.results.txt', 'w') as fopen:
        fopen.write(signalp41_data)
    
    ## DeepSig
    cmd = f"docker run --gpus all --rm -v {OUTPUT_PREFIX}:/data bolognabiocomp/deepsig -f /data/{input_filename} -o /data/deepsig.results -k euk -m json"
    subprocess.run(shlex.split(cmd))
    

    ## Parse results
    with open(f'{OUTPUT_PREFIX}/phobius.results.txt', 'r') as fopen:
        phobius_data = ''.join(fopen.readlines())
    phobius_results = parse_phobius(phobius_data, seqs)
    deep_tmhmm_results = parse_deeptmhmm(f'{OUTPUT_PREFIX}/biolib_results/predicted_topologies.3line')
    signalp_6_results = parse_signalp6(f'{OUTPUT_PREFIX}/output.json', seqs)
    signalp_5_results = parse_signalp5(f'{OUTPUT_PREFIX}/results_summary.signalp5', seqs)
    with open(f'{OUTPUT_PREFIX}/signalp41.results.txt', 'r') as fopen:
        signalp41_data = ''.join(fopen.readlines())
    signalp_41_results = parse_signalp41(signalp41_data, seqs)
    deepsig_results = parse_deepsig(f'{OUTPUT_PREFIX}/deepsig.results', seqs)

    deepsig_results = pd.DataFrame(deepsig_results).T
    signalp_6_results = pd.DataFrame(signalp_6_results).T
    signalp_5_results = pd.DataFrame(signalp_5_results).T
    signalp_41_results = pd.DataFrame(signalp_41_results).T
    phobius_results = pd.DataFrame(phobius_results).T
    deep_tmhmm_results = pd.DataFrame(deep_tmhmm_results).T
    string_representations = pd.DataFrame()
    string_representations['Deepsig'] = deepsig_results['string_rep']
    string_representations['SignalP 6slow'] = signalp_6_results['string_rep']
    string_representations['SignalP 5b'] = signalp_5_results['string_rep']
    string_representations['SignalP 4.1'] = signalp_41_results['string_rep']
    string_representations['Phobius'] = phobius_results['string_rep']
    string_representations['DeepTMHMM'] = deep_tmhmm_results['string_rep']
    
    vtx_ids = []
    df_ordr = []
    for k, v in seqs.items():
        vtx_ids.append(k)
        df_ordr.append(v)
    string_representations = string_representations.loc[df_ordr].copy()
    string_representations.index = vtx_ids
    string_representations['Sequence'] = [seqs[i] for i in string_representations.index]

    string_representations.to_csv(f'{OUTPUT_PREFIX}/sequence_features_strings.csv')
    
    scores = pd.DataFrame()
    scores['Deepsig'] = deepsig_results['score']
    scores['SignalP 6slow'] = signalp_6_results['score']
    scores['SignalP 5b'] = signalp_5_results['score']
    scores['SignalP 4.1'] = signalp_41_results['score']
    scores['Phobius'] = -1
    
    vtx_ids = []
    df_ordr = []
    for k, v in seqs.items():
        vtx_ids.append(k)
        df_ordr.append(v)
    scores = scores.loc[df_ordr].copy()
    scores.index = vtx_ids
    scores['Sequence'] = [seqs[i] for i in scores.index]
    
    scores.to_csv(f'{OUTPUT_PREFIX}/sequence_features_scores.csv')
        
    output_fmt = '6 qaccver saccver stitle score qcovs length pident gaps evalue'
    subprocess.run(shlex.split(f"docker run  --rm -v {OUTPUT_PREFIX}:/data -v /efs/databases/blast:/db ncbi/blast tblastn -outfmt '{output_fmt}' -db /db/mouse.rna.fna -query /data/{input_filename} -num_threads 32 -out /data/tblastn.results.csv"))
    
    cmd = f'docker run --gpus all -it -v {os.path.dirname(os.path.abspath(__file__))}:/opt/openfold -v {OUTPUT_PREFIX}:/data 328315166908.dkr.ecr.us-west-2.amazonaws.com/esmfold:latest python /opt/openfold/run_batch_fasta.py /data/{input_filename} /data/esmfold.jsonlines'
    subprocess.run(shlex.split(cmd))


if __name__ == '__main__':
    args, unknown = parser.parse_known_args()

    INPUT_FASTA_FILE = args.input#os.path.abspath(os.path.join(OUTPUT_DIR, 'protein_data', 'protein_tools_input.fasta'))
    OUTPUT_PREFIX = args.output_prefix#os.path.abspath(os.path.join(OUTPUT_DIR, 'protein_data'))
    run_protein_feature_tools(INPUT_FASTA_FILE, OUTPUT_PREFIX, number_threads=8)