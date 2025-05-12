"""
Note, unlike other scripts, this script currently has to be manually changed to process different folders.
"""

import os
import re
import copy
import requests
import shutil
import json
from multiprocessing import Pool

TARGET_METRIC_THRESHOLD = 0.4
TARGET_METRIC = 'ipTM'
RESULTS_FOLDER = 'output_raw_results/GEF_1500AA_2000AA_AF3_SCREEN'
RESULTS_CSV = RESULTS_FOLDER + '_results.csv'
# TODO: need to implement for AF3
# CREATE_REFINED_RUN_SCRIPT = False
# REFINED_RUN_FOLDER = RESULTS_FOLDER + '_REFINED'
# REFINED_RUN_SCRIPT = REFINED_RUN_FOLDER + '_RUN.sh'

def main():
    results = []

    with Pool(os.cpu_count() // 2) as pool:
        results = pool.map(process_model_directory, [os.path.join(RESULTS_FOLDER, file) for file in os.listdir(RESULTS_FOLDER) if os.path.isdir(os.path.join(RESULTS_FOLDER, file))])

    results_sorted = copy.deepcopy(results)
    results_sorted.sort(key=lambda result: max(prediction[TARGET_METRIC] for prediction in result[1]), reverse=True)

    candidate_results = [result for result in results_sorted if max(prediction[TARGET_METRIC] for prediction in result[1]) > TARGET_METRIC_THRESHOLD]
    candidate_folders = [result[0] for result in candidate_results]

    print("\n\n\n\n=================== ALL RESULTS ===================")
    for result in results_sorted:
        # print(result)
        uniprot_id = result[0].split('/')[1].split('_')[0]
        # uniprot edge case
        if len(uniprot_id) == 2:
            uniprot_id = result[0].split('/')[1].split('_')[1]
        print(f'UNIPROT LINK: https://www.uniprot.org/uniprotkb/{uniprot_id}/entry')

    with open(RESULTS_CSV, 'w') as f:
        f.write('model,fasta,seed,sample,uniprot link,ipTM,pLDDT,pTM\n')
        for result in results_sorted:
            uniprot_id = result[0].split('/')[1].split('_')[0]
            # uniprot edge case
            if len(uniprot_id) == 2:
                uniprot_id = result[0].split('/')[1].split('_')[1]
            for prediction in result[1]:
                f.write(f'{prediction["model"]},{result[0]}.fasta,{prediction["seed"]},{prediction["sample"]},https://www.uniprot.org/uniprotkb/{uniprot_id}/entry,{prediction["ipTM"]},{prediction["pLDDT"]},{prediction["pTM"]}\n')

            # best_prediction = max(result[1], key=lambda prediction: prediction[TARGET_METRIC])
            # print(best_prediction)
            # f.write(f'{result[0]}.fasta,https://www.uniprot.org/uniprotkb/{uniprot_id}/entry,{best_prediction["ipTM"]},{best_prediction["pLDDT"]},{best_prediction["pTM"]}\n')

    print(f'\n\n\n\n=================== CANDIDATE RESULTS ({TARGET_METRIC} > {TARGET_METRIC_THRESHOLD}) ===================')
    for result in candidate_results:
        print(f'{TARGET_METRIC}: {max(prediction[TARGET_METRIC] for prediction in result[1])}')
        uniprot_id = result[0].split('/')[1].split('_')[0]
        # uniprot edge case
        if len(uniprot_id) == 2:
            uniprot_id = result[0].split('/')[1].split('_')[1]
        print_uniprot_details(uniprot_id)
        print(f'Full result details: {result}\n')

def process_model_directory(model_path):
    print(f"Processing {model_path}...", flush=True)
    model_scores = []
    
    for seed in os.listdir(model_path):
        if not seed.startswith('seed-'):
            continue
        seed_path = os.path.join(model_path, seed)
        if not os.path.isdir(seed_path):
            continue

        full_seed_match = re.search(r"seed-(\d+)_sample-(\d+)", seed) 
        seed = full_seed_match.group(1)
        sample = full_seed_match.group(2)

        confidence_file = [f for f in os.listdir(seed_path) if f.endswith('confidences.json')]
        if len(confidence_file) == 0:
            print(f"Warning: No confidences.json found in {seed_path}. Skipping this seed.")
            continue

        confidence_json_data = json.load(open(os.path.join(seed_path, confidence_file[0])))
        plddt = sum(confidence_json_data['atom_plddts']) / len(confidence_json_data['atom_plddts']) / 100.0

        summary_confidences_file = [f for f in os.listdir(seed_path) if f.endswith('summary_confidences.json')]
        if len(summary_confidences_file) == 0:
            print(f"Warning: No summary_confidences.json found in {seed_path}. Skipping this seed.")
            continue
        
        summary_confidences_json_data = json.load(open(os.path.join(seed_path, summary_confidences_file[0])))

        seed_model_path = os.path.join(seed_path, [f for f in os.listdir(seed_path) if f.endswith('.cif')][0])

        model_scores.append({
            'model': seed_model_path,
            'seed': seed,
            'sample': sample,
            'pLDDT': plddt,
            'pTM': summary_confidences_json_data['ptm'],
            'ipTM': summary_confidences_json_data['iptm'],
        })

    return model_path, model_scores

def print_uniprot_details(uniprot_id):
    uniprot_details = get_uniprot_details(uniprot_id)
    if uniprot_details is not None:
        for key, value in uniprot_details.items():
            print(f'{key}: {value}')
        print(f'UNIPROT LINK: https://www.uniprot.org/uniprotkb/{uniprot_id}/entry')
    
    return uniprot_details

def get_uniprot_details(uniprot_id):
    """
    Fetch details about a UniProt ID using the UniProt API.
    Args:
        uniprot_id (str): The UniProt ID to fetch.
    Returns:
        dict: The filtered details of the UniProt entry.
    """
    base_url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}"
    
    try:
        response = requests.get(base_url, headers={"Accept": "application/json"})
        response.raise_for_status()  # Raise an error for HTTP errors
        json = response.json()
        id = json['uniProtkbId']	
        try:
            full_name = json['proteinDescription']['recommendedName']['fullName']['value'] 
        except:
            full_name = ''
        if full_name == '':
            try:
                full_name = json['proteinDescription']['submissionNames'][0]['fullName']['value']
            except:
                full_name = ''
        organism = json['organism']['scientificName']
        

        return {'id': id, 'full_name': full_name, 'organism': organism}
    except requests.exceptions.RequestException as e:
        # print(f"An error occurred: {e}")
        return None


def archive_best_results(source_paths, destination_name):
    """
    Copy multiple folders to a new location.
    
    Args:
        source_paths (list): List of paths to folders that need to be backed up
        destination_name (str): Name of the destination folder.
    
    Returns:
        str: Path to the created tar.gz file
    """
    
    print("Creating archive of candidate results...")
    # Create the destination folder
    os.makedirs(destination_name, exist_ok=True)
    
    # Copy each source folder to the destination
    for source_path in source_paths:
        prediction_path = source_path + '_prediction'
        fasta_path = source_path + '.fasta'
        if not os.path.exists(prediction_path):
            raise FileNotFoundError(f"Path not found: {prediction_path}")
        
        # Get the base folder name
        folder_name = os.path.basename(os.path.normpath(prediction_path))
        destination_path = os.path.join(destination_name, folder_name)
        
        # Copy the folder and its contents
        shutil.copytree(prediction_path, destination_path, dirs_exist_ok=True)
        shutil.copy(fasta_path, destination_name)

if "__name__" == "__main__":
    main()