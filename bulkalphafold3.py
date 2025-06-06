import os
import re
import requests
import json

from files_helper import get_sequences_from_fasta
from constants import BAIT_FILENAME, SUPERFOLDER_TO_FASTA_AND_FOLDER, TEMPLATE_FILE, NUMBER_OF_SEEDS, MODEL_WEIGHTS_FOLDER, DATABASE_FOLDER

# TODO: Move to config/*.yaml file and constants.py
MAX_COMBINED_SEQ_LENGTH = 3600

def process_folder(input_fasta, output_folder):
    bash_script_file = output_folder + "_RUN.sh"
    bait_uniprot_ids = []  # list of uniprot ids
    bait_sequences = []  # list of tuples, (sequence id, sequence)
    prey_uniprot_ids = []  # list of uniprot ids
    prey_sequences = []  # list of tuples (sequence id, sequence)

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    bait_uniprot_ids, bait_sequences = get_sequences_from_fasta(BAIT_FILENAME)
    prey_uniprot_ids, prey_sequences = get_sequences_from_fasta(input_fasta)

    for uniprot_id in bait_uniprot_ids + prey_uniprot_ids:
        file = f"{output_folder}/{uniprot_id}.fasta"
        if not os.path.exists(file):
            url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
            response = requests.get(url)
            if response.status_code == 200:
                print("Downloaded " + url)
                with open(file, "wb") as f:
                    f.write(response.content)
            else:
                print("Failed to download: " + url)
        else:
            print("Already downloaded, skipping " + uniprot_id)

    # input bait protein file were uniprot ids, so set bait_sequences to the fetched uniprot sequences
    if len(bait_sequences) == 0:
        for bait_uniprot_id in bait_uniprot_ids:
            with open(f"{output_folder}/{bait_uniprot_id}.fasta", "r") as f:
                bait_lines = f.readlines()
            bait_sequences.append((bait_uniprot_id, "".join(bait_lines[1:])))

    # input prey protein file were uniprot ids, so set bait_sequences to the fetched uniprot sequences
    if len(prey_sequences) == 0:
        for prey_uniprot_id in prey_uniprot_ids:
            with open(f"{output_folder}/{prey_uniprot_id}.fasta", "r") as f:
                prey_lines = f.readlines()
            prey_sequences.append((prey_uniprot_id, "".join(prey_lines[1:])))

    template_json = json.load(open(TEMPLATE_FILE, "r"))
    template_json["modelSeeds"] = list(range(NUMBER_OF_SEEDS))
    assert template_json["sequences"][0]["protein"]["id"] == "PREY", "Template file should have PREY as the first sequence id"
    assert template_json["sequences"][1]["protein"]["id"] == "BAIT", "Template file should have BAIT as the second sequence id"

    # create combined input json files for alphafold3
    for prey_protein, prey_content in prey_sequences:
        for bait_protein, bait_content in bait_sequences:
            total_seq_length = -1
            af3_input_json = f"{output_folder}/{prey_protein}_{bait_protein}_input.json"
            with open(af3_input_json, "w") as af3_input_json_file:
                af3_input_json_content = template_json.copy()
                af3_input_json_content["name"] = f"{prey_protein}_{bait_protein}"
                af3_input_json_content["sequences"][0]["protein"]["sequence"] = prey_content.strip()
                af3_input_json_content["sequences"][1]["protein"]["sequence"] = bait_content.strip()
                af3_input_json_file.write(json.dumps(af3_input_json_content, indent=4))

            total_seq_length = len(prey_content) + len(bait_content)
            if total_seq_length > MAX_COMBINED_SEQ_LENGTH:
                print(f"Skipping {prey_protein}_{bait_protein} due to length {total_seq_length}")
                continue

    # create a bash script that runs all the alphafold3 jobs
    with open(bash_script_file, "w") as f:
        results_folder_abs = os.path.abspath(output_folder)
        output_command = f"""docker run -it --detach \
    --name af3_run_{output_folder} \
    --volume {results_folder_abs}:/root/af_input/{output_folder} \
    --volume {results_folder_abs}:/root/af_output/{output_folder} \
    --volume {MODEL_WEIGHTS_FOLDER}:/root/models \
    --volume {DATABASE_FOLDER}:/root/public_databases \
    --gpus all \
    alphafold3 \
    python run_alphafold.py \
    --input_dir=/root/af_input/{output_folder} \
    --model_dir=/root/models \
    --output_dir=/root/af_output/{output_folder}
"""
        f.write(output_command)


if __name__ == "__main__":
    for superfolder, folders in SUPERFOLDER_TO_FASTA_AND_FOLDER.items():
        for fasta, folder in folders:
            print(f"Processing {fasta} in {folder}", flush=True)
            process_folder(fasta, folder)