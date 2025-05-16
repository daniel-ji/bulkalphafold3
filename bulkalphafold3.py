import os
import re
import requests
import json

BAIT_FILENAME = "input_fasta/bait_LRRK2_ROC_COR.fasta"
PREY_FILENAME = "input_fasta/uniprotkb_Guanine_Nucleotide_Exchange_F_2024_11_08_1500AA_2000AA.fasta"
RESULTS_FOLDER = "GEF_1500AA_2000AA_AF3_SCREEN"
DATABASE_FOLDER = "$HOME/public_databases"
MODEL_WEIGHTS_FOLDER = "$HOME/"
BASH_SCRIPT_FILE = RESULTS_FOLDER + "_RUN.sh"
MAX_ID_LENGTH = 40
MAX_COMBINED_SEQ_LENGTH = 3600
TEMPLATE_FILE = "input_fasta/alphafold3_input_template.json"
NUMBER_OF_SEEDS = 2

bait_proteins = []  # list of uniprot ids
bait_sequences = []  # list of tuples, (sequence id, sequence)
prey_proteins = []  # list of uniprot ids
prey_sequences = []  # list of tuples (sequence id, sequence)

if not os.path.exists(RESULTS_FOLDER):
    os.makedirs(RESULTS_FOLDER)

with open(BAIT_FILENAME) as f:
    bait_proteins = [line.rstrip() for line in f if line.strip()]
    # provided bait file was actual sequences
    if bait_proteins[0][0] == ">":
        id = ""
        running_sequence = ""
        for line in bait_proteins:
            if line[0] == ">":
                if running_sequence != "":
                    assert id != ""
                    print("Found bait sequence " + id)
                    bait_sequences.append(id, running_sequence)
                id = re.sub(r"\W", "_", line[1 : MAX_ID_LENGTH + 1])
                running_sequence = ""
            else:
                running_sequence += line

        if running_sequence != "":
            assert id != ""
            print("Found bait sequence " + id)
            bait_sequences.append((id, running_sequence))
        bait_proteins = []

with open(PREY_FILENAME) as f:
    prey_proteins = [line.rstrip() for line in f if line.strip()]
    # provided prey file was actual sequences
    if prey_proteins[0][0] == ">":
        id = ""
        running_sequence = ""
        for line in prey_proteins:
            if line[0] == ">":
                if running_sequence != "":
                    assert id != ""
                    print("Found prey sequence " + id)
                    prey_sequences.append((id, running_sequence))
                id = re.sub(r"\W", "_", line[1 : MAX_ID_LENGTH + 1])
                running_sequence = ""
            else:
                running_sequence += line

        if running_sequence != "":
            assert id != ""
            print("Found prey sequence " + id)
            prey_sequences.append((id, running_sequence))
        prey_proteins = []

for protein in bait_proteins + prey_proteins:
    file = f"{RESULTS_FOLDER}/{protein}.fasta"
    if not os.path.exists(file):
        url = f"https://rest.uniprot.org/uniprotkb/{protein}.fasta"
        response = requests.get(url)
        if response.status_code == 200:
            print("Downloaded " + url)
            with open(file, "wb") as f:
                f.write(response.content)
        else:
            print("Failed to download: " + url)
    else:
        print("Already downloaded, skipping " + protein)

# input bait protein file were uniprot ids, so set bait_sequences to the fetched uniprot sequences
if len(bait_sequences) == 0:
    for bait_protein in bait_proteins:
        with open(f"{RESULTS_FOLDER}/{bait_protein}.fasta", "r") as f:
            bait_lines = f.readlines()
        bait_sequences.append((bait_protein, "".join(bait_lines[1:])))

# input prey protein file were uniprot ids, so set bait_sequences to the fetched uniprot sequences
if len(prey_sequences) == 0:
    for prey_protein in prey_proteins:
        with open(f"{RESULTS_FOLDER}/{prey_protein}.fasta", "r") as f:
            prey_lines = f.readlines()
        prey_sequences.append((prey_protein, "".join(prey_lines[1:])))

template_json = json.load(open(TEMPLATE_FILE, "r"))
template_json["modelSeeds"] = list(range(NUMBER_OF_SEEDS))
assert template_json["sequences"][0]["protein"]["id"] == "PREY", "Template file should have PREY as the first sequence id"
assert template_json["sequences"][1]["protein"]["id"] == "BAIT", "Template file should have BAIT as the second sequence id"

# create combined input json files for alphafold3
for prey_protein, prey_content in prey_sequences:
    for bait_protein, bait_content in bait_sequences:
        total_seq_length = -1
        af3_input_json = f"{RESULTS_FOLDER}/{prey_protein}_{bait_protein}_input.json"
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
with open(BASH_SCRIPT_FILE, "w") as f:
    results_folder_abs = os.path.abspath(RESULTS_FOLDER)
    output_command = f"""docker run -it --detach \
--name af3_run_{RESULTS_FOLDER} \
--volume {results_folder_abs}:/root/af_input/{RESULTS_FOLDER} \
--volume {results_folder_abs}:/root/af_output/{RESULTS_FOLDER} \
--volume {MODEL_WEIGHTS_FOLDER}:/root/models \
--volume {DATABASE_FOLDER}:/root/public_databases \
--gpus all \
alphafold3 \
python run_alphafold.py \
--input_dir=/root/af_input/{RESULTS_FOLDER} \
--model_dir=/root/models \
--output_dir=/root/af_output/{RESULTS_FOLDER}
"""
    f.write(output_command)
