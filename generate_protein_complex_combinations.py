# TODO: Prevent writing to existing folder
# TODO: Support for fasta sequences with long IDs (that are possibily repetitive when truncated)
import json
import itertools
import os
import re

from files_helper import get_sequences_from_fasta
from constants import CURRENT_PIPELINE
assert CURRENT_PIPELINE == "complex", "This script is only for the complex pipeline"
from constants import CONFIG_FILE, CURRENT_PIPELINE, MAX_ID_LENGTH, MODEL_WEIGHTS_FOLDER, DATABASE_FOLDER, NUMBER_OF_SEEDS, TEMPLATE_FILE, INPUT_FASTA, FIXED_PROTEINS, OUTPUT_FOLDER, MAX_COMBINED_SEQ_LENGTH, CUSTOM_PREDICTIONS

DIGIT_TO_WORD = {
    '0': 'zero', '1': 'one', '2': 'two', '3': 'three', '4': 'four',
    '5': 'five', '6': 'six', '7': 'seven', '8': 'eight', '9': 'nine'
}

def clean_id(id):
    id = re.sub(r"\W", "_", id)
    # Replace all numbers with their digit-wise English words
    id = re.sub(r'\d+', lambda x: ''.join(DIGIT_TO_WORD[d] for d in x.group()), id)
    return id[:MAX_ID_LENGTH].upper()

def generate_protein_complex_combinations():
    """
    Generate all possible combinations of protein complexes, given a list of proteins and a set of fixed proteins (that must be included in each complex).

    Args:
    input_fasta (str): Path to the input FASTA file containing protein sequences. Split by '>' to get individual protein sequences.
    fixed_proteins_file (str): Path to a file containing a list of fixed proteins that must be included in each complex.
    """
    _, sequences = get_sequences_from_fasta(INPUT_FASTA)
    sequence_ids, sequence_fastas = zip(*sequences)

    fixed_proteins = [protein.strip()[0:MAX_ID_LENGTH] for protein in FIXED_PROTEINS if protein.strip()]

    if set(fixed_proteins) - set(sequence_ids):
        raise ValueError("Some fixed proteins are not present in the input FASTA file.")

    if len(set(fixed_proteins)) != len(fixed_proteins):
        raise ValueError(f"Fixed proteins file contains duplicate entries after truncation to {MAX_ID_LENGTH} characters. Please ensure all fixed proteins after truncation are unique.")
    
    unfixed_proteins = [p for p in sequence_ids if p not in fixed_proteins]

    combinations = []
    for r in range(0 if len(unfixed_proteins) > 0 else 1, len(unfixed_proteins) + 1):
        combinations.extend(itertools.combinations(unfixed_proteins, r))

    combinations = [tuple((tuple(fixed_proteins) + combo)) for combo in combinations]

    print(f"Generated {len(combinations)} combinations of protein complexes; from {len(unfixed_proteins)} unfixed proteins and {len(fixed_proteins)} fixed proteins.")

    custom_predictions = [tuple(custom_prediction) for custom_prediction in CUSTOM_PREDICTIONS] if CUSTOM_PREDICTIONS else []

    if custom_predictions:
        overlap = set(custom_predictions) & set(combinations)
        print(f"Adding {len(custom_predictions)} custom predictions to the combinations (found {len(overlap)} overlaps with generated combinations).")
        non_overlap_combinations = [combo for combo in custom_predictions if combo not in combinations]
        combinations.extend(non_overlap_combinations)
        print(f"{len(combinations)} total combinations after adding custom predictions.")

    if not os.path.exists(OUTPUT_FOLDER):
        os.makedirs(OUTPUT_FOLDER)

    template_json = json.load(open(TEMPLATE_FILE, "r"))
    template_json["modelSeeds"] = list(range(NUMBER_OF_SEEDS))

    # Create FASTA files and AF3 input JSON files for each combination
    for i, combo in enumerate(combinations):
        complex_name = "_".join(combo)

        total_seq_length = sum(len(sequence_fastas[sequence_ids.index(id)]) for id in combo)
        if total_seq_length > MAX_COMBINED_SEQ_LENGTH:
            print(f"Skipping combination {complex_name}: Total sequence length {total_seq_length} exceeds maximum allowed length {MAX_COMBINED_SEQ_LENGTH}.")
            continue
        else: 
            print(f"Processing combination {complex_name}: with total sequence length {total_seq_length}")

        fasta_content = "\n".join([f">{id}\n{sequence_fastas[sequence_ids.index(id)]}" for id in combo])
        with open(f"{OUTPUT_FOLDER}/{total_seq_length}_{complex_name}.fasta", "w") as f:
            f.write(fasta_content)

        af3_json_file = f"{OUTPUT_FOLDER}/{total_seq_length}_{complex_name}_input.json"
        af3_json = template_json.copy()
        af3_json["name"] = f"{complex_name}"
        af3_json["sequences"] = [{"protein": {"id": clean_id(id), "sequence": sequence_fastas[sequence_ids.index(id)]}} for id in combo]
        with open(af3_json_file, "w") as f:
            f.write(json.dumps(af3_json, indent=4))

    # Generate a Docker command to run AlphaFold3 for each complex
    output_bash_script = f"{OUTPUT_FOLDER}/run_alphafold3.sh"

    with open(output_bash_script, "w") as bash_file:
            # create a bash script that runs all the alphafold3 jobs
        results_folder_abs = os.path.abspath(OUTPUT_FOLDER)
        cleaned_output_folder_name = OUTPUT_FOLDER.replace("/", "_").replace(" ", "_")
        output_command = f"""docker run -it --detach \
    --name af3_run_{cleaned_output_folder_name} \
    --volume {results_folder_abs}:/root/af_input/{OUTPUT_FOLDER} \
    --volume {results_folder_abs}:/root/af_output/{OUTPUT_FOLDER} \
    --volume {MODEL_WEIGHTS_FOLDER}:/root/models \
    --volume {DATABASE_FOLDER}:/root/public_databases \
    --gpus all \
    alphafold3 \
    python run_alphafold.py \
    --input_dir=/root/af_input/{OUTPUT_FOLDER} \
    --model_dir=/root/models \
    --output_dir=/root/af_output/{OUTPUT_FOLDER}
"""
        bash_file.write(output_command)
    


# Example usage: python generate_protein_complex_combinations.py --input-fasta protein_complex_input/PxdA.fasta --fixed-proteins-file protein_complex_input/PxdA.fixed.txt --output-folder output_protein_complex/PxdA
if __name__ == "__main__":
    print(f"Using config file: {CONFIG_FILE}")
    if CURRENT_PIPELINE != "complex":
        raise ValueError("This script is only for the complex pipeline (individual protein complex prediction). Please ensure you are using the correct config file and have current_pipeline set to 'complex'.")
    generate_protein_complex_combinations()