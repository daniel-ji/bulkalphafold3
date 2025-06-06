import json
import itertools
import argparse
import os

from files_helper import get_sequences_from_fasta
from constants import MAX_ID_LENGTH, MODEL_WEIGHTS_FOLDER, DATABASE_FOLDER

NUMBER_OF_SEEDS = 5
TEMPLATE_FILE = "input_fasta/alphafold3_input_template_protein_complexes.json"

def generate_protein_complex_combinations(input_fasta, fixed_proteins_file, output_folder):
    """
    Generate all possible combinations of protein complexes, given a list of proteins and a set of fixed proteins (that must be included in each complex).

    Args:
    input_fasta (str): Path to the input FASTA file containing protein sequences. Split by '>' to get individual protein sequences.
    fixed_proteins_file (str): Path to a file containing a list of fixed proteins that must be included in each complex.
    """
    _, sequences = get_sequences_from_fasta(input_fasta)
    sequence_ids, sequence_fastas = zip(*sequences)

    if not os.path.exists(fixed_proteins_file):
        fixed_proteins_file = []

    with open(fixed_proteins_file, "r") as f:
        fixed_proteins = [line.strip()[0:MAX_ID_LENGTH] for line in f if line.strip()]

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

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    template_json = json.load(open(TEMPLATE_FILE, "r"))
    template_json["modelSeeds"] = list(range(NUMBER_OF_SEEDS))

    # Create FASTA files and AF3 input JSON files for each combination
    for i, combo in enumerate(combinations):
        complex_name = "_".join(combo)
        fasta_content = "\n".join([f">{id}\n{sequence_fastas[sequence_ids.index(id)]}" for id in combo])
        with open(f"{output_folder}/{complex_name}.fasta", "w") as f:
            f.write(fasta_content)

        af3_json_file = f"{output_folder}/{complex_name}_input.json"
        af3_json = template_json.copy()
        af3_json["name"] = f"{complex_name}"
        af3_json["sequences"] = [{"protein": {"id": id, "sequence": sequence_fastas[sequence_ids.index(id)]}} for id in combo]
        with open(af3_json_file, "w") as f:
            f.write(json.dumps(af3_json, indent=4))

    # Generate a Docker command to run AlphaFold3 for each complex
    output_bash_script = f"{output_folder}/run_alphafold3.sh"

    with open(output_bash_script, "w") as bash_file:
            # create a bash script that runs all the alphafold3 jobs
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
        bash_file.write(output_command)
    


# Example usage: python generate_protein_complex_combinations.py --input-fasta protein_complex_input/PxdA.fasta --fixed-proteins-file protein_complex_input/PxdA.fixed.txt --output-folder output_protein_complex/PxdA
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate all possible combinations of protein complexes from a FASTA file and a list of fixed proteins.")
    parser.add_argument("--input-fasta", type=str, help="Path to the input FASTA file containing protein sequences.", required=True)
    parser.add_argument("--output-folder", type=str, help="Folder where the generated protein complex combinations will be saved.", required=True)
    parser.add_argument("--fixed-proteins-file", type=str, help="Path to a file containing a list of fixed proteins that must be included in each complex.")

    args = parser.parse_args()

    generate_protein_complex_combinations(args.input_fasta, args.fixed_proteins_file, args.output_folder)