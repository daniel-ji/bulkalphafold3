# Running a PPI Complex Screen with BulkAlphaFold3 (One Bait Protein and Multiple Prey Proteins)

This guide will walk you through the steps to run a PPI complex screen using BulkAlphaFold3, which is designed to analyze AlphaFold3 predictions between a single bait protein and multiple prey proteins. For this guide, we will use the LRRK2 ROC-COR protein as the bait and a set of GTPase-Activating Proteins (GAPs) and Guanine Nucleotide Exchange Factors (GEFs) as prey proteins.

Follow the [environment setup guide](setting_up_environment.md) to set up AlphaFold3 on a LambdaLabs instance. 

1. **Create the bait FASTA file**
Create a FASTA file for the **single** bait protein (e.g., LRRK2 ROC-COR domains). It should be formatted as a regular FASTA file with the bait protein sequence, for example:

```
>LRRK2_ROC_COR
IIRFLQQRLKKAVPYNRMKLMIVGNT...
```

2. **Create the prey FASTA file**

Create a FASTA file(s) for the prey proteins (e.g., GAPs and GEFs). Each prey protein should be separated by a new entry in the FASTA file(s), like so:

```
>GAP1
MSTNPKPQRKTKQK...
>GAP2
MSTNPKPQRKTKQK...
>GEF1
MSTNPKPQRKTKQK...
>GEF2
MSTNPKPQRKTKQK...
```

Alternatively, the prey protein file can be a newline-separated list of Uniprot IDs, which will be automatically resolved to the corresponding FASTA sequences. For example:

```
O43182
Q5ZEZ3
Q92888
...
Q9Y2X3
```

3. **Providing a model file for clash calculations**

If you want to calculate clashes between the predicted prey proteins' binding location with a larger bait protein region (e.g., the RCKW region of LRRK2), you will need to provide a model file for the bait protein (PDB or mmCIF format). This model file should contain the full structure of the bait protein, as well as the additional region that will potentially clash with the prey proteins and can be used to filter out infeasible predictions. For example, the LRRK2 RCKW region can be provided to filter out prey protein predictions that clash with the kinase and WD40 domains of LRRK2. 

3. **Create the configuration file and modify constants.py**

After creating the bait and prey FASTA files, create a configuration file (e.g., `configs/pulldown/lrrk2_roc_cor.yml`) that defines the parameters for the screen. This file should have the following values:

#### Compute specific parameters
- `current_pipeline` (required): `pulldown` (this is the pipeline type for running a PPI complex screen with one bait protein and multiple prey proteins).
- `number_of_seeds` (default: `2`): The number of seeds to initialize for each prediction. Each seed will have 5 samples, so `5 * number_of_seeds` models will be generated for each prediction.
- `max_combined_seq_length` (default: `10000`): The maximum combined sequence length of the bait and prey proteins. This is used to limit the size of the predictions to fit within the GPU memory. Note that environment variables can be set to increase this value by using shared memory (see the [environment setup guide](setting_up_environment.md#predicting-larger-sequences)).
- `process_count` (default: `total_cpus // 2`): The number of CPU processes to run in parallel for downstream data processing. 

#### AlphaFold3 prediction parameters

- `bait_filename` (required): The path to the bait FASTA file.
- `superfolder_to_fasta_and_folder` (required): A dictionary mapping the superfolder names to a list of tuples, where each tuple contains the path to the prey FASTA file and the corresponding output folder name. This allows you to specify multiple prey protein files and their corresponding output folders for the predictions, and how they should be grouped together in the output. 

#### Downstream analysis parameters

- `plddt_sliding_window` (required): An integer representing the size of the sliding window to use for smoothing the pLDDT scores by residue. This is used to better filter out low-confidence predictions and/or disordered regions for downstream analysis and visualization (calculating binding domains, clashes, overall quality of the predictions, etc.). The value should be an odd integer, and a value of `-1` will disable the sliding window smoothing.
- `all_plddt_windows` (required): A list of integers representing the pLDDT windows to use for filtering the predictions. The values in this list will be used to smooth out the pLDDT scores by residue. Note that this is a bit of an ugly implementation, where every time the `plddt_sliding_window` is changed, the `all_plddt_windows` list should be updated to include the new value. 
- `prediction_threshold_metric` (default: `ipTM`): The metric to use for filtering the predictions. Can be `ipTM`, `pLDDT`, or `pTM`.
- `prediction_threshold_metric_value` (default: `0.4`): The value of the metric to use for filtering the predictions.
- `min_contacts_thresholds` (required): A list of integers representing the minimum number of contacts between the prey protein and the bait protein domains to consider a prediction as a potential binding domain.
- `domains_to_residues` (required): A dictionary mapping the domain names to a list of integers representing the residue ranges for each domain in the bait protein. This is used to calculate the binding domains of the prey proteins to the bait protein. The residue ranges should be 1-indexed and inclusive, and the values should be in the format `[start, end]`.
- `clashes_model` (optional): The path to the model file for the bait protein that will be used to calculate clashes with the prey proteins. This should be a PDB or mmCIF file that contains the full structure of the bait protein, including the additional region that will potentially clash with the prey proteins.
- `max_clashes_threshold` (required): The maximum number of clashes allowed between the prey protein and the bait protein for downstream visualization and analysis.

For example:

```yaml
current_pipeline: pulldown
number_of_seeds: 2
max_combined_seq_length: 3600
process_count: 6

bait_filename: input_fasta/bait_LRRK2_ROC_COR.fasta
superfolder_to_fasta_and_folder:
  GAP:
    - [input_fasta/uniprotkb_GTPase_Activating_Protein_AND_2024_11_09_0AA_600AA.fasta, GAP_0AA_600AA_AF3_SCREEN]
    - [input_fasta/uniprotkb_GTPase_Activating_Protein_AND_2024_11_09_600AA_1000AA.fasta, GAP_600AA_1000AA_AF3_SCREEN]
    - [input_fasta/uniprotkb_GTPase_Activating_Protein_AND_2024_11_09_1000AA_1500AA.fasta, GAP_1000AA_1500AA_AF3_SCREEN]
    - [input_fasta/uniprotkb_GTPase_Activating_Protein_AND_2024_11_09_1500AA_2000AA.fasta, GAP_1500AA_2000AA_AF3_SCREEN]
  GEF:
    - [input_fasta/uniprotkb_Guanine_Nucleotide_Exchange_F_2024_11_08_0AA_600AA.fasta, GEF_0AA_600AA_AF3_SCREEN]
    - [input_fasta/uniprotkb_Guanine_Nucleotide_Exchange_F_2024_11_08_600AA_1000AA.fasta, GEF_600AA_1000AA_AF3_SCREEN]
    - [input_fasta/uniprotkb_Guanine_Nucleotide_Exchange_F_2024_11_08_1000AA_1500AA.fasta, GEF_1000AA_1500AA_AF3_SCREEN]
    - [input_fasta/uniprotkb_Guanine_Nucleotide_Exchange_F_2024_11_08_1500AA_2000AA.fasta, GEF_1500AA_2000AA_AF3_SCREEN]

plddt_sliding_window: 11
all_plddt_windows: [-1, 11, 25]
prediction_threshold_metric: ipTM
prediction_threshold_metric_value: 0.4

min_contacts_thresholds: [1, 20, 80]
domains_to_residues:
  ROC:
    - 1 # 1317 - 1316
    - 205 # 1521 - 1316
  "COR-A":
    - 206 # 1522 - 1316
    - 355 # 1671 - 1316
  "COR-B":
    - 356 # 1672 - 1316
    - 557 # 1873 - 1316

clashes_model: "LRRK2_RCKW.pdb"
max_clashes_threshold: 1000

```

4. **bulkalphafold3.py**

Run `python bulkalphafold3.py` script to generate the AlphaFold3 one-liner submit script(s) to Docker to run AF3 predictions. This script will create bash scripts in the `run_scripts/` directory that can be executed to run the predictions on the LambdaLabs instance. 

5. **Running the AlphaFold3 predictions**

Run the generated bash script(s) in the `run_scripts/` directory to start the AlphaFold3 predictions (can take from several hours to a few days depending on the number of prey proteins and the GPU resources available).

6. **process_results.py**

Run `python process_results.py` to process the results from the AlphaFold3 predictions. This script will iterate through the results folder and create a CSV file in `output_raw_results/` with the results of the predictions (ipTM, pLDDT, pTM). It will also print out well-scoring models and some summary information.

7. **get_binding_domain.py**

Run `python get_binding_domain.py` to calculate the binding domains of the prey proteins to the bait protein. This script will take the modified mmCIF files and use ChimeraX to determine the number of contacts and buried area between the prey protein and various domains of the bait protein. It will create a new CSV file in `output_binding_domain/` with the results.

8. **calculate_all_clashes.py**

If you provided a clashes model in the configuration file, run `python calculate_all_clashes.py` to calculate the clashes between the prey proteins and the bait protein. This script will take the modified mmCIF files and use ChimeraX to determine the number of clashes between the prey protein and the provided larger bait protein region (e.g., the RCKW region of LRRK2). It will create a new CSV file in `output_all_clashes/` with the results.

9. **merge_results.py**

Run `python merge_results.py` to merge the results from the previous scripts (`process_results.py`, `get_binding_domain.py`, and `calculate_all_clashes.py`) into a single CSV file in `output_merged_results/`. This will create a comprehensive overview of the predictions, binding domains, and clashes.

10. **analyze_results.py**

Run `python analyze_results.py` to analyze the merged results and filter the predictions based on the prediction threshold metric (ipTM, pLDDT, or pTM) and the clashes score. This script will create a CSV & TSV file of good candidates (filtered by specified metrics and thresholds) and plots for visualization in the `output_analyze_results/` directory.
