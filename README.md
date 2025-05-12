# BulkAlphaFold3

A simple pipeline to run AF3 PPI screens. Designed specifically for LRRK2 ROC-COR screens against GTPase-Activating Proteins (GAPs) and Guanine Nucleotide Exchange Factors (GEFs).

**Still in development.** This is a work in progress and not a complete solution (~80% complete). Please check back for updates. Documentation is also a work in progress.

## Requirements
- AlphaFold3
- Python 3.8+
- BioPython
- numpy
- pandas
- matplotlib
- seaborn

## Setting Up AlphaFold3 

See the [AlphaFold3 documentation](https://github.com/google-deepmind/alphafold3) for installation instructions. Install with Docker. Note that the AlphaFold3 installation is not included in this repository and you will need to obtain model weights and install it separately.

AF3 requires a powerful GPU to run. In our case (and possibly in yours), we did not have access to such a GPU, so we resorted to using LambdaLabs to run the predictions.

## Overall Workflow

## Script Run Order
1. `bulkaphafold3.py` - A script to generate an AlphaFold3 one-liner submit script to Docker to run AF3 predictions on the bait protein (LRRK2 ROC-COR) and the prey proteins (GAPs or GEFs). By default created in the `run_scripts` directory.
2. `process_results.py` - After the AF3 predictions run on the Docker container, this script will iterate through the results folder and create a csv file in `output_raw_results` with the results of the predictions (ipTM,pLDDT,pTM). It will also print out well-scoring models and some summary information.
3. `modify_mmcif_plddt.py` - This script will modify the output prediction (mmcif) files to smooth out the pLDDT scores by residue and with an optional sliding window for further downstream analysis. Creates a new modified mmcif file in the same directory as the original output file.
4. `get_all_binding_domains.py` - This script will take the modified mmcif files and use ChimeraX to determine the number of contacts and buried area between the prey protein and various domains of the bait protein. In doing so, we estimate the binding domain of the prey protein to the bait protein.
5. `calculate_all_clashes.py` - This script will take the modified mmcif files and use ChimeraX to determine the number of clashes between the prey protein and the LRRK2 RCKW (inflexible) region. This is used a filtering criterion to determine if the prey protein is a good candidate for further analysis.
6. `merge_results.py` - This script will take the output from the previous scripts (`process_results.py`, `get_all_binding_domains.py`, and `calculate_all_clashes.py`) and merge them into a single csv file.
7. `analyze_results.py` - This script will take the merged csv file and filters the results based on the ipTM and clashes scores, creating a list of good candidates and plots for visualization.


## Directory Breakdown


For questions or issues, please contact Daniel at daji@ucsd.edu.