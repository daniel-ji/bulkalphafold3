""" 
A script to get the predicted binding domain of a protein to LRRK2. Uses the chimerax `buriedarea` command to calculate the interacting region between the prey protein and LRRK2. 
""" 

import os 
import re 
import subprocess 
import matplotlib.pyplot as plt 
import pandas as pd 
import multiprocessing as mp
from files_helper import get_model_files
from constants import FOLDERS

PLDDT_SLIDING_WINDOW = 11
USE_MULTIPROCESSING = True  
PROCESS_COUNT = 11

# FULL LRRK2
# DOMAIN_TO_RESIDUES = { 
#     "ARM": [1, 704], 
#     "ANK": [705, 799], 
#     "LRR": [800, 1331], 
#     "ROC": [1332, 1521], 
#     "COR-A": [1522, 1671], 
#     "COR-B": [1672, 1878], 
#     "KINASE": [1879, 2139], 
#     "WD40": [2140, 2527] 
# } 
 
# for just ROC-COR 
DOMAIN_TO_RESIDUES = { 
    "ROC": [1317-1316, 1521-1316],
    # "COR": [1522-1316, 1873-1316],
    "COR-A": [1522-1316, 1671-1316], 
    "COR-B": [1672-1316, 1873-1316], 
    # "ALL": [1317-1316, 1873-1316] 
}

def get_binding_domain_area(model_file): 
    rand_suffix = os.urandom(4).hex()
 
    contacts_script_file = f"{os.path.basename(model_file)}_{rand_suffix}_contacts_script.temp.cxc" 
    contacts_script = f""" 
open {model_file} 
delete #1/PREY @@bfactor<40 
sel #1/PREY
contacts sel restrict both
~sel
""" 
 
    for residues in DOMAIN_TO_RESIDUES.values(): 
        contacts_script += f"measure buriedarea /PREY withAtoms2 /BAIT:{residues[0]}-{residues[1]}\n"
        contacts_script += f"sel #1/BAIT:{residues[0]}-{residues[1]}\n" 
        contacts_script += f"contacts sel restrict both\n"
        contacts_script += f"sel add #1/PREY\n"
        contacts_script += f"contacts sel restrict both\n"
        contacts_script += f"~sel\n"
                 
    contacts_script += "exit" 
 
    with open(contacts_script_file, "w") as f: 
        f.write(contacts_script) 
 
    command = "chimerax --nogui " + contacts_script_file 
    process = subprocess.run(command, shell=True, capture_output=True) 
 
    stdout = process.stdout.decode().strip()
    buriedarea_matches = re.findall(r"INFO:\nBuried area between /PREY and /BAIT:(\d+)-(\d+) = ([\d\.\-e]+)", stdout) 
    # first match is all of prey protein contacts
    # the next len(DOMAIN_TO_RESIDUES) * 2 matches are the contacts within the bait protein itself and the bait protein + prey protein
    stdout = stdout.replace("INFO:\nNo contacts", "INFO:\n0 contacts")
    contacts_matches = re.findall(r"INFO:\n([0-9]+ contacts)\n", stdout)

    os.remove(contacts_script_file)
    
    failed = False
    
    # get contacts count
    interdomain_contacts = []
    if len(contacts_matches) == 1 + len(DOMAIN_TO_RESIDUES) * 2:
        prey_protein_contacts = int(contacts_matches[0].split()[0])
        for i in range(1, len(contacts_matches), 2):
            only_domain_contacts = int(contacts_matches[i].split()[0])
            domain_and_prey_contacts = int(contacts_matches[i + 1].split()[0])
            interdomain_contacts.append(domain_and_prey_contacts - only_domain_contacts - prey_protein_contacts)
    else:
        failed = True
    
    # get buried area
    buriedareas = []
    if len(buriedarea_matches) == len(DOMAIN_TO_RESIDUES):
        buriedareas = [float(match[2]) for match in buriedarea_matches]
    else:
        failed = True

    # edge case: when the entire structure has a sliding window plDDT < 40, the program will fail, but just return 0 for the buried area and contacts
    PLDDT_TOO_LOW_MATCH = re.search(r"Executing: select #1/PREY\nINFO:\nNothing selected", stdout)
    if PLDDT_TOO_LOW_MATCH:
        buriedareas = [0 for _ in DOMAIN_TO_RESIDUES.keys()]
        interdomain_contacts = [0 for _ in DOMAIN_TO_RESIDUES.keys()]
        failed = False

    output_text = f"Processing {model_file}\n"
    if failed:
        output_text += "=================================================\n"
        output_text += f"Failed to calculate buried area or contacts for {model_file}\n"
        output_text += stdout + "\n"
        output_text += process.stderr.decode().strip() + "\n"
        output_text += "=================================================\n"
        print(output_text, flush=True)
        return [model_file] + [None for _ in DOMAIN_TO_RESIDUES.keys()] + [None for _ in DOMAIN_TO_RESIDUES.keys()]
    else:
        result = list(zip(DOMAIN_TO_RESIDUES.keys(), buriedareas, interdomain_contacts))
        output_text += f"{result}\n"
        print(output_text, flush=True)
        return [model_file] + [area for _, area, _ in result] + [contacts for _, _, contacts in result]
 
def process_model_parallel(model_files):
    with mp.Pool(processes=PROCESS_COUNT) as pool:
        results = pool.map(get_binding_domain_area, model_files)
    
    return results

if __name__ == "__main__":
    for folder in FOLDERS:
        model_files = get_model_files(folder, residue_sliding_window=PLDDT_SLIDING_WINDOW)
        if USE_MULTIPROCESSING:
            results = process_model_parallel(model_files)
        else:
            results = [get_binding_domain_area(model_file) for model_file in model_files]
        binding_domain_area = [key + "_area" for key in DOMAIN_TO_RESIDUES.keys()]
        binding_domain_contacts = [key + "_contacts" for key in DOMAIN_TO_RESIDUES.keys()]
        output_file = f"binding_domain_area_plddt_filtered_sliding_window_{PLDDT_SLIDING_WINDOW}.csv"
        result_df = pd.DataFrame(results, columns=["model"] + binding_domain_area + binding_domain_contacts)
        result_df.to_csv(f"{folder}_{output_file}", index=False)
        print(f"Saved {folder}_{output_file}")
    print("Done!\n")