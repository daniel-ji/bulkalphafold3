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
from constants import FOLDERS, DOMAINS_TO_RESIDUES, PLDDT_SLIDING_WINDOW, USE_MULTIPROCESSING, PROCESS_COUNT
from get_binding_domain_combinations import process_binding_domain_combinations

OUTPUT_FOLDER = "output_binding_domain/"

def get_binding_domain_stats(model_file):
    rand_suffix = os.urandom(4).hex()

    contacts_script_file = f"{os.path.basename(model_file)}_{rand_suffix}_contacts_script.temp.cxc"
    contacts_script = f""" 
open {model_file} 
delete #1/PREY @@bfactor<40 
sel #1/PREY
contacts sel restrict both
~sel
"""

    for residues in DOMAINS_TO_RESIDUES.values():
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
    # the next len(DOMAINS_TO_RESIDUES) * 2 matches are the contacts within the bait protein itself and the bait protein + prey protein
    stdout = stdout.replace("INFO:\nNo contacts", "INFO:\n0 contacts")
    contacts_matches = re.findall(r"INFO:\n([0-9]+ contacts)\n", stdout)

    os.remove(contacts_script_file)

    failed = False

    # get contacts count
    interdomain_contacts = []
    if len(contacts_matches) == 1 + len(DOMAINS_TO_RESIDUES) * 2:
        prey_protein_contacts = int(contacts_matches[0].split()[0])
        for i in range(1, len(contacts_matches), 2):
            only_domain_contacts = int(contacts_matches[i].split()[0])
            domain_and_prey_contacts = int(contacts_matches[i + 1].split()[0])
            interdomain_contacts.append(domain_and_prey_contacts - only_domain_contacts - prey_protein_contacts)
    else:
        failed = True

    # get buried area
    buriedareas = []
    if len(buriedarea_matches) == len(DOMAINS_TO_RESIDUES):
        buriedareas = [float(match[2]) for match in buriedarea_matches]
    else:
        failed = True

    # edge case: when the entire structure has a sliding window plDDT < 40, the program will fail, but just return 0 for the buried area and contacts
    PLDDT_TOO_LOW_MATCH = re.search(r"Executing: select #1/PREY\nINFO:\nNothing selected", stdout)
    if PLDDT_TOO_LOW_MATCH:
        buriedareas = [0 for _ in DOMAINS_TO_RESIDUES.keys()]
        interdomain_contacts = [0 for _ in DOMAINS_TO_RESIDUES.keys()]
        failed = False

    output_text = f"Processing {model_file}\n"
    if failed:
        output_text += "=================================================\n"
        output_text += f"Failed to calculate buried area or contacts for {model_file}\n"
        output_text += stdout + "\n"
        output_text += process.stderr.decode().strip() + "\n"
        output_text += "=================================================\n"
        print(output_text, flush=True)
        return [model_file] + [None for _ in DOMAINS_TO_RESIDUES.keys()] + [None for _ in DOMAINS_TO_RESIDUES.keys()]
    else:
        result = list(zip(DOMAINS_TO_RESIDUES.keys(), buriedareas, interdomain_contacts))
        output_text += f"{result}\n"
        print(output_text, flush=True)
        return [model_file] + [area for _, area, _ in result] + [contacts for _, _, contacts in result]


def process_model_parallel(model_files):
    with mp.Pool(processes=PROCESS_COUNT) as pool:
        results = pool.map(get_binding_domain_stats, model_files)

    return results


if __name__ == "__main__":
    if not os.path.exists(OUTPUT_FOLDER):
        os.makedirs(OUTPUT_FOLDER)

    for folder in FOLDERS:
        model_files = get_model_files(folder, residue_sliding_window=PLDDT_SLIDING_WINDOW)
        if USE_MULTIPROCESSING:
            results = process_model_parallel(model_files)
        else:
            results = [get_binding_domain_stats(model_file) for model_file in model_files]
        binding_domain_area = [key + "_area" for key in DOMAINS_TO_RESIDUES.keys()]
        binding_domain_contacts = [key + "_contacts" for key in DOMAINS_TO_RESIDUES.keys()]
        output_file = f"binding_domain_plddt_window_{PLDDT_SLIDING_WINDOW}.csv"
        result_df = pd.DataFrame(results, columns=["model"] + binding_domain_area + binding_domain_contacts)
        result_df.to_csv(f"{OUTPUT_FOLDER}{folder}_{output_file}", index=False)
        process_binding_domain_combinations(folder)
        
        print(f"Saved {folder}_{output_file}")
    print("Done!\n")
