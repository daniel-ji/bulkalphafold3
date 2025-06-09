"""
This script generates all possible combinations of binding domains for a given set of binding domains.
For each prey protein, it will create a new column for each combination of binding domains,
and assign a value of 1 if the combination is present in the binding domains of that prey protein (the number of contacts surpass a threshold).
For example, if the binding domains are "A", "B", and "C", the script will generate the following combinations:
- A, B, C, AB, AC, BC, ABC
"""

import itertools
import pandas as pd
from constants import CURRENT_PIPELINE
assert CURRENT_PIPELINE == "pulldown", "This script is only for the pulldown pipeline"
from constants import DOMAINS_TO_RESIDUES, FOLDERS, PLDDT_SLIDING_WINDOW, MIN_CONTACTS_THRESHOLDS


def generate_binding_domain_combinations(binding_domains):
    """
    Generate all possible combinations of binding domains.
    """
    combinations = []
    for r in range(1, len(binding_domains) + 1):
        combinations.extend(itertools.combinations(binding_domains, r))

    # reverse the order of combinations to get the most specific ones first
    return list(reversed(combinations))


def process_binding_domain_combinations(folder):
    """
    Process a folder to generate binding domain combinations.
    """
    binding_domain_combinations = generate_binding_domain_combinations(DOMAINS_TO_RESIDUES.keys())
    # append "_contacts" to get df column names
    binding_domain_combinations_columns = [tuple(f"{domain}_contacts" for domain in combination) for combination in binding_domain_combinations]

    binding_domain_file = f"output_binding_domain/{folder}_binding_domain_plddt_window_{PLDDT_SLIDING_WINDOW}.csv"

    binding_domain_df = pd.read_csv(binding_domain_file)

    # create new columns for each combination of binding domains
    for i, combination in enumerate(binding_domain_combinations_columns):
        for threshold in MIN_CONTACTS_THRESHOLDS:
            column_name = "_".join(binding_domain_combinations[i]) + f"_MIN_{threshold}"
            # check if all domains in the combination meet the threshold for each row
            mask = (binding_domain_df[list(combination)] >= threshold).all(axis=1)
            # and check that any of the (superset) domains with this domain as a subset do not meet the threshold (keep columns mutually exclusive)
            for j in range(i):
                superset_combination = binding_domain_combinations_columns[j]
                if set(combination).issubset(set(superset_combination)):
                    mask &= (binding_domain_df[list(superset_combination)] < threshold).any(axis=1)
            binding_domain_df[column_name] = mask.astype(int)

    binding_domain_df.to_csv(binding_domain_file, index=False)


if __name__ == "__main__":
    for folder in FOLDERS:
        print(f"Processing {folder}", flush=True)
        process_binding_domain_combinations(folder)
