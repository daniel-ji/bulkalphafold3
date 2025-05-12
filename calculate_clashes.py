# install chimerax for calculating clashes if needed

import os
import subprocess
import re
import argparse


def calculate_clashes(input_file, reference_file):
    """
    Calculate protein structure clashes between input and reference files using ChimeraX.
    Returns the number of unique clashes between the structures.
    """
    clash_script_file = os.path.basename(input_file) + "_clash_script.temp.cxc"

    clash_script = f"""
open ./{input_file}
delete #1/PREY @@bfactor<40
clashes ignoreHiddenModels true
open ./{reference_file}
match #1 to #2
select #1/BAIT
delete atoms sel
delete bonds sel
clashes ignoreHiddenModels true
hide #1 models
clashes ignoreHiddenModels true
exit
    """

    with open(clash_script_file, "w") as f:
        f.write(clash_script)

    command = "chimerax --nogui " + clash_script_file
    process = subprocess.run(command, shell=True, capture_output=True)

    # print(process.stdout.decode().strip())
    matches = re.findall(r"STATUS:\n([0-9]+ clashes\n|No clashes)", process.stdout.decode().strip())
    if len(matches) != 3:
        print(f"ERROR: Clash calculation failed for {input_file} and {reference_file}")
        return (-1, -1, -1, -1)

    os.remove(clash_script_file)

    input_only_clashes = 0 if "No clashes" in matches[0] else int(matches[0].split()[0])
    reference_only_clashes = 0 if "No clashes" in matches[2] else int(matches[2].split()[0])
    both_clashes = 0 if "No clashes" in matches[1] else int(matches[1].split()[0])

    between_clashes = both_clashes - input_only_clashes - reference_only_clashes
    print(f"Clashes between {input_file} and {reference_file}: {between_clashes}", flush=True)

    return (input_only_clashes, reference_only_clashes, both_clashes, between_clashes)


def main():
    parser = argparse.ArgumentParser(description="Calculate protein structure clashes using ChimeraX")
    parser.add_argument("input", help="Input model file")
    parser.add_argument("reference", help="Reference model file")

    args = parser.parse_args()

    print("Calculating clashes...")
    _, _, _, clash_count = calculate_clashes(args.input, args.reference)
    print(f"Clashes between {args.input} and {args.reference}: {clash_count}")

    return args.input, args.reference, clash_count


if __name__ == "__main__":
    main()
