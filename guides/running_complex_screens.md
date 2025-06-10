# Running Single Protein Complex Prediction Screens

As an additional feature of BulkAlphaFold3 aside from running AlphaFold3 predictions on a large set of prey proteins against a bait protein, we also provide scripts to run individual protein complex screens, predicting all possible combinations of proteins in a complex. This is primarily intended to see what AlphaFold3 predicts in terms of a protein complex forming between multiple proteins. It is worth noting that AlphaFold3 is not necessarily robust in predicting protein complexes (especially for large complexes or complexes with many proteins), but it can still potentially provide useful insights into possible interactions.

This is done through the `generate_protein_complex_combinations.py` script. The script takes as input a configuration file that defines various parameters for the protein complex screen as well as a FASTA file that contains the sequences of the proteins in the complex. Notably, the configuration file can optionally define a fixed number of proteins that are always included in any prediction, which is useful for the case you desire to have multiple proteins always as part of the complex. There is also an option to define custom predictions, which allows you to specify specific combinations of proteins to predict, in addition to the automatically generated combinations based on the fixed proteins and the maximum combined sequence length.

The script will generate AlphaFold3 input json files for running predictions all possible combinations of the proteins in the complex (constrained on always including the fixed number of proteins, if specified). It will also create a bash one-liner script to run the predictions via Docker, intended to be run on the environment that has been setup by following the [setting up environment guide](setting_up_environment.md).

An example can be seen in the `configs/complex/PxdA.yml` configuration file and `examples/complex/PxdA.fasta`, which contains a complex with the PxdA protein and six other proteins it potentially forms a complex with:

```yaml
current_pipeline: complex
number_of_seeds: 5
max_combined_seq_length: 10000
process_count: 10

fixed_proteins:
- PxdA
- AN4585
- AN3602
- AN4965

custom_predictions:
- [AN4585, AN3602, AN4965, AN3669, AN5898, AN7628]

input_fasta: examples/complex/PxdA.fasta
output_folder: output_protein_complex/PxdA
```

After modifying the configuration file to your needs, and setting the configuration file in `constants.py`, you can just run the script with `python generate_protein_complex_combinations.py`. This will generate the AlphaFold3 input json files in the specified output folder, as well as a bash script to run the predictions.