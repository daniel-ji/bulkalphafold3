# BulkAlphaFold3

A set of pipelines to run AF3 PPI screens and individual protein complex screens and do downstream analysis of the results. Developed under the Leschziner Lab.

**Still in development.** This is a work in progress and not a complete solution (~80% complete). Please check back for updates. Documentation is also a work in progress.

Planned features:
- Create a main script to run all scripts in order, with command line arguments.
- Parameterize constants & hardcoded values in the scripts.
- Add more documentation and examples.

## Requirements
- AlphaFold3
- Python 3.8+
- BioPython
- numpy
- pandas
- matplotlib
- seaborn

## Setting Up AlphaFold3 

For BulkAlphaFold3, we used the Docker container to run the predictions and have a guide on [setting up the environment](guides/setting_up_environment.md) specifically for LambdaLabs, a cloud computing service that provides powerful GPUs for running AlphaFold3 predictions (we are not affiliated with LambdaLabs, but we found it to be a convenient service for running AF3 predictions since we did not have access to a powerful enough GPU). For AlphaFold3 setup in general, see the [AlphaFold3 documentation](https://github.com/google-deepmind/alphafold3). 

## Running BulkAlphaFold3

See the [guides](guides/) directory for more detailed guides on running the two primary use cases of BulkAlphaFold3: running a PPI complex screen with one bait protein and multiple prey proteins, and running individual protein complex screen (predicting all possible combinations of proteins in a complex).

Pulldown screen guide: [Running a PPI Complex Screen with BulkAlphaFold3](guides/running_pulldown_screens.md)

Individual protein complex screen guide: [Running Single Protein Complex Prediction Screens](guides/running_complex_screens.md)

For questions or issues, please contact Daniel at daji@ucsd.edu.
