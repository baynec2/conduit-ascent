# Conduit: A Modular Metaproteomics Analysis Platform

Conduit is a scalable and modular workflow management system for metaproteomics data analysis, designed to seamlessly integrate with metagenomic data if available. Built using Snakemake, it provides a robust pipeline for processing Data Independent Acquisition (DIA) mass spectrometry data with particular emphasis on metaproteomics applications.

## Features

*Note: Conduit is currently a work in progress and all features are not yet available.* 

- **Search Space Definition**: Choose the best way to define the search space for *your* experiment.
    - **User Defined Taxa**: Know what taxa are in your sample? Great! Just give Conduit the NCBI taxa IDs, we will do the rest.
    **Works in Progress:** 
    - **User Defined Proteomes**: Know the specific proteome IDs in your sample? Also amazing, we can work with those too. 
    - **Reference Based Metagenomic**: Don't know what taxa are in your sample, but were able to get metagenomic sequencing data? No problem! We can use it to define the search space.
    - **MAGs**: Have metagenomic sequencing data, but reference based metagenomic data not specific enough for your microbial community? No problem, you create the MAGs and we will then use those to define the search space.
    - **16S**: Don't know what is in your sample, can't get your hands on any metagenomic sequencing data, but have 16S sequencing data? Definitely the worst option, but probably better than nothing. We will do our best to use it.
- **Metadata Handling**: Easily associate your sample data with the information about what each sample actually is. Less time spent on data wrangling means more time to figure out what it all means.
- **DIA-NN Integration**: Automated, spectral library free processing of DIA data.
- **Taxonomic Analysis**: Multi-level taxonomic classification of proteins
- **Functional Analysis**: GO term, KEGG pathway, and subcellular location annotation.
- **R Integration**: Direct integration with R to enable easy statistical analysis, plotting, and more. 
- **Integration with a Dedicated GUI**: Want to explore your data quickly without having to write any code? Conduit-GUI has you covered. 
- **Container Support**: Full containerization via Apptainer. Run it on any reasonable Linux machine with ease.
- **Scalable**: Conduit runs on a single machine or scales to HPC clusters. Tip: metaproteomics loves compute — the more, the better.
- **Open Source**: Conduit is open source, allowing users to customize the pipeline to their specific needs.

## Dependencies

### Required Software
- Snakemake (≥7.0.0)
- R
- conduitR
- Python libraries
   * glob
   * pandas 
   * shutil
   * pdb  
   * logging
   * datetime
- DIA-NN (2.1.0)
- Apptainer/Singularity (≥1.1.0)

Note: All other dependencies (R, Python, DIA-NN, etc.) are automatically handled through Apptainer containers. You only need to ensure that Snakemake and Apptainer are installed on your *Linux* system.

If using Windows (not recommended), you can still use Conduit without containers.

Conduit will not run on MACOS.

### Hardware Requirements
- 64GB+ RAM recommended.
- Multi-core processor (8+ cores recommended)
- Sufficient storage space for raw data and results. Raw files from an Astral MS are usually around ~7 GB each, so it gets to TB scale very quickly.

## Quick Start

1. Make sure you have the required dependencies:

**If you are using Apptainer (Recommended but only works on Linux):**

You will need to install Snakemake and Apptainer. I recommend using conda to manage your dependencies.

```bash
conda create -n conduit
conda activate conduit
conda install -c bioconda snakemake
conda install -c bioconda apptainer
```

**If you are not using Apptainer (Not recommended but will work on Windows or Linux):**  

You will need to install the following dependencies manually:  

- Snakemake (≥7.0.0)
- R
- conduitR
- Python libraries
   * glob
   * pandas 
   * shutil
   * pdb  
   * logging
   * datetime
- DIA-NN (2.1.0)
- Apptainer/Singularity (≥1.1.0)

2. Clone the repository:
```bash
git clone https://github.com/baynec2/conduit.git
cd conduit
```

3. Set up an experiment directory following the example structure:

```
experiments/your_experiment/
├── input/
│   ├── raw_files/        # Your .raw files
│   ├── ncbi_taxa_ids.txt     # NCBI organism IDs
│   └── sample_annotation.txt
└── config/
    └── snakemake.yaml    # Experiment configuration

Your config file should be a YAML file named `snakemake.yaml` that looks like the following:

```yaml
# Define the name of the experiment (this will be used to find the correct experiment directory)
experiment: 
# Config file to generate the spectral library for DIA-NN (relative to main Snakefile)
generate_diann_spectral_library_config: config/generate_diann_spectral_library.cfg
# Config file to run DIA-NN search (relative to main Snakefile)
run_diann_config: config/run_diann.cfg
# Define the method used to define the search space. 
# Options: 
# - "ncbi_taxonomy_id": Use NCBI taxonomy IDs to define the search space.
# - "uniprot_proteome_id": Use UniProt proteome IDs to define the search space.
# - "MAG": Use MAGs to define the search space.
# - "metagenomic_profiling": Use metagenomic profiling to define the search space.
# - "16S": Use 16S rRNA gene sequences to define the search space.
search_space_method: ncbi_taxonomy_id
# Define the sample annotation filepath relative to the experiment directory
sample_annotation: input/sample_annotation.txt
# Define the output directory relative to the experiment directory
output_dir: output/
```

Make sure that the experiment in the config file matches the name of the experiment directory. A template is provided in the top level config folder.

4. If you plan to use Apptainer, you will need to build the diann and conduitR apptainer images. This only needs to be done once, future runs of Conduit will automatically use the images you have built. You can build the images by running the following commands:

```bash
snakemake build_conduitR_apptainer build_diann_apptainer --configfile experiments/your_experiment/config/snakemake.yaml
```

5. Run the workflow:
```bash
snakemake --configfile experiments/your_experiment/config/snakemake.yaml --use-apptainer
```

6. Use `experiments/your_experiment/output/conduit.rds` as the input to [Conduit-GUI](https://github.com/baynec2/conduit-GUI) to explore your data. Or use [conduitR](https://github.com/baynec2/conduitR) if you prefer a more bespoke data analysis experience.

7. Repeat steps 3-5 for each experiment you want to run (keeping a reproducible record of your experiments in the experiments subdirectory).

## Project Structure

Below you will find a high level annotated diagram displaying the general structure of Conduit:

```
conduit/
├── modules/                          # Snakemake modules. Each module should have snakefile and associated scripts. 
│   ├── search_space/                 # Defining the search space
│   │   ├── user_specified            # with ncbi_taxonomy and proteome_id 
│   │   ├── proteotyping              # with species specific peptides
│   │   ├── metagenomic_profiling     # with metagenomic profiling (reference based)
│   │   ├── mags                      # with metagenome assembled genome
│   │   └── 16s                       # with 16S data
│   ├── diann/                        # Identification and Quantification with DIA-NN
│   ├── annotation/                   # Protein / taxonomic annotation
│   │   ├── user_specified            # with ncbi_taxonomy and proteome_id 
│   │   ├── proteotyping              # with species specific peptides
│   │   ├── metagenomic_profiling     # with metagenomic profiling (reference based)
│   │   ├── mags                      # with metagenome assembled genome
│   │   └── 16s                       # with 16S data
│   ├── matrices/                     # Matrix processing
│   └── r_integration                 # Integration into R
├── apptainer/                        # Apptainer configurations. Contains .def files
├── config/                           # Default configurations (copied to each experiment config if not altered)
├── images/                           # Contains diagrams of the workflow or any other images associated with Conduit.
├── experiments/                      # Experiment directories
│   └── example/                      # Example experiment (for additional experiments, create a new directory and config file.)
│       ├── config                    # Configuration files for example experiment
│       ├── logs                      # Log files for each rule. 
│       ├── output                    # Outputs from the workflow
│       └── input                     # Inputs into the workflow
│           ├── sample_annotation.txt # Sample annotation file matching names of .raw MS files. 
│           ├── ncbi_taxa_ids.txt     # Text file NCBI organism IDs defining the search space (if using ncbi_taxonomy)
│           ├── proteome_id.txt       # Text file containing proteome ids (if using proteome id workflow)
│           ├── fastq                 # Fastq input for metagenomic profiling or 16S
│           ├── mags                  # MAG inputs.
│           └── raw_files             # Thermo Mass Spectrometry .raw files. 
└── tests/                            # Test suite
```

## Configuration

The main configuration is done through `snakemake.yaml`. Key parameters include:

```yaml
# Define the name of the experiment (this will be used to find the correct experiment directory)
experiment: 
# Config file to generate the spectral library for DIA-NN (relative to main Snakefile)
generate_diann_spectral_library_config: config/generate_diann_spectral_library.cfg
# Config file to run DIA-NN search (relative to main Snakefile)
run_diann_config: config/run_diann.cfg
# Define the method used to define the search space. 
# Options: 
# - "ncbi_taxonomy_id": Use NCBI taxonomy IDs to define the search space.
# - "uniprot_proteome_id": Use UniProt proteome IDs to define the search space.
# - "MAG": Use MAGs to define the search space.,
# - "metagenomic_profiling": Use metagenomic profiling to define the search space.
# - "16S": Use 16S rRNA gene sequences to define the search space.
search_space_method: ncbi_taxonomy_id
# Define the sample annotation filepath relative to the experiment directory
sample_annotation: input/sample_annotation.txt
```
You can use the configuration file in `config/snakemake.yaml` as a template file for your experiments.

## Output Structure

The workflow generates several key outputs:

1. **Database Resources**
   - Protein FASTA files
   - Taxonomy information

2. **DIA-NN Results**
   - Spectral library
   - Protein/peptide quantification

3. **Detected Protein Annotations**
   - GO term annotations
   - KEGG pathway annotations
   - Subcellular location annotations

4. **Processed Matrices**
   - Taxonomic level matrices
   - Functional annotation matrices
   - Protein-level matrices
   - Peptide-level matrices
   - Precursor-level matrices

5. **R Objects**
   - QFeatures object
   - Analysis metrics
   - Conduit object (accepted by Conduit-GUI/conduitR)

## Troubleshooting

### Common Issues

1. **Issues getting NCBI taxonomy**: If you are having issues getting taxonomy with the ncbi_taxonomy_id search space method, try waiting for a bit and then retrying. I have seen that occasionally the NCBI API will be down.

2. **Configuration Errors**: Make sure your `snakemake.yaml` file is properly formatted and all required fields are filled in.

### Getting Help

- Check the [Issues](https://github.com/baynec2/conduit/issues) page for known problems
- Create a new issue with detailed error messages and your configuration
- Include relevant log file from the `logs/` directory when reporting problems

## Contributing

Contributions are welcome! Please feel free to submit an issue or pull request.

### Contributing New Modules

We welcome contributions of new modules to expand Conduit's capabilities! Here's how to contribute:

#### New Conduit Module Development Guide

##### 1. Module Structure and Conventions

All Snakemake modules in Conduit should adhere to the following structure:

- **Directory Naming:**  
  Each module should reside in its own directory describing the task it performs (e.g., `search_space/`, `diann/`, `annotation/`).

- **Alternative Modules:**  
  If multiple alternative modules exist for accomplishing the same task, they should be placed in subdirectories under the main task directory.  
  *Example:*  
  ```
  defining_search_space/
        ├── user_specified            # with ncbi_taxonomy and proteome_id 
        ├── proteotyping              # with species specific peptides
        ├── metagenomic_profiling     # with metagenomic profiling (reference based)
        ├── mags                      # with metagenome assembled genome
        └── 16s                       # with 16S data
  ```
- **Snakemake file:**  
  Each module must contain a snakemake file with its rules. This should be named `module_name.smk`.

- **Scripts:**  
  A `scripts/` subdirectory should contain any scripts required by that module.

- **Logging:**  
  Logs should be written to the main experiment's `logs/` directory.

- **No `rule all` in Modules:**  
  `rule all` should only be defined in the main Conduit Snakefile, never within modules.

##### 2. Search Space and Annotation Modules

Each search space module **must have a corresponding annotation module**, because different approaches to defining the search space will require different strategies for annotating proteins.

##### 3. Required Outputs: Search Space Modules

Each search space module must produce the following files, located in a folder named `database_resources/`.  
If a particular file cannot be generated (e.g., unclear what to use as a proteome_id), it must still be produced with `NA` values as appropriate.

**Required Output Files:**

| Filename                          | Contents                                                                 |
|------------------------------------|--------------------------------------------------------------------------|
| `database.fasta`                | Fasta file containing all protein sequences to be included in search space with UniProt-style headers. |
| `proteome_ids.txt`              | .txt file containing the proteome IDs that were used in the experiment.  |
| `taxonomy.txt`                  | .txt file with all taxonomy information in the sample.                   |
| `protein_info.txt`              | .txt file containing protein information.                                |
| `taxonomic_tree_of_database.pdf`| PDF file showing an image of what the search space looks like from a taxonomic lens. |
| `database.predicted.speclib`    | Spectral library produced by DIA-NN (automatically produced via the shared diann module). |
| `README.md`                        | Readme containing metrics about the database resources.                  |

##### 4. Required Outputs: Annotation Modules
Each annotation module must produce the following files. If a particular annotation type cannot be retrieved (e.g., no subcellular prediction for MAGs), the file must still be generated with appropriate placeholder content.

| Filename                        | Contents                                                                                   |
|----------------------------------|--------------------------------------------------------------------------------------------|
| `detected_protein_info.txt`      | Contains `protein_info.txt` filtered to only contain the proteins that were actually detected by DIA-NN. (Automatically generated via the shared diann module.) |
| `detected_protein.fasta`         | The content from `detected_protein_info.txt` in fasta format. (Automatically generated via the shared diann module.) |
| `annotated_protein_info.txt`     | Adds annotations to the detected proteins.                                                  |
| `go_annotations.txt`             | File containing GO annotations in long format. Each protein should have its GO terms listed.|
| `subcellular_locations.txt`      | File containing subcellular location predictions.                                           |
| `kegg_annotations.txt`           | File containing KEGG pathway annotations.                                                   |

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

Conduit would not be possible without the great work of many people. 

- Snakemake developers
- DIA-NN developers
- R Bioconductor community
- R tidyverse community
- R developers
- UniProt consortium
- NCBI taxonomy database maintainers
- KEGG database maintainers
- Gene Ontology consortium

