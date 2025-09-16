# Running Conduit on Barnacle2 (Knight Lab HPC)

These instructions describe how to run the **Conduit** workflow on Barnacle2 using Snakemake and Singularity for Knight lab members only. However, they may be adapted to other HPC systems using SLURM.
Note that if trying to run Conduit for a tutorial's sake you will also need files within directories of `conduit/experiments/example/input/database_resources` and `conduit/experiments/example/input/raw_files`. These files are a bit hefty, so please contact baynec2 directly for them.
---

## 1. Login to Barnacle2
Use your UCSD credentials to connect:

```bash
ssh <username>@barnacle2.ucsd.edu
```

## 2. Clone the repository and create a output directory for SLURM files

```bash
git clone https://github.com/baynec2/conduit.git
cd conduit
mkdir slurm_out
```
## 3. Install Snakemake in your BASE environment
Barnacle2 currently has Singularity available for all users at a system-wide level. For this reason, we want to use the Singularity within the base environment. In order to run Conduit, we will also need to install Snakemake (and any Snakemake-related dependencies, likely datrie and wheel which can be pip installed). Do NOT create and install this within a separate conda environment (as listed in the README). This will not work because we need to use the Singularity on Barnacle2. All other packages and environments are managed through the Singularity container.

```bash
# First check that your Singularity works for you on Barnacle2, currently on version 3.6.4 as of 9/16/25
singularity --version

# Pip install Snakemake and other needed dependencies
pip install snakemake
pip install wheel
pip install datrie
```

## 4. Add in the following lines to the beginning of `conduit/modules/annotation/ncbi_taxonomy/scripts get_annotations_from_uniprot.R`
Easiest/quickest way to do this is to `vim get_annotations_from_uniprot.R`, make your changes after typing `i`, and then `:x` and press enter to save and quit the file.

```
# Get Detected Proteins Annotation From Uniprot # previous line
################################################################################ # previous line

slurm_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", 1)) # *** add this line here ***
options(parallelly.maxWorkers.localhost = slurm_cores) # *** add this line here ***

# Open the log file to write both stdout and stderr # line continued
````

## 5. Download and put the SLURM script `run_conduit_barnacle2.sh` from the Conduit GitHub repo and place in your `conduit/` directory.
You will need to update --email in the SBATCH header with your own email. You may also change any notification preferences as well.

```bash
# Submit your job to run the Conduit Snakemake workflow
sbatch run_conduit_barnacle2.sh
# Check your job via
squeue --me
# Check the Snakemake logged outputs in slurm_out .err file
cat slurm_out/*.err
```


With --cpus-per-task=16 and --mem=64G (pre-set in the SBATCH header), the test inputs took approximately 50 minutes to run on Barnacle2

