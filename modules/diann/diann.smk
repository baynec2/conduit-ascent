import glob
import os
import sys

EXPERIMENT_DIR = config["experiment_dir"]
RUN_DIR = config["run_dir"]

sys.path.insert(0, os.path.join(workflow.basedir, "modules", "_shared"))
from diann_staging import (
    list_raw_files,
    list_samples,
    raw_path_for_sample,
    stage3_symlink_commands,
)

RAW_FILEPATHS = list_raw_files(EXPERIMENT_DIR)
SAMPLES = list_samples(EXPERIMENT_DIR)

DIANN_OUT = os.path.join(RUN_DIR, "diann_output")
DIANN_QUANTS = os.path.join(DIANN_OUT, "quant_files")
#################################################################################
# Generating Spectral Library
#################################################################################
rule generate_diann_spectral_library:
    input:
        fasta = os.path.join(RUN_DIR,"database_resources/database.fasta"),
        config_file = os.path.join(RUN_DIR,"config/diann_spectral_library_base.cfg")
    output:
        os.path.join(RUN_DIR,"database_resources/database.predicted.speclib")
    params:
        out_lib = lambda w, output: os.path.splitext(os.path.splitext(output[0])[0])[0]
    log: os.path.join(RUN_DIR,"logs/diann/generate_diann_spectral_library.log")
    container:
        config["containers"]["diann"]
    threads: workflow.cores
    shell:
        """
        diann --cfg {input.config_file} \
        --fasta {input.fasta} \
        --out-lib {params.out_lib} \
        --met-excision \
        --cut "K*,R*" \
        --missed-cleavages 1 \
        --min-pep-len 7 \
        --max-pep-len 30 \
        --threads {threads} >> {log} 2>&1
        """
################################################################################
# Running DIANN
#
# Standard mode: three-stage manual MBR split (library-gen, per-raw, combine).
# Without --reanalyse, stage 1 is single-pass + library-write; stages 2 and 3
# parallelise the per-file conventional search and combine. Wins on resumability,
# incremental sample addition, and cheap iteration on quant-only settings.
#
# InfinDIA mode: monolithic. DIA-NN's --pre-search internally runs all three
# phases (pre-search, first pass, second pass) and there is no documented or
# undocumented flag to stop after the empirical library is written (verified
# empirically against DIA-NN 2.5.0 on 2026-05-08). Splitting into 3 stages would
# pay the second-pass cost twice — once internally, once in stage 2 — for no
# scientific benefit. So we let DIA-NN produce the final report directly.
################################################################################
if config.get("diann_search_mode", "standard") == "standard":

    rule run_diann_build_empirical_lib:
        input:
            raw_files_dir = os.path.join(EXPERIMENT_DIR,"input/ms_files"),
            spectral_library = os.path.join(RUN_DIR,"database_resources/database.predicted.speclib"),
            fasta = os.path.join(RUN_DIR,"database_resources/database.fasta"),
            config_file = os.path.join(RUN_DIR,"config/run_diann.cfg")
        output:
            empirical_lib = os.path.join(DIANN_OUT, "empirical.parquet")
        params:
            out_lib = os.path.join(DIANN_OUT, "empirical"),
        log: os.path.join(RUN_DIR,"logs/diann/run_diann_build_empirical_lib.log")
        container:
            config["containers"]["diann"]
        threads: workflow.cores
        shell:
            """
            mkdir -p $(dirname {log}) $(dirname {output.empirical_lib})
            # --gen-spec-lib + --out-lib gives us the empirical library;
            # --rt-profiling writes empirically-aligned RTs into it (quality
            # knob, no extra pass). No --reanalyse: it would only add a
            # second pass that re-searches the raws against this library,
            # exactly what stage 2 does per-file.
            diann --cfg {input.config_file} \
            --fasta {input.fasta} \
            --dir {input.raw_files_dir} \
            --lib {input.spectral_library} \
            --gen-spec-lib \
            --rt-profiling \
            --out-lib {params.out_lib} \
            --threads {threads} --verbose 1 >> {log} 2>&1
            """

    rule run_diann_search_one_raw:
        input:
            empirical_lib = os.path.join(DIANN_OUT, "empirical.parquet"),
            fasta = os.path.join(RUN_DIR,"database_resources/database.fasta"),
            config_file = os.path.join(RUN_DIR,"config/run_diann.cfg"),
            raw = lambda w: raw_path_for_sample(EXPERIMENT_DIR, w.sample)
        output:
            quant = os.path.join(DIANN_QUANTS, "{sample}.quant")
        params:
            tmpdir = lambda w: os.path.join(DIANN_OUT, "quant_tmp", w.sample)
        log: os.path.join(RUN_DIR,"logs/diann/run_diann_search_one_raw.{sample}.log")
        container:
            config["containers"]["diann"]
        threads: min(8, workflow.cores)
        shell:
            """
            mkdir -p $(dirname {log})
            rm -rf {params.tmpdir} && mkdir -p {params.tmpdir} $(dirname {output.quant})
            diann --cfg {input.config_file} \
            --f {input.raw} \
            --lib {input.empirical_lib} \
            --fasta {input.fasta} \
            --temp {params.tmpdir} \
            --out {params.tmpdir}/per_run_report \
            --threads {threads} --verbose 1 >> {log} 2>&1
            mv {params.tmpdir}/*.quant {output.quant}
            rm -rf {params.tmpdir}
            """

    rule run_diann_combine:
        input:
            quants = expand(
                os.path.join(DIANN_QUANTS, "{sample}.quant"),
                sample=SAMPLES
            ),
            empirical_lib = os.path.join(DIANN_OUT, "empirical.parquet"),
            fasta = os.path.join(RUN_DIR,"database_resources/database.fasta"),
            config_file = os.path.join(RUN_DIR,"config/run_diann.cfg"),
            raw_files_dir = os.path.join(EXPERIMENT_DIR,"input/ms_files")
        output:
            diann_stats = os.path.join(DIANN_OUT, "diann.stats.tsv"),
            diann_parquet = os.path.join(DIANN_OUT, "diann.parquet")
        params:
            tmpdir = os.path.join(DIANN_OUT, "quant_combine_tmp"),
            out_prefix = os.path.join(DIANN_OUT, "diann"),
            symlink_cmds = stage3_symlink_commands(
                os.path.join(DIANN_OUT, "quant_combine_tmp"),
                RAW_FILEPATHS,
                DIANN_QUANTS
            )
        log: os.path.join(RUN_DIR,"logs/diann/run_diann_combine.log")
        container:
            config["containers"]["diann"]
        threads: workflow.cores
        shell:
            """
            mkdir -p $(dirname {log})
            rm -rf {params.tmpdir} && mkdir -p {params.tmpdir}
            {params.symlink_cmds}
            diann --cfg {input.config_file} \
            --dir {input.raw_files_dir} \
            --lib {input.empirical_lib} \
            --fasta {input.fasta} \
            --temp {params.tmpdir} \
            --use-quant \
            --out {params.out_prefix} \
            --threads {threads} --verbose 1 >> {log} 2>&1
            rm -rf {params.tmpdir}
            """

else:  # infinidia — monolithic, see comment above

    rule run_diann_monolithic:
        input:
            raw_files_dir = os.path.join(EXPERIMENT_DIR,"input/ms_files"),
            fasta = os.path.join(RUN_DIR,"database_resources/database.fasta"),
            config_file = os.path.join(RUN_DIR,"config/run_diann.cfg")
        output:
            diann_stats = os.path.join(DIANN_OUT, "diann.stats.tsv"),
            diann_parquet = os.path.join(DIANN_OUT, "diann.parquet")
        params:
            out_prefix = os.path.join(DIANN_OUT, "diann"),
            tmpdir = os.path.join(DIANN_OUT, "monolithic_quant_files"),
        log: os.path.join(RUN_DIR,"logs/diann/run_diann_monolithic.log")
        container:
            config["containers"]["diann"]
        threads: workflow.cores
        shell:
            # --temp keeps DIA-NN's per-run .pre.quant + .quant files inside the
            # run directory; without it DIA-NN writes them next to the raw files
            # in the user's input dir. We KEEP the .quant files after the run
            # completes — they're not just throwaway intermediates: a future
            # `--use-quant` invocation could re-combine them with different FDR /
            # normalisation settings without redoing the per-file search.
            """
            mkdir -p $(dirname {log}) $(dirname {output.diann_parquet})
            rm -rf {params.tmpdir} && mkdir -p {params.tmpdir}
            diann --cfg {input.config_file} \
            --fasta {input.fasta} \
            --dir {input.raw_files_dir} \
            --temp {params.tmpdir} \
            --pre-search --pre-filter \
            --gen-spec-lib \
            --rt-profiling \
            --out {params.out_prefix} \
            --threads {threads} --verbose 1 >> {log} 2>&1
            """
################################################################################
# Extracting Detected Proteins
################################################################################
rule extract_detected_proteins:
  input:
    protein_info_df=os.path.join(RUN_DIR,"database_resources/protein_info.txt"),
    protein_info_fasta =os.path.join(RUN_DIR,"database_resources/database.fasta"),
    diann_parquet=os.path.join(RUN_DIR,"diann_output/diann.parquet")
  output:
    detected_protein_info_df = os.path.join(RUN_DIR,"database_resources/detected_protein_resources/detected_protein_info.txt"),
    detected_protein_info_fasta = os.path.join(RUN_DIR,"database_resources/detected_protein_resources/detected_protein.fasta")
  log: os.path.join(RUN_DIR,"logs/diann/extract_detected_proteins.log")
  container: config["containers"]["conduitr"]
  script:
    "scripts/extract_detected_proteins.R"
