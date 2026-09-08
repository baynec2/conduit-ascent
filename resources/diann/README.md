# Put your DIA-NN installation here

Conduit does not distribute DIA-NN. Its licence permits **one copy for backup
purposes** and forbids renting, leasing, lending, or sublicensing, so the binary
cannot ship inside a Conduit container image. You download it yourself and
accept its terms directly — free for academic use; commercial use requires a
licence from the authors.

    https://github.com/vdemichev/DiaNN/releases
    (DIA-NN-<version>-Academia-Linux.zip)

Unzip it here and make the CLI executable:

    unzip DIA-NN-2.5.0-Academia-Linux.zip -d resources/diann
    chmod +x resources/diann/diann-2.5.0/diann-linux

That produces `resources/diann/diann-2.5.0/`, which is the default value of
`diann_path` in `config/snakemake.yaml` — nothing else to configure. The
directory is bind-mounted into the `diann_runtime` image (Debian/Ubuntu + .NET 8
+ libgomp, and no DIA-NN) and `diann-linux` is executed from the mount.

Installing somewhere else is fine; point `diann_path` at it and the workflow
adds the bind for you:

    export CONDUIT_DIANN_PATH=/opt/diann-2.5.0

Tested versions: 2.3.0 and 2.5.0. Below 2.2 is rejected — those builds silently
ignore `--pre-search` / `--pre-filter`.

See "Obtaining DIA-NN" in the repository README for the full story, including
the two other supported shapes (a self-contained executable, or your own
container image).

The contents of this directory are gitignored apart from this file.
