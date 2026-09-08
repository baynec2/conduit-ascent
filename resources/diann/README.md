# Put your DIA-NN installation here

Conduit does not distribute DIA-NN. Its licence permits **one copy for backup
purposes** and forbids renting, leasing, lending, or sublicensing, so the binary
cannot ship inside a Conduit container image. You download it yourself and
accept its terms directly — free for academic use; commercial use requires a
licence from the authors.

    https://github.com/vdemichev/DiaNN/releases
    (DIA-NN-<version>-Academia-Linux.zip)

Careful with that page: the newest *release* shown is 2.0 (January 2025), but that is
not the newest DIA-NN. Every build since is attached as an asset to that same `2.0`
tag, so the download URL always says `2.0` whatever version you want. Look at the
Assets list, not the release title.

    curl -fL -o /tmp/diann.zip \
      https://github.com/vdemichev/DiaNN/releases/download/2.0/DIA-NN-2.5.0-Academia-Linux.zip
    unzip /tmp/diann.zip -d resources/diann && rm /tmp/diann.zip
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

Conduit checks that whatever you install runs, accepts every flag it passes, and
writes the report columns it reads. Those are structural checks only: a new DIA-NN
release can change scoring, FDR or quantification and pass all of them. DIA-NN
releases often and we do not test every version, so pin one version for the duration
of a study and treat an upgrade as a change worth measuring.

See "Obtaining DIA-NN" in the repository README for the full story, including
the two other supported shapes (a self-contained executable, or your own
container image).

The contents of this directory are gitignored apart from this file.
