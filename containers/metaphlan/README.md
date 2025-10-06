## MetaPhlAn Container

Container image providing [MetaPhlAn](https://github.com/biobakery/metaphlan) via Bioconda, with `bowtie2` and `samtools` preinstalled. A convenience entrypoint supports database installation.

### Build

```bash
docker build -f containers/metaphlan/Dockerfile -t metaphlan:4 containers/metaphlan
```

### Prepare a database directory (once)

```bash
mkdir -p metaphlan_db
docker run --rm -it \
  -v "${PWD}:/work" \
  -v "${PWD}/metaphlan_db:/db/metaphlan" \
  metaphlan:4 download-db
```

This downloads the MetaPhlAn index into `./metaphlan_db` on the host.

### Run MetaPhlAn

Examples (mount input/output working dir and the database dir):

```bash
# Profile a FASTQ (single-end)
docker run --rm -it \
  -v "${PWD}:/work" \
  -v "${PWD}/metaphlan_db:/db/metaphlan" \
  metaphlan:4 input.fastq.gz --input_type fastq --nproc $(nproc) -o profile.txt

# If you prefer to call metaphlan explicitly
docker run --rm -it \
  -v "${PWD}:/work" \
  -v "${PWD}/metaphlan_db:/db/metaphlan" \
  metaphlan:4 metaphlan input.fastq.gz --input_type fastq -o profile.txt
```

The container defaults to working in `/work`. Mount your current directory to `/work` for convenience.

### Notes
- Database path inside the image: `/db/metaphlan` (controlled by `METAPHLAN_DIR`/`METAPHLAN_BOWTIE2_DB`).
- Use `download-db` subcommand to install or update the database within the mounted directory.
- Show version: `docker run --rm metaphlan:4 version`.


