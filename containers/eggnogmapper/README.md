## eggNOG-mapper Container

Container image providing [eggNOG-mapper v2](https://github.com/eggnogdb/eggnog-mapper) via Bioconda, for sequence-based functional annotation using orthologous groups.

### Build

```bash
docker build -f containers/eggnogmapper/Dockerfile -t eggnogmapper:2.1.12 containers/eggnogmapper
```

### Download the eggNOG database (once)

```bash
mkdir -p eggnog_db
docker run --rm -it \
  -v "${PWD}/eggnog_db:/opt/eggnog_db" \
  eggnogmapper:2.1.12 \
  download_eggnog_data.py --data_dir /opt/eggnog_db -y
```

This downloads the full eggNOG database (~50GB) into `./eggnog_db` on the host. To download only the bacteria database (~8GB), add `--taxids 2`:

```bash
docker run --rm -it \
  -v "${PWD}/eggnog_db:/opt/eggnog_db" \
  eggnogmapper:2.1.12 \
  download_eggnog_data.py --data_dir /opt/eggnog_db --taxids 2 -y
```

### Run eggNOG-mapper

```bash
docker run --rm -it \
  -v "${PWD}:/work" \
  -v "${PWD}/eggnog_db:/opt/eggnog_db" \
  eggnogmapper:2.1.12 \
  emapper.py \
    -i proteins.fasta \
    --itype proteins \
    --output emapper_results \
    --data_dir /opt/eggnog_db \
    --cpu $(nproc)
```

### Notes

- Database path inside the container: `/opt/eggnog_db` (controlled by `EGGNOG_DATA_DIR`).
- Output file of interest: `emapper_results.emapper.annotations` — tab-separated with columns for eggNOG OGs, COG category, GO terms, KEGG KO/pathways, Pfam, CAZy, and more.
- Show version: `docker run --rm eggnogmapper:2.1.12 emapper.py --version`.
