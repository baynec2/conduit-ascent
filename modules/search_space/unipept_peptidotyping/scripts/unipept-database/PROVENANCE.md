# Provenance — vendored `unipept-database`

This directory is a **partial, unmodified copy** of the Unipept project's
database-build tooling, redistributed here under its MIT license.

| | |
|---|---|
| Upstream | <https://github.com/unipept/unipept-database> |
| Commit | [`7f4e5951e25c`](https://github.com/unipept/unipept-database/commit/7f4e5951e25c) (2 July 2025) |
| License | MIT — see [`LICENSE`](LICENSE), copyright (c) 2023 Universiteit Gent |
| Files | 39 — 38 script/README files plus upstream's `LICENSE`, all byte-identical to upstream |
| Modifications | None |

## Why it is vendored rather than fetched

`build_sequence_index` runs `scripts/generate_umgap_tables.sh` to turn a UniProtKB
release into the `sequences.tsv.lz4` / `taxons.tsv.lz4` pair that every peptidotyping
rule reads (see `modules/search_space/_shared/unipept_resources.smk`). Pinning the
tooling in-tree makes the index build reproducible against a fixed version of the
scripts rather than against whatever upstream `master` happens to hold on the day a
user builds their index.

## Scope of the copy

Only the `scripts/` tree and `README.md` were taken. Upstream's CI configuration,
issue templates and sample data were not copied. `scripts/rust-utils/target/` in this
repository is local build output, not upstream content.

## Verifying this claim

Every file listed above was compared against upstream by git blob hash. To re-check:

```bash
BASE=modules/search_space/unipept_peptidotyping/scripts/unipept-database
curl -s "https://api.github.com/repos/unipept/unipept-database/git/trees/7f4e5951e25c?recursive=1" \
  | python3 -c "
import sys, json, subprocess, os
base = os.environ['BASE']
upstream = {t['path']: t['sha'] for t in json.load(sys.stdin)['tree'] if t['type'] == 'blob'}
tracked = subprocess.run(['git', 'ls-files', base], capture_output=True, text=True).stdout.split()
checked = 0
for path in tracked:
    rel = os.path.relpath(path, base)
    # target/ is local build output; PROVENANCE.md is this file, not upstream's.
    if '/target/' in path or rel == 'PROVENANCE.md':
        continue
    checked += 1
    local = subprocess.run(['git', 'hash-object', path], capture_output=True, text=True).stdout.strip()
    if upstream.get(rel) != local:
        print('DIFFERS:', rel)
print('checked', checked, 'files against upstream')
"
```

If a file is ever modified locally, record the change here — MIT permits modification,
but the "unmodified" claim above must stop being made.

## Attribution

Unipept is developed at Universiteit Gent. If you use the peptidotyping search space,
please cite the Unipept publications listed at <https://unipept.ugent.be/publications>.
