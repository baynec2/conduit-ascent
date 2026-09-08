#!/usr/bin/env python3
"""subset_mzml.py — extract a scan range from an indexed mzML without altering CV terms.

Why this exists: msconvert in our env rejects newer PSI-MS CV terms emitted by
ThermoRawFileParser (e.g. MS:1003378 "Orbitrap Astral"), so it can't be used to
subset Astral mzMLs. Stripping the CVs to make msconvert happy in turn breaks
DIA-NN. This script bypasses msconvert entirely: it uses pyteomics' on-disk
spectrum offset index to copy bytes for the requested scan range, then writes
a fresh indexList and closing tags.

Usage:
    python subset_mzml.py <input.mzML> <output.mzML> <scan_start> <scan_end>

Notes:
- scan_start/end are inclusive, 1-based scan numbers (as in the spectrum
  idRef "controllerType=0 controllerNumber=1 scan=N").
- Output is valid indexed mzML; fileChecksum is recomputed.
- Chromatogram list is preserved unchanged.
"""

import hashlib
import re
import sys
from pyteomics.mzml import MzML


def _read_chromatogram_list_bytes(path, header_end_offset):
    """Return raw bytes for <chromatogramList>...</chromatogramList> block, or b'' if none."""
    with open(path, "rb") as f:
        f.seek(header_end_offset)
        tail = f.read()  # everything from end of spectrumList onward
    # Find chromatogramList opening tag
    m = re.search(rb"<chromatogramList[^>]*>", tail)
    if not m:
        return b""
    start = m.start()
    end_m = re.search(rb"</chromatogramList>", tail)
    if not end_m:
        return b""
    return tail[start:end_m.end()]


def _read_spectrum_blob(path, off, end):
    with open(path, "rb") as f:
        f.seek(off)
        return f.read(end - off)


def subset_mzml(input_path, output_path, scan_start, scan_end):
    reader = MzML(input_path, use_index=True)
    index = reader._offset_index["spectrum"]

    scan_re = re.compile(r"scan=(\d+)")
    spectrumref_re = re.compile(rb'spectrumRef="[^"]*scan=(\d+)"')
    mslevel_re = re.compile(rb'accession="MS:1000511"\s+value="(\d+)"')

    # Build sorted index of all spectra: (scan_number, idref, byte_offset).
    all_scans = []
    for idref, offset in index.items():
        m = scan_re.search(idref)
        if not m:
            continue
        all_scans.append((int(m.group(1)), idref, offset))
    all_scans.sort()

    # offset_by_scan + next_offset for fast blob access.
    sorted_byte_offsets = [off for _, _, off in all_scans]
    next_offset = {sorted_byte_offsets[i]: sorted_byte_offsets[i + 1]
                   for i in range(len(sorted_byte_offsets) - 1)}
    with open(input_path, "rb") as f:
        f.seek(sorted_byte_offsets[-1])
        chunk = f.read(20 * 1024 * 1024)
        slist_end_in_chunk = chunk.find(b"</spectrumList>")
        if slist_end_in_chunk < 0:
            raise SystemExit("Could not locate </spectrumList> after final spectrum")
        spectrumlist_end_offset = sorted_byte_offsets[-1] + slist_end_in_chunk
    next_offset[sorted_byte_offsets[-1]] = spectrumlist_end_offset

    # Naive in-range set first.
    in_range_scans = {n for n, _, _ in all_scans if scan_start <= n <= scan_end}

    # DIA MS2 spectra reference their precursor MS1 by scan number via
    # spectrumRef="...scan=N". If we drop those precursor scans, DIA-NN
    # dereferences them and segfaults. Pre-extend the keep set with every
    # referenced precursor scan that exists in the source file.
    source_scan_set = {n for n, _, _ in all_scans}
    extended = set(in_range_scans)
    for n, _, off in all_scans:
        if n not in in_range_scans:
            continue
        blob = _read_spectrum_blob(input_path, off, next_offset[off])
        # Skip non-MS2 quickly (only MS2 has precursor refs).
        ms_level_m = mslevel_re.search(blob)
        if not ms_level_m or ms_level_m.group(1) != b"2":
            continue
        for ref_m in spectrumref_re.finditer(blob):
            ref_scan = int(ref_m.group(1))
            if ref_scan in source_scan_set:
                extended.add(ref_scan)

    items = [(n, idref, off) for n, idref, off in all_scans if n in extended]
    if not items:
        raise SystemExit(f"No spectra found in scan range [{scan_start},{scan_end}]")
    items.sort()
    if len(items) != len(in_range_scans):
        added = len(items) - len(in_range_scans)
        print(f"Extended subset to include {added} additional MS1 precursor scans "
              f"that MS2 spectra reference (avoids DIA-NN segfault on missing refs)")

    # The mzML "header" (everything up to and including the <spectrumList ...>
    # opening tag, but excluding any actual <spectrum> elements) runs from byte 0
    # to the offset of the file's very first <spectrum>. Note this is NOT the
    # first KEPT spectrum — using the kept offset would copy every preceding
    # spectrum into the output header verbatim.
    first_in_file_offset = min(v for v in index.values())
    with open(input_path, "rb") as f:
        header = f.read(first_in_file_offset)

    # Update the spectrumList count to match what we're keeping.
    new_count = len(items)
    header = re.sub(
        rb'<spectrumList count="\d+"',
        f'<spectrumList count="{new_count}"'.encode(),
        header,
        count=1,
    )

    # Read selected spectra in order and renumber their `index` attribute.
    # The original `index="N"` was the position in the original spectrumList; if
    # we leave it, DIA-NN allocates `count` slots (the new count) and then tries
    # to write to position N >= count, causing a segfault. Re-anchor each
    # spectrum's index to its 0-based position in the output spectrumList.
    spectra_blobs = []
    spectrum_index_re = re.compile(rb'(<spectrum [^>]*?)index="\d+"')
    with open(input_path, "rb") as f:
        for new_idx, (_, _, off) in enumerate(items):
            end = next_offset[off]
            f.seek(off)
            blob = f.read(end - off)
            blob = spectrum_index_re.sub(
                rb'\1index="' + str(new_idx).encode() + b'"',
                blob,
                count=1,
            )
            spectra_blobs.append(blob)

    # Skip chromatogramList in the subset — it's optional per the mzML spec
    # and DIA-NN doesn't require it. Including a partially-valid one
    # (chromatograms reference scans we may have dropped) is a likely cause
    # of DIA-NN segfaults; safer to omit.
    chrom_bytes = b""
    _ = _read_chromatogram_list_bytes  # silence unused-function warning

    # Assemble the new mzML body (without indexList yet — we'll compute its
    # offset after writing the body, then patch the indexListOffset at the end).
    out_chunks = [header]
    spectra_offsets_in_output = []
    running = len(header)
    for blob in spectra_blobs:
        spectra_offsets_in_output.append(running)
        out_chunks.append(blob)
        running += len(blob)
    out_chunks.append(b"</spectrumList>\n      ")
    running += len(out_chunks[-1])
    chrom_offsets_in_output = []
    if chrom_bytes:
        # Find chromatogram offsets within chrom_bytes for the index.
        # idRef is the `id` attribute of each <chromatogram> element.
        for m in re.finditer(rb'<chromatogram[^>]*id="([^"]+)"[^>]*>', chrom_bytes):
            chrom_offsets_in_output.append((m.group(1).decode(), running + m.start()))
        out_chunks.append(chrom_bytes)
        running += len(chrom_bytes)
    out_chunks.append(b"\n    </run>\n  </mzML>\n  ")
    running += len(out_chunks[-1])

    # Build the new indexList.
    index_chunks = [b'<indexList count="2">\n    <index name="spectrum">\n']
    for (_, idref, _), out_off in zip(items, spectra_offsets_in_output):
        index_chunks.append(
            f'      <offset idRef="{idref}">{out_off}</offset>\n'.encode()
        )
    index_chunks.append(b"    </index>\n    <index name=\"chromatogram\">\n")
    for idref, off in chrom_offsets_in_output:
        index_chunks.append(
            f'      <offset idRef="{idref}">{off}</offset>\n'.encode()
        )
    index_chunks.append(b"    </index>\n  </indexList>\n")
    index_list_blob = b"".join(index_chunks)

    index_list_offset = running
    out_chunks.append(index_list_blob)
    running += len(index_list_blob)

    out_chunks.append(
        f"  <indexListOffset>{index_list_offset}</indexListOffset>\n  ".encode()
    )
    running += len(out_chunks[-1])

    # Compute SHA-1 over everything written so far + the <fileChecksum> opening
    # tag, per the indexed-mzML spec.
    body = b"".join(out_chunks) + b"<fileChecksum>"
    sha1 = hashlib.sha1(body).hexdigest()
    out_chunks.append(f"<fileChecksum>{sha1}</fileChecksum>\n</indexedmzML>\n".encode())

    with open(output_path, "wb") as f:
        f.write(b"".join(out_chunks))

    print(f"Wrote {len(items)} spectra ({sum(len(b) for b in spectra_blobs) / 1e6:.1f} MB of spectra) → {output_path}")
    print(f"Total output size: {sum(len(c) for c in out_chunks) / 1e6:.1f} MB")


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print(__doc__)
        sys.exit(1)
    subset_mzml(sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]))
