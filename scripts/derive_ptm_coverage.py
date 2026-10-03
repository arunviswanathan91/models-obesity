"""Export eight coverage-plot counts from the original completed PTM report ZIP.

Usage: python derive_ptm_coverage.py Protein_PTM_reports.zip output.csv
Uses the exact filters in the original fig_coverage; no model fitting or NC input.
"""
import csv
import hashlib
import io
import json
import sys
import zipfile
from collections import Counter
from pathlib import Path


def derive(archive, output):
    archive, output = Path(archive), Path(output)
    with zipfile.ZipFile(archive) as z:
        members = {name: z.read(name) for name in (
            'coverage_exclusions.csv', 'eligible_features.csv', 'fit_status.csv')}
    tables = {name: list(csv.DictReader(io.StringIO(raw.decode('utf-8-sig'))))
              for name, raw in members.items()}
    coverage = tables['coverage_exclusions.csv']
    stages = [
        ('source_features', coverage),
        ('signature_catalogue_matched', [r for r in coverage if
            r['exclusion_reason'] != 'outside original genes and named targets']),
        ('paired_coverage_eligible', [r for r in tables['eligible_features.csv']
            if r['included'].lower() == 'true']),
        ('primary_numerical_checks_pass', [r for r in tables['fit_status.csv']
            if r['variant'] == 'primary' and r['numerical_checks_pass'].lower() == 'true']),
    ]
    with output.open('w', newline='') as h:
        writer = csv.DictWriter(h, fieldnames=['stage_order', 'stage', 'assay', 'count'])
        writer.writeheader()
        for index, (stage, records) in enumerate(stages, 1):
            counts = Counter(r['assay'] for r in records)
            for assay in ('phospho', 'glyco'):
                writer.writerow(dict(stage_order=index, stage=stage, assay=assay, count=counts[assay]))
    provenance = {'archive_sha256': hashlib.sha256(archive.read_bytes()).hexdigest(),
                  'member_sha256': {n: hashlib.sha256(b).hexdigest() for n, b in members.items()},
                  'filters': 'Original fig_coverage stage filters; counts only',
                  'derived_csv_sha256': hashlib.sha256(output.read_bytes()).hexdigest()}
    output.with_suffix('.provenance.json').write_text(json.dumps(provenance, indent=2))


if __name__ == '__main__':
    derive(*sys.argv[1:])
