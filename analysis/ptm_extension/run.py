"""Run the final PTM sensitivity and calibration analyses from frozen settings."""
from pathlib import Path
import argparse
import hashlib
import json
import os


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage', choices=['sensitivity', 'calibration', 'both'], default='both')
    parser.add_argument('--input-zip', type=Path, required=True)
    parser.add_argument('--signature-zip', type=Path, required=True)
    parser.add_argument('--primary-report', type=Path, required=True,
                        help='Directory containing the frozen all_effects_95HDI.csv.')
    parser.add_argument('--output', type=Path, required=True,
                        help='Separate output directory for this curated source version.')
    parser.add_argument('--work', type=Path, required=True)
    args = parser.parse_args()

    for name in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS']:
        os.environ[name] = '1'
    os.environ['PYTENSOR_FLAGS'] = 'blas__ldflags='

    import run_common as common
    import run_ptm

    here = Path(__file__).resolve().parent
    cfg = json.loads((here / 'settings.json').read_text())
    versions = dict(line.split('==') for line in (here / 'requirements.txt').read_text().splitlines()
                    if '==' in line)
    if common.versions() != versions:
        raise RuntimeError(f'Package versions differ from the fitted analysis: {common.versions()}')

    source = {name: (here / name).read_text() for name in
              ['ptm_core.py', 'run_common.py', 'run_ptm.py', 'run.py']}
    cfg.update(source_hash=hashlib.sha256(json.dumps(source, sort_keys=True).encode()).hexdigest(),
               input_zip=str(args.input_zip.resolve()),
               signature_zip=str(args.signature_zip.resolve()),
               report=str(args.primary_report.resolve()),
               output=str(args.output.resolve()), work=str(args.work.resolve()))
    for key, expected in [('input_zip', 'input_sha256'), ('signature_zip', 'signature_sha256')]:
        if common.sha(Path(cfg[key])) != cfg[expected]:
            raise ValueError(f'Frozen input checksum mismatch: {key}')
    if common.sha(Path(cfg['report']) / 'all_effects_95HDI.csv') != cfg['effects_sha256']:
        raise ValueError('Frozen primary effects checksum mismatch')

    output = Path(cfg['output'])
    output.mkdir(parents=True, exist_ok=True)
    configuration = output / 'configuration.json'
    if configuration.exists() and json.loads(configuration.read_text()) != cfg:
        raise ValueError('Output configuration differs; select a separate output directory.')
    common.write_json(configuration, cfg)
    data, features, arrays = run_ptm.load(cfg)
    if args.stage in ['sensitivity', 'both']:
        run_ptm.sensitivity(cfg, data, features, arrays)
    if args.stage in ['calibration', 'both']:
        run_ptm.calibration(cfg, data, features, arrays)


if __name__ == '__main__':
    main()
