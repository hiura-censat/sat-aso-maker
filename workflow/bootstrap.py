#!/usr/bin/env python3
"""Create an immutable, separate workspace for one Snakemake run."""
import argparse
import hashlib
import json
import os
import re
import shutil
import sys
from pathlib import Path

import yaml


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', required=True)
    parser.add_argument('--run', required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    config_path = (root / args.config).resolve()
    config = yaml.safe_load(config_path.read_text())
    sys.path.insert(0, str(root / 'scripts'))
    from workflow_settings import settings
    os.environ['SAT_ASO_CONFIG'] = str(config_path)
    settings()
    run_id = config['run_id']
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', run_id):
        raise ValueError('run_id must contain only letters, digits, _, . and -')
    output = (root / args.run).absolute()
    if output != root / 'runs' / run_id:
        raise ValueError('run directory must equal runs/{run_id}')
    source = (root / config['source_root']).resolve()
    if not (source / 'results/data').is_dir():
        raise FileNotFoundError(source / 'results/data')
    for name in ('python', 'samtools'):
        if not Path(config[name]).is_file():
            raise FileNotFoundError(config[name])
    local = root / '.local'
    if not (local / 'bin/jellyfish').is_file() or not (local / 'bin/jellyfish-prefix').is_file():
        raise FileNotFoundError('local Jellyfish binaries are required')
    scripts = sorted(p for p in (root / 'scripts').iterdir() if p.is_file())
    snapshot = {
        'config_sha256': digest(config_path),
        'source_root': str(source),
        'scripts_sha256': {p.name: digest(p) for p in scripts},
    }
    old = output / '.bootstrap.json'
    if old.exists():
        if json.loads(old.read_text()) != snapshot:
            raise RuntimeError('run_id already exists with different config or scripts; choose a new run_id')
        return
    if output.exists() and any(output.iterdir()):
        raise RuntimeError('run directory already contains files; choose a new run_id')
    output.mkdir(parents=True, exist_ok=True)
    source_link = output.parent / 'VallePrep_v0.0.0'
    if source_link.is_symlink():
        if source_link.resolve() != source:
            raise RuntimeError('runs/VallePrep_v0.0.0 points at a different source')
    elif source_link.exists():
        raise RuntimeError('runs/VallePrep_v0.0.0 already exists and is not a symlink')
    else:
        source_link.symlink_to(source, target_is_directory=True)
    (output / '.local').symlink_to(local, target_is_directory=True)
    (output / 'scripts').mkdir()
    for path in scripts:
        shutil.copy2(path, output / 'scripts' / path.name)
    shutil.copy2(config_path, output / 'workflow_config.yaml')
    old.write_text(json.dumps(snapshot, indent=2) + '\n')


if __name__ == '__main__':
    main()
