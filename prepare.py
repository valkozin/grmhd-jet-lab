#!/usr/bin/env python3
"""Restore bundled simulation data using only the Python standard library."""
from pathlib import Path
import hashlib
import zipfile

ROOT = Path(__file__).resolve().parent
DIST = ROOT / 'dist'


def checked(path, expected):
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != expected:
        raise RuntimeError(f'Checksum mismatch: {path.name}')


def main():
    archive = DIST / 'grmhd-repro.zip'
    expected = '120b47c7040899d4362860fb2ab11f3be14558615662dc54e46c7291903bfe68'
    if not archive.exists():
        temporary = archive.with_suffix('.tmp')
        try:
            with temporary.open('wb') as out:
                for i in range(1, 4):
                    out.write((DIST / f'grmhd-repro.zip.part{i:03}').read_bytes())
            checked(temporary, expected)
            temporary.replace(archive)
        finally:
            temporary.unlink(missing_ok=True)
    checked(archive, expected)
    data = DIST / 'visualization-data.zip'
    checked(data, '3c9c023633c3ed0456f120783a9b5bf0c716efaec351e82d61d829cd95d1ccb6')
    with zipfile.ZipFile(data) as z:
        for info in z.infolist():
            destination = (DIST / info.filename).resolve()
            if not destination.is_relative_to(DIST.resolve()):
                raise RuntimeError('Unsafe archive path')
        z.extractall(DIST)
    print('Data restored and SHA-256 verified.')
    print('Run: python3 -m http.server 8000 --directory dist')


if __name__ == '__main__':
    main()
