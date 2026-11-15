"""Verify the published file manifest without third-party dependencies."""
from pathlib import Path
import json,hashlib
root=Path(__file__).resolve().parent
manifest=json.loads((root/'MANIFEST.json').read_text())
for entry in manifest['files']:
    p=root/entry['file']
    if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest()!=entry['sha256']:
        raise SystemExit('Missing or changed file: '+entry['file'])
print('Verified',len(manifest['files']),'files.')
