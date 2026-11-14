"""Check the distributed bytes against MANIFEST.json; not a proof checker."""
from pathlib import Path
import hashlib
import json

root = Path(__file__).resolve().parent
records = json.loads((root / 'MANIFEST.json').read_text())
errors = []
for row in records:
    p = root / row['file']
    if not p.is_file():
        errors.append((row['file'], 'missing'))
        continue
    data = p.read_bytes()
    if len(data) != row['bytes']:
        errors.append((row['file'], 'size differs'))
    if hashlib.sha256(data).hexdigest() != row['sha256']:
        errors.append((row['file'], 'digest differs'))
if errors:
    for item in errors:
        print(*item, sep=': ')
    raise SystemExit(1)
print(f'Integrity verified for {len(records)} files.')
