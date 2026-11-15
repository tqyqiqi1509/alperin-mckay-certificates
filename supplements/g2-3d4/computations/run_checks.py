#!/usr/bin/env python3
"""Execute the four exact checks and compare their complete JSON results."""
import json
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parent
CHECKS = (
    "verify_g2_torus",
    "verify_three_module",
    "verify_triality_roots",
    "verify_triality_chevalley",
)


def main():
    output = ROOT / "rerun_results"
    output.mkdir(exist_ok=True)
    for name in CHECKS:
        result = subprocess.run(
            [sys.executable, str(ROOT / "checks" / (name + ".py"))],
            text=True,
            capture_output=True,
            check=True,
        )
        actual = json.loads(result.stdout)
        expected = json.loads(
            (ROOT / "checks" / (name + ".json")).read_text(encoding="utf-8")
        )
        (output / (name + ".json")).write_text(result.stdout, encoding="utf-8")
        if actual != expected:
            raise SystemExit(name + ": output differs from the recorded result")
        print(name + ": passed; complete JSON matches the recorded result")
    print("All four checks passed. See README.txt for their mathematical scope.")


if __name__ == "__main__":
    main()
