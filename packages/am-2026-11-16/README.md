# Computations for The Alperin–McKay conjecture

Qiyue Tang — 16 November 2026.

This package contains the finite calculations used in the accompanying manuscript: source programs, their necessary inputs, exact results, and a reproducible runner. The [calculation index](INDEX.md) identifies the mathematical use of each group of files.

## Contents

- `computations/`: one version of each required program and its local dependencies.
- `data/chevie/`: four recorded generic-character exports, with their generating GAP3 programs included.
- `results/`: exact JSON results, finite-group outputs, and execution receipts from the packaged programs.
- `SOURCE_VERSIONS.json`: original and packaged source hashes; line-ending normalization is recorded.
- `MANIFEST.json` and `verify.py`: integrity check for every distributed file.

Intermediate pickles, Python caches, research transcripts, manuscript PDFs, superseded variants, and the GAP/CHEVIE installations are not distributed. The runner regenerates intermediate data in its output directory.

## Requirements

The tested Python environment is Python 3.13, SymPy 1.14.0 and NumPy 2.5.1. Install the Python dependencies with:

```sh
python -m pip install -r requirements.txt
python verify.py
```

Group computations use GAP4 4.16.1 with CTblLib 1.3.11, `perfgrp`, and `smallgrp`. Generic-character computations use [GAP3/CHEVIE](https://github.com/jmichel7/gap3-jm), revision `64365f30757f5374c86511b377d7d98eb15c1e02`. GAP3 and GAP4 are different runtimes; the CHEVIE generators here are GAP3 programs.

## Reproduce

Run commands from this package directory. Each invocation creates a new isolated output tree and preserves the distributed results.

```sh
# Python checks, using the supplied generic-character exports:
python run.py --suite python

# Regenerate generic-character exports and check their valuations:
python run.py --suite gap3 --gap3-command '["/path/to/gap3-jm/bin/gap.linux","-m","512m","-l","/path/to/gap3-jm/lib/","-h","/path/to/gap3-jm/doc/","-b"]'

# GAP4 checks, including both mixed Python/GAP pipelines:
python run.py --suite gap4 --gap4-command '["gap","-q","-b","--quitonbreak"]'
```

Use `--suite all` with both runtime options to run the entire sequence. Runtime options are JSON argument arrays, so no shell interpolation is needed. `--output PATH` selects a new results directory. The runner checks exit status, completion markers where supplied, and the final E7 coefficient counts and E8 Phi4 Sylow result. Do not use Python's `-O` option.

The two `.wl` files are standalone Wolfram expressions. They are retained because the manuscript names them. Their phase and GGGR assertions are also reproduced by the supplied exact Python programs; a Wolfram installation is not needed for the default runner. With WolframScript, evaluate a file using `wolframscript -code 'Print[Get["computations/e67_source_regular_phases.wl"]]'` (replace the filename for the GGGR expression).

## Reading the results

The uncorrected external-A1 lift in the quartic pipeline fails its prescribed line condition. This is an intermediate calculation: the lift-scan program finds the torus correction, and the character and field checks use that corrected lift. Its failed intermediate row is preserved in the result data.

The E8 Phi4 calculation is complete only after Python generates the GAP program, GAP supplies the actual order-five element, and Python checks it on the integral lattice. The runner enforces this order and requires the final relative-group order 46080. Likewise, the E7 coefficient program is regenerated from its actual carrier before GAP checks 22 subgroup classes, 47 pairs and 260 ordinary-character choices.

The finite calculations verify their stated matrix, lattice, character and group assertions. Their mathematical interpretation and the deduction of the manuscript's theorems are given in the paper.

Atlas permutation inputs retain their [source provenance](computations/exceptional_covers/external_data/atlas_generators_provenance.json). Published character tables and CHEVIE data remain attributed to the sources cited in the manuscript.
