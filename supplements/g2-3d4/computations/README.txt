COMPUTATIONAL SUPPLEMENT
The inductive Alperin--McKay condition for G2(q) and 3D4(q)
14 November 2026

RUNNING THE CHECKS

Use Python 3.9 or later with SymPy 1.14.0 installed. From this directory:

    python3 run_checks.py

The runner executes four programs, checks their assertions, parses their JSON
output and compares the complete result with the recorded JSON. A failed
assertion, failed program or changed result gives a nonzero exit status.
Fresh results are written to rerun_results/. The three programs other than
verify_triality_roots.py use only the Python standard library.

CONTENTS AND SCOPE

1. checks/verify_g2_torus.py
   Constructs the twelve integral matrices for the G2 torus action, checks
   their relations and root evaluations, and checks reduction modulo 3 and 5.
   The manuscript supplies the uniform argument for every odd prime.

2. checks/verify_three_module.py
   Checks the triality action on C_(3^(a+1)) x C_(3^a) for a=1,2,3,4;
   computes fixed points, commutator subgroups and the relevant character
   labels; and checks recorded field-return and unipotent-degree examples.

3. checks/verify_triality_roots.py
   Uses SymPy for an independent D4 root-lattice and Smith-normal-form
   calculation, including the fixed torus and recorded integer parameters.

4. checks/verify_triality_chevalley.py
   Starts with integral eight-dimensional Chevalley root matrices. It
   reconstructs u, v and n0, compares their signed basis images with all
   three columns printed in Appendix A, checks the defining relations and
   the complement order (also after reduction modulo 2), derives the
   coroot action, checks all 24 roots, and checks 15 specified prime powers,
   a=1,2,3,4 and 56 field-return examples. Positive and inverse character
   conventions are checked for the final return congruence. The q=4
   nonprincipal-block example has three fixed labels and three two-cycles
   under F_2^3, with the inverse-p character convention used in the paper.

Each .py file has a corresponding .json containing output actually produced
by that program. environment.json records the software versions used for
this delivery. The JSON outputs are deterministic and compared in full.

RELATION TO THE PROOF

The integer matrix relations hold after reduction in every characteristic.
Appendix A proves their lift to the simply connected group using the action
of triality on its centre. The manuscript proves the torus formulas,
3-adic congruences and root restrictions uniformly. The finite module and
prime-power checks provide separate reproducibility checks of those formulas.

The character parametrizations, block classification and inductive criterion
are the cited mathematical inputs listed in the article. The exceptional
covering groups G2(3) and G2(4) use Thomas Breuer's published verification:
https://www.math.rwth-aachen.de/~Thomas.Breuer/ctblocks/doc/overview.html
This supplement reproduces the arithmetic programs developed with the
manuscript; Breuer's GAP calculations are a separately cited external result.

PROVENANCE

The first three programs and their recorded JSON were recovered unchanged
from Source_and_Verification_20261114.zip and rerun successfully. The fourth
program was independently developed during the current arithmetic review,
then augmented with an explicit comparison to the printed signed table and
the inverse-p convention. provenance.json gives the original and delivered
checksums. The original uploaded manuscript and the historical ZIP are
preserved in the complete source-and-audit archive.
