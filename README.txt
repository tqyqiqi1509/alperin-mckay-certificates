FINITE CALCULATIONS ACCOMPANYING
A proof of the Alperin--McKay conjecture
Qiyue Tang -- 14 November 2026

Purpose
=======
The article contains the mathematical arguments and the deductions from
the finite calculations. This directory supplies the exact procedures,
their finite inputs, and recorded outputs. It is not a formal proof
assistant development. A recorded success string or an output count is
not a substitute for the mathematical identification made in the article.

Contents and arithmetic
=======================
Python calculations use Python 3, SymPy and, in a few root calculations,
NumPy with integer matrices. Rational, integer and cyclotomic identities
are exact. Ordinary finite-group calculations use GAP 4.16.1 with CTblLib
1.3.11 and, for the small perfect-group models, the Perfect Groups
library (perfgrp). The subgroup identifiers used in some optional
diagnostics also require the Small Groups library.

Files beginning with RequirePackage("chevie") are GAP 3/CHEVIE programs.
They must be run with the CHEVIE implementation associated with GAP 3,
not by passing them unmodified to GAP 4. The generic-character data used
by the Python calculations are supplied as text exports. Their source
programs are supplied alongside them. The pinned export used for the
E6/E8 degree valuation calculation is from jmichel7/gap3-jm commit
64365f30757f5374c86511b377d7d98eb15c1e02. The article cites CHEVIE and
the primary mathematical sources for these generic data.

Run commands from this supplement directory. Python files generally
locate their inputs relative to __file__. The exceptional-cover GAP
programs use the relative paths given here. Several calculations write
their output beside their source. Preserve an untouched copy if the
original recorded outputs are needed for byte comparison.

Main calculation index
======================
Classical groups and genuine triality:
  computations/d4_trivial_triality_root_08.py
  computations/classical_triality_lattice.py

F4, uniform mixed-torus and source structure:
  computations/f4_mixed_torus_uniform_root_08.py
  computations/f4_structural_72/construct.py

Exceptional rational Levi lattices and actual field returns:
  computations/exceptional_69/levi_lattices.py
  computations/exceptional_69/e6_t2_returns.py

E6:
  computations/e6_root_independent_40.py
  computations/e6_a5a1_e3_40.py
  computations/e6_a5a1_quadratic_40.py
  computations/e6_a5a1_bad3_43.py
  computations/e6_a5a1_hook3_43.py
  computations/e6_complex_parabolics_37.py and .g
  computations/e6_nonmatching_three_37.py and .g
  computations/e6_independent_and_order3_46.py and .g
  computations/e6_principal_three_root_08.py
  computations/e6_d4_independent_32.py
  computations/e6_matching54/
  computations/e6_nonisolated_65/
  computations/e6_whole_74/

E7:
  computations/e7_diagonal_and_toral_41.py
  computations/e7_euclidean_carrier_21.py
  computations/e7_good_principal_toral_13.py
  computations/e7_order_three_six/
  computations/e7_d6_carrier/
  computations/e7_56/
  computations/e7_whole_76/

E8 isolated sources of semisimple orders five and six:
  computations/e8_55/independent_geometry.py
  computations/e8_55/gluing_and_cores.py
  computations/e8_55/independent_tits_carrier.py
  computations/e8_55/source_manifest.py
  computations/e8_55/normalizers.g

E8 order-three sources:
  computations/e8_60/geometry.py
  computations/e8_60/gluing_and_defect.py
  computations/e8_60/carrier.py
  computations/e8_60/normalizers.g
  computations/e8_three_five_75/

E8 involution, order-four, and bad-prime sources:
  computations/e8_58/independent_geometry.py
  computations/e8_58/centralizer_details.py
  computations/e8_58/independent_pinning.py
  computations/e8_58/even_gggr.py
  computations/e8_bad_68/construct.py
  computations/e8_classical_sources_independent_27.py
  computations/e8_e6theta_source_27.py
  computations/e8_e7a1_e4_30.py
  computations/e8_d4_source_structure.py

E8 principal and unipotent degree calculations:
  computations/e8_principal_three_coxeter_11.py
  computations/e8_phi4_sylow5.py and .g
  computations/e8_unipotent_degree_valuation.py
  computations/e8_unipotent_degrees_reexport.g
  computations/e8_five_series_14.g

Exceptional Schur covers, actual inclusions, outer actions and blocks:
  computations/exceptional_covers/central_block_quotients.g
  computations/exceptional_covers/full_outer_blocks.g
  computations/exceptional_covers/l34_actual_outer.g
  computations/exceptional_covers/l34_all_odd_projective.g
  computations/exceptional_covers/l34_fourfold_actual.g
  computations/exceptional_covers/o8_weyl_normalizers.g
  computations/exceptional_covers/u6_fi22_actual_normalizer.g
  computations/exceptional_covers/e6_three_actual_outer.g
  computations/exceptional_covers/sixfold_five_normalizer.g
  computations/exceptional_covers/local_matrix_carriers.g
  computations/exceptional_covers/inertia_module_filters.g

Dependencies and examples
=========================
The E8 directories 55, 58 and 60 have their geometric JSON input included.
To regenerate it, run the geometry program before its carrier, gluing,
centralizer, pinning or GGGR program. Run source_manifest.py in e8_55
after its geometry, gluing and carrier outputs exist. In e7_56, the
coefficient program imports the independent carrier and writes the GAP
program that checks the ordinary coefficients. In e6_matching54,
active_actions.py reads both the geometry JSON and the supplied exact
GAP output. Regenerate the latter with independent_geometry.g if needed.

Examples:
  python computations/f4_mixed_torus_uniform_root_08.py
  python computations/exceptional_69/e6_t2_returns.py
  python computations/e8_55/source_manifest.py
  python computations/e8_unipotent_degree_valuation.py
  gap -q computations/exceptional_covers/sixfold_five_normalizer.g
  gap3 -q computations/e8_unipotent_degrees_reexport.g

Use a GAP installation with the named packages installed. In a deliberately
minimal GAP installation the option --bare -A may be used before explicit
LoadPackage calls. A startup or package-loading error is a failed run,
even if the surrounding shell happens to return zero. Inspect the actual
calculation output and its assertions. A GAP break loop is not success.

External finite input
=====================
The two permutation generators for the official 3510-point Fi22:2 action
are included in computations/exceptional_covers/external_data/. The
adjacent atlas_generators_provenance.json records their original ATLAS
URLs and SHA-256 hashes. The other displayed affine, root-permutation,
and matrix models are constructed by the supplied programs.

Packaging changes
=================
The exceptional-cover Read paths are relative. Two E7 calculation
directories have descriptive names under computations/. Mathematical
inputs and assertions have not changed. The E8 degree-valuation program
omits a redundant comparison with an older export; it still checks the
complete supplied 30-row E6 and 166-row E8 data and the Ennola valuations.
The E6 principal-three program uses its explicit Weyl word as its full
input. Its old upstream file hash is labelled as recorded provenance;
the calculation no longer tries to open an unrelated historical checkout.
Both adjusted Python programs were executed after these packaging edits.
The triality program's output name follows its packaged filename; that
path-only adjustment was also executed and checked.

Reproduction and review status
=============================
The package contains both programs and recorded outputs. It does not claim
that every supplied export was generated again during the final manuscript
review. The separate audit report identifies newly executed checks,
comparisons with recorded outputs, and literature inputs. In particular,
generic CHEVIE exports and published Small Groups identifications are
distinguished from the independently executed actual-group calculations.

Integrity
=========
MANIFEST.json records each distributed file's size and SHA-256 hash, and
also retains the earlier source hash where available. Run

  python verify_manifest.py

before modifying or regenerating the outputs. This checks file integrity,
not the mathematical validity of a certificate.
