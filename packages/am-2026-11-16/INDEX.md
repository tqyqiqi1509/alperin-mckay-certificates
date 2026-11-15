# Calculation index

The sections below follow the mathematical uses in the manuscript. [Recorded results](results/) and [run instructions](README.md#reproduce) accompany the sources.

## Type D

The actual triality stabilizers and their integral root kernels.

[odd_d_round4_check.py](computations/odd_d_round4_check.py).

## E6 roots and carriers

E6 roots and embedded centralizers; cubic carriers; the 27-weight D4 pinning; cubic centralizer-index valuations.

[e6_full_root_geometry_round3.py](computations/e6_full_root_geometry_round3.py), [e6_carrier_round3_checks.py](computations/e6_carrier_round3_checks.py), [e6_D4_carrier_round3_checks.py](computations/e6_D4_carrier_round3_checks.py), [e6_cubic_index_round3_checks.py](computations/e6_cubic_index_round3_checks.py).

## E6 source phases

The Cartan determinant, two-A2 Gram matrix and rational phases. The Python file checks the same exact pairings. Root closure and carrier relations are in the preceding group.

[e67_source_regular_phases.wl](computations/e67_source_regular_phases.wl), [e67_source_regular_phases.py](computations/e67_source_regular_phases.py).

## E7 source and Levi lattices

The 56-dimensional integral carrier, selected source data, and the separate exhaustive 127-proper-Levi table.

[e7_uniform_carrier_11.py](computations/e7_uniform_carrier_11.py), [e7_global20_rootdata.py](computations/e7_global20_rootdata.py), [levi_lattices.py](computations/e7_whole_76/levi_lattices.py).

## E7 ordinary coefficients

The natural carrier generates the GAP input. All 22 subgroup classes, 47 normal pairs and 260 ordinary choices are checked. The .g file is generated and retained for direct inspection.

[independent_carrier.py](computations/e7_56/independent_carrier.py), [independent_coefficients.py](computations/e7_56/independent_coefficients.py), [independent_coefficients.g](computations/e7_56/independent_coefficients.g).

## E8 cubic root action

Integral signed-root action and the pinned cubic A2 carrier.

[e8_tits_engine_root17.py](computations/e8_tits_engine_root17.py), [m69_cubic_pinning.py](computations/m69_cubic_pinning.py).

## Quartic carriers

Run in this order. The lift scan corrects the external A1 lift. The overlap check proves the finite central kernel is exactly the listed two-torsion, not the entire Sylow two-subgroup of the rational center.

[e8_quartic19.py](computations/e8_quartic19.py), [e8_quartic19_phase.py](computations/e8_quartic19_phase.py), [e8_quartic19_liftscan.py](computations/e8_quartic19_liftscan.py), [e8_quartic19_structural.py](computations/e8_quartic19_structural.py), [e8_quartic19_internal_source.py](computations/e8_quartic19_internal_source.py), [e8_quartic19_characters.py](computations/e8_quartic19_characters.py), [e8_quartic19_field.py](computations/e8_quartic19_field.py), [quartic_centre_overlap.py](computations/quartic_centre_overlap.py).

## D4 source and GGGR

Actual D4 graph identification, A3(22) truncated induction and Springer class 5111, exact even grading, trace functional and Smith factors. The Python GGGR checks are an independent exact alternative to the Wolfram expression.

[r4_d4_identity.py](computations/r4_d4_identity.py), [r4_d4_springer.g](computations/r4_d4_springer.g), [m69_gggr_dual_functional.wl](computations/m69_gggr_dual_functional.wl), [cubic_functional_exact.py](computations/cubic_functional_exact.py), [cubic_functional_smith.py](computations/cubic_functional_smith.py).

## Generic source catalogues

GAP3/CHEVIE generators of actual source records, E6/E7 exceptional records, cyclic/zero tails and the 166 E8 unipotent degrees. Recorded exports are in data/chevie.

[e8_strict21_actual_sources.g](computations/e8_strict21_actual_sources.g), [e8_strict21_exceptional_catalog.g](computations/e8_strict21_exceptional_catalog.g), [e8_strict21_cyclic_tails.g](computations/e8_strict21_cyclic_tails.g), [e8_nonprincipal22_generic.g](computations/e8_nonprincipal22_generic.g).

## All-prime valuations

Parse the generic exports and verify exact lattice/valuation conditions, including excluded and non-Sylow source rows.

[e8_strict21_parse_actual.py](computations/e8_strict21_parse_actual.py), [e8_strict21_check_tails.py](computations/e8_strict21_check_tails.py), [e66_full_J_valuation.py](computations/e66_full_J_valuation.py), [e8_nonprincipal22_valuations.py](computations/e8_nonprincipal22_valuations.py).

## E8 principal and nonprincipal lattice checks

Coxeter Smith forms at 3 and 5; full-E8-lattice nonprincipal five-cycle checks; the D4/F4 source geometry and centralizer/color counts.

[e8_principal_three_coxeter_11.py](computations/e8_principal_three_coxeter_11.py), [geometry.py](computations/e8_three_five_75/geometry.py), [e8_d4_source_structure.py](computations/e8_d4_source_structure.py).

## E8 Phi4 Sylow five

Generate GAP, compute its actual order-five element, then rerun Python on the observed output. The final relative group has order 46080. The .g file is generated and retained for direct inspection.

[e8_phi4_sylow5.py](computations/e8_phi4_sylow5.py), [e8_phi4_sylow5.g](computations/e8_phi4_sylow5.g).

## Cover block and quotient tables

Central colors, block rows, quotient fusions, faithful local colors and the seven-local affine model.

[c60_independent_block_and_kernel.g](computations/c60_independent_block_and_kernel.g), [round2_colours.g](computations/round2_colours.g), [round2_smallgrp.g](computations/round2_smallgrp.g), [round2_seven.g](computations/round2_seven.g), [full_outer_blocks.g](computations/exceptional_covers/full_outer_blocks.g), [central_block_quotients.g](computations/exceptional_covers/central_block_quotients.g).

## Actual cover inclusions

Actual inclusions, outer actions, inertia restrictions and projective obstructions for the exceptional covers.

[l34_actual_outer.g](computations/exceptional_covers/l34_actual_outer.g), [l34_all_odd_projective.g](computations/exceptional_covers/l34_all_odd_projective.g), [l34_fourfold_actual.g](computations/exceptional_covers/l34_fourfold_actual.g), [o8_weyl_normalizers.g](computations/exceptional_covers/o8_weyl_normalizers.g), [u6_fi22_actual_normalizer.g](computations/exceptional_covers/u6_fi22_actual_normalizer.g), [e6_three_actual_outer.g](computations/exceptional_covers/e6_three_actual_outer.g).

## Local cover models

Small matrix carriers, inertial-quotient possibilities and the sixfold five-local model. The required helper is actual_c2_pairing.g; Fi22 permutation inputs and their provenance are in exceptional_covers/external_data.

[local_matrix_carriers.g](computations/exceptional_covers/local_matrix_carriers.g), [inertia_module_filters.g](computations/exceptional_covers/inertia_module_filters.g), [sixfold_five_normalizer.g](computations/exceptional_covers/sixfold_five_normalizer.g).
