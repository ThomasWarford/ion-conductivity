"""Short checks for build_runs.py, run before submitting anything.

    conda run -n dft python test_build_runs.py
"""
from __future__ import annotations

import build_runs as br

# (system, temp) -> expected atom count, cross-checked against
# sp_validation/runs/perf_2500K/<system>/POSCAR
EXPECTED_NATOMS = {
    ("Cl-Li", 1000): 100,
    ("Li-P", 2500): 100,
    ("Ge-Li-Sn", 1500): 108,
    ("Li-Nb-S", 2000): 80,
}


def test_nsw_formula():
    # NSW = int(0.85 * 604800 / s_per_step)
    assert br.nsw_for(9.38) == int(0.85 * 604800 / 9.38) == 54805
    assert br.nsw_for(32.11) == int(0.85 * 604800 / 32.11) == 16009
    assert br.nsw_for(29.15) == int(0.85 * 604800 / 29.15) == 17635


def test_manifest_dedup():
    base_cases, extra_cases = br.build_manifest()
    assert len(base_cases) == 26, len(base_cases)
    assert len(extra_cases) == 20, len(extra_cases)
    assert base_cases.isdisjoint(extra_cases)
    # every extra case is one of the two chosen systems
    assert {system for system, _, _ in extra_cases} == set(br.EXTRA_SYSTEMS)


def test_no_random_seed():
    # RANDOM_SEED hangs vasp_gam indefinitely on this cluster's build (confirmed by
    # isolating it as the only changed setting from an otherwise-working INCAR) --
    # must never appear in a generated INCAR, regardless of value.
    incar = br.incar_for(1500, ispin=2, nsw=1000)
    assert "RANDOM_SEED" not in incar


def test_structure_loading():
    for (system, temp), expected in EXPECTED_NATOMS.items():
        matches = sorted(
            (br.DATA / f"temperature={temp}K" / f"chemical_system={system}").glob(
                "*_last_quarter.traj"
            )
        )
        assert len(matches) == 1, (system, temp, matches)
        struct = br.load_final_structure(system, temp)
        assert len(struct) == expected, (system, temp, len(struct))


def test_incar_hygiene():
    incar = br.incar_for(1500, ispin=2, nsw=1000)
    assert incar["LAECHG"] is False
    assert incar["LCHARG"] is False
    assert incar["LORBIT"] is None
    assert incar["TEBEG"] == incar["TEEND"] == 1500
    assert incar["ISPIN"] == 2
    assert incar["NSW"] == 1000


def main() -> None:
    tests = [
        test_nsw_formula,
        test_manifest_dedup,
        test_no_random_seed,
        test_structure_loading,
        test_incar_hygiene,
    ]
    for test in tests:
        test()
        print(f"PASS {test.__name__}")
    print(f"\n{len(tests)}/{len(tests)} passed")


if __name__ == "__main__":
    main()
