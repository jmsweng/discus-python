"""Tests for the pre-registered Gate 2 experiment tooling (synthetic only; no engine data).

These check that the registered inputs are reproducible, that the oracle's fast lattice-sum path
equals a brute-force atom sum, that the package transform reproduces the independent real-space
prediction, and that analyze.py passes correct data and rejects corrupted data. They do NOT show
that DISCUS passes Gate 2; only the real engine run analysed by analyze.py can do that.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("h5py")
G2 = Path(__file__).resolve().parents[1] / "benchmarks" / "discus_6_23_01" / "gate2_prereg"
sys.path.insert(0, str(G2))
import oracle as O  # noqa: E402

PRED = json.loads((G2 / "predictions.json").read_text())


def test_registered_files_unchanged_and_reproducible():
    for name, h in PRED["code_sha256"].items():
        assert O.sha256_file(G2 / name) == h, f"{name} changed after registration"
    for name, h in PRED["input_sha256"].items():
        assert O.sha256_file(G2 / "inputs" / name) == h
    for model in ("block", "random"):
        sig = O.build_sigma(model)
        text = O.stru_text(sig, f"gate2_{model}_seed{O.SEEDS[model]}")
        assert (G2 / "inputs" / f"gate2_{model}.stru").read_text() == text
        assert np.array_equal(O.parse_stru(G2 / "inputs" / f"gate2_{model}.stru"), sig)


def test_atom_lines_parse_with_discus_fixed_column_reader():
    """DISCUS reads 'format numbers' atom lines with a fixed-column Fortran format (atom_line.f90).
    gate2-prereg-v1 failed exactly this on the real engine (error -49)."""
    for model in ("block", "random"):
        sig = O.build_sigma(model)
        names, pos = O.atom_list(sig)
        lines = [l for l in (G2 / "inputs" / f"gate2_{model}.stru").read_text().splitlines()
                 if l[:2] in ("CU", "AU")]
        assert len(lines) == 4096
        for line, nm, p in zip(lines, names, pos):
            got = O.discus_numbers_read(line)
            assert got == (nm, p[0], p[1], p[2], 0.0)
    # positive control: the archived, DISCUS-written C1 file
    c1 = [l for l in (G2.parent / "gate1_c1_dyadic" / "inputs" / "C1_frozen.stru").read_text().splitlines()
          if l[:2] in ("GA", "AS")]
    assert all(O.discus_numbers_read(l) for l in c1)
    # negative control: the v1 layout (two blanks after the name) must be rejected
    # literal first atom line of the v1 file (commit c3ab156), which DISCUS rejected
    v1 = ("CU          0.000000,       0.000000,       0.000000,      0.000000,"
          "       1,       0,       0,   1.000000, _,   0,   0,   0")
    with pytest.raises(ValueError):
        O.discus_numbers_read(v1)


def test_structure_is_exactly_stoichiometric_and_ordered_for_celltoindex():
    for model in ("block", "random"):
        sig = O.build_sigma(model)
        assert sig.sum() == 0 and sig.shape == (16, 16, 16)
        names, pos = O.atom_list(sig)
        # DISCUS celltoindex: iatom = ((iz*ny + iy)*nx + ix)*ncatoms + isite, x fastest
        i = np.arange(len(pos))
        assert np.array_equal(pos[:, 0], i % 16) and np.array_equal(pos[:, 1], (i // 16) % 16)
        assert np.array_equal(pos[:, 2], i // 256)


@pytest.mark.parametrize("mode", ["physics", "engine"])
def test_lattice_sum_forward_model_equals_brute_force(mode):
    sig = O.build_sigma("block")
    It, Id, g2 = O.forward(sig, mode)
    rng = np.random.default_rng(1)
    idx = rng.integers(0, O.N_H, size=(200, 3))
    hkl = O.H_AXIS[idx]
    F = O.direct_F(sig, hkl, mode)
    np.testing.assert_allclose(np.abs(F) ** 2, It[tuple(idx.T)], rtol=1e-10, atol=1e-6)
    # diffuse = g^2 S exactly for this exactly stoichiometric, one-site-per-cell crystal
    S = np.abs(np.fft.ifftn(sig.astype(float)) * sig.size) ** 2
    k = np.mod(np.rint(O.H_AXIS * 16).astype(int), 16)
    np.testing.assert_allclose(Id, g2 * S[np.ix_(k, k, k)], rtol=1e-9, atol=1e-6)


def test_pair_sums_direct_match_fft_and_analytic_within_block_values():
    sig = O.build_sigma("block")
    C = O.pair_sums_direct(sig)
    Cf = np.rint(np.real(np.fft.ifftn(np.abs(np.fft.fftn(sig)) ** 2))).astype(np.int64)
    assert np.array_equal(C, Cf)
    assert C[0, 0, 0] == sig.size


def test_pipeline_transform_matches_independent_prediction():
    for model in ("block", "random"):
        m = PRED["models"][model]
        for mode in ("physics", "engine"):
            assert m[mode]["pipeline_on_oracle_vs_independent_rho_max_abs_over_rho0"] < 1e-12
            assert m[mode]["pipeline_on_oracle_imag_over_real"] < 1e-12


def test_analysis_passes_clean_and_rejects_corruptions(tmp_path):
    import selftest
    report = selftest.run_selftest(tmp_path)
    bad = [r for r in report if not r[3]]
    assert not bad, bad
