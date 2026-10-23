"""Tests for the DISCUS 6.23.01 common/nexus adapter.

Two kinds of evidence, kept apart on purpose:

* REAL: ``tests/fixtures/discus_6.23.01/C3_total.h5`` / ``C4_total.h5`` were written by a real
  DISCUS_SUITE 6.23.01 run (macros ``C3.mac``/``C4.mac`` beside them; monoclinic P1 two-atom cell,
  a,b,c = 5,6,7 Å, beta = 105 deg, 8x8x8 cells, unequal axis sizes 33/41/49).
* SYNTHETIC: files written here with the same dataset names to exercise validation and
  to check the transform against the existing Cartesian ``delta_pdf_from_diffuse``.
"""
from pathlib import Path

import h5py
import numpy as np
import pytest

from discuspdf.discus_output import (delta_pdf_fractional, evaluate_real_space, read_discus_volume)
from discuspdf.volume import delta_pdf_from_diffuse

FIX = Path(__file__).parent / "fixtures" / "discus_6.23.01"


def write_discus_like(path, values, *, corner, inc, cell=(5.0, 5.0, 5.0), angles=(90., 90., 90.),
                      bragg="bragg_subtracted", axes="hkl", unit="basecell_fractional",
                      content="intensity", reciprocal="reciprocal", data_axes=(1, 2, 3)):
    with h5py.File(path, "w") as f:
        g = f.create_group("entry/data")
        s = lambda k, v: g.create_dataset(k, data=np.bytes_(v))
        s("audit_conform_dict_name", "Disorder unified data"); s("audit_creation_method", "synthetic-test")
        s("coordinate_unit", unit); s("data_type_axes", axes); s("data_type_reciprocal", reciprocal)
        s("data_type_content", content); s("data_type_with_bragg", bragg)
        g.create_dataset("data_values", data=values)
        g.create_dataset("data_dimension", data=np.array(values.shape, dtype="int32"))
        g.create_dataset("data_axes", data=np.array(data_axes, dtype="int32"))
        g.create_dataset("data_corner", data=np.array(corner, float))
        g.create_dataset("data_increment_vector", data=np.array(inc, float))
        g.create_dataset("unit_cell_lengths", data=np.array(cell, float))
        g.create_dataset("unit_cell_angles", data=np.array(angles, float))


def centered_blob(n=33, d=0.125):
    q = (np.arange(n) - n // 2) * d
    x, y, z = np.meshgrid(q, q, q, indexing="ij")
    return np.exp(-((x - 1)**2 + y**2 + z**2) / 0.3) + np.exp(-((x + 1)**2 + y**2 + z**2) / 0.3) + 0.01


# ------------------------------ REAL DISCUS files ------------------------------
def test_real_monoclinic_volume_metadata_and_axis_order():
    v = read_discus_volume(FIX / "C3_total.h5")
    assert v.values.shape == (33, 41, 49)          # (na, no, nt) as set in C3.mac
    assert np.allclose(v.corner, (-4, -5, -6)) and np.allclose(np.diag(v.increments), 0.25)
    assert np.allclose(v.cell_lengths, (5, 6, 7)) and np.allclose(v.cell_angles, (90, 105, 90))
    assert v.with_bragg == "bragg_present" and v.is_axis_aligned()
    assert np.unravel_index(v.values.argmax(), v.values.shape) == (16, 20, 24)   # hkl = 000
    assert v.provenance["metadata_asserted"]["coordinate_unit"] == "basecell_fractional"


def test_real_patterson_peak_gives_correct_angstrom_length_in_nonorthogonal_cell():
    v = read_discus_volume(FIX / "C3_total.h5")
    truth = np.array([0.30, 0.20, 0.10])            # atom O minus atom C in mono2.cll
    g = [np.linspace(t - 0.1, t + 0.1, 21) for t in truth]
    pts = np.stack(np.meshgrid(*g, indexing="ij"), -1).reshape(-1, 3)
    rho = evaluate_real_space(v, pts).real
    ustar = pts[rho.argmax()]
    assert np.allclose(ustar, truth, atol=0.006)
    res = delta_pdf_fractional(v, require_diffuse=False)
    length = res.cartesian_length(ustar)
    expected = 1.9070                                # sqrt(u.G.u) with beta = 105 deg
    assert length == pytest.approx(expected, abs=0.02)
    orthogonal_wrong = float(np.linalg.norm(ustar * v.cell_lengths))   # ignoring beta
    assert abs(orthogonal_wrong - expected) > 0.1    # the test can tell the two metrics apart


def test_real_fft_axis_spacing_is_1_over_n_dh_not_header_value():
    v = read_discus_volume(FIX / "C3_total.h5")
    res = delta_pdf_fractional(v, require_diffuse=False)
    assert np.diff(res.u_axes[0])[0] == pytest.approx(1 / (33 * 0.25))   # 0.1212, DISCUS 3DPDF header says 0.125


def test_real_skewed_grid_orientation_and_refusal():
    v = read_discus_volume(FIX / "C4_total.h5")
    assert np.allclose(v.increments[:, 1], (0.025, 0.25, 0.0))   # increment vector k is COLUMN k
    assert not v.is_axis_aligned()
    with pytest.raises(ValueError, match="skewed"):
        v.hkl_axes()
    with pytest.raises(ValueError, match="skewed"):
        delta_pdf_fractional(v, require_diffuse=False)


def test_bragg_present_refused_by_default_and_noncubic_has_no_cartesian_q_view():
    v = read_discus_volume(FIX / "C3_total.h5")
    with pytest.raises(ValueError, match="bragg_present"):
        delta_pdf_fractional(v)
    with pytest.raises(ValueError, match="orthogonal"):
        v.to_reciprocal_volume()


# --------------------------------- SYNTHETIC ---------------------------------
def test_fractional_transform_matches_cartesian_transform_in_cubic_cell(tmp_path):
    p = tmp_path / "cubic.h5"
    a, d, n = 5.0, 0.125, 33
    write_discus_like(p, centered_blob(n, d), corner=[-(n // 2) * d] * 3, inc=np.eye(3) * d, cell=(a, a, a))
    v = read_discus_volume(p)
    new = delta_pdf_fractional(v)
    old = delta_pdf_from_diffuse(v.to_reciprocal_volume())
    assert np.allclose(new.complex_volume, old.complex_volume, rtol=1e-10, atol=1e-12 * np.abs(old.complex_volume).max())
    assert all(np.allclose(a * u, r) for u, r in zip(new.u_axes, old.r_axes))


def test_orthorhombic_cartesian_q_axes_use_each_cell_length(tmp_path):
    p = tmp_path / "ortho.h5"
    n, d = 17, 0.25
    write_discus_like(p, centered_blob(n, d), corner=[-(n // 2) * d] * 3, inc=np.eye(3) * d, cell=(4.0, 5.0, 8.0))
    q = read_discus_volume(p).to_reciprocal_volume().q_axes
    assert [round(float(a[1] - a[0]), 9) for a in q] == [round(2 * np.pi * d / c, 9) for c in (4.0, 5.0, 8.0)]


@pytest.mark.parametrize("kw,msg", [
    (dict(axes="uvw"), "uvw"), (dict(unit="angstrom^-1"), "units"),
    (dict(reciprocal="patterson"), "reciprocal space"), (dict(content="3dpdf"), "content"),
    (dict(data_axes=(2, 1, 3)), "data_axes"), (dict(bragg="maybe"), "bragg"),
])
def test_ambiguous_or_wrong_metadata_fails_with_diagnostic(tmp_path, kw, msg):
    p = tmp_path / "bad.h5"
    n, d = 9, 0.25
    write_discus_like(p, centered_blob(n, d), corner=[-(n // 2) * d] * 3, inc=np.eye(3) * d, **kw)
    with pytest.raises(ValueError, match=msg):
        read_discus_volume(p)


def test_uncentered_grid_is_refused(tmp_path):
    p = tmp_path / "shifted.h5"
    n, d = 9, 0.25
    write_discus_like(p, centered_blob(n, d), corner=[0.0] * 3, inc=np.eye(3) * d)
    with pytest.raises(ValueError, match="centered"):
        delta_pdf_fractional(read_discus_volume(p))


def test_missing_group_and_non_hdf5_are_rejected(tmp_path):
    q = tmp_path / "empty.h5"
    with h5py.File(q, "w") as f:
        f.create_dataset("x", data=1)
    with pytest.raises(ValueError, match="common/nexus"):
        read_discus_volume(q)
    t = tmp_path / "t.txt"; t.write_text("nope")
    with pytest.raises(ValueError, match="Not an HDF5"):
        read_discus_volume(t)
