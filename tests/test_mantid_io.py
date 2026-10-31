"""Mantid MDHisto import/export (Package B, Feature 18). No Mantid needed: every file here was written
by Mantid 7.0.0 (tests/fixtures/mantid_7_0_0, see its README), and expected values were fixed in
expected.json before Mantid ran. roundtrip/ holds files that discuspdf wrote and Mantid loaded and
saved again."""
from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pytest

h5py = pytest.importorskip("h5py")

from discuspdf import Lattice, ReciprocalGrid, ScatteringVolume
from discuspdf.mantid_io import (MantidFormatError, parse_axis_name, read_mantid_mdhisto,
                                 write_mantid_mdhisto)

FIX = Path(__file__).parent / "fixtures" / "mantid_7_0_0"
RT = FIX / "roundtrip"
EXPECTED = json.loads((FIX / "expected.json").read_text())
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks" / "mantid_7_0_0"))
from write_for_mantid import skewed_volume  # noqa: E402


def at(sv: ScatteringVolume, hkl):
    """Index of the grid point closest to hkl, and its distance (r.l.u.)."""
    d = np.linalg.norm(sv.grid.hkl() - np.asarray(hkl, float), axis=-1)
    i = np.unravel_index(np.argmin(d), d.shape)
    return i, float(d[i])


# --------------------------------------------------------------------------- reading Mantid's files
@pytest.mark.parametrize("name", ["aligned_hkl.nxs", "projected_hh0.nxs"])
def test_events_land_at_their_hkl(name):
    sv = read_mantid_mdhisto(FIX / name)
    exp = EXPECTED["files"][name]
    assert list(sv.grid.shape) == exp["shape"]
    for h, k, l, n in exp["events_hkl"]:
        i, dist = at(sv, (h, k, l))
        assert dist < 1e-9                                 # events were placed at bin centres
        assert sv.intensity[i] == n and sv.variance[i] == n
    assert sv.intensity.sum() == sum(e[3] for e in exp["events_hkl"])
    assert all("name matches transform" in c for c in sv.provenance["consistency_checks"])


def test_projection_basis_from_transform():
    sv = read_mantid_mdhisto(FIX / "projected_hh0.nxs")
    V = sv.grid.axis_vectors
    assert np.allclose(V, [[0.25, 0.25, 0], [-0.25, 0.25, 0], [0, 0, 0.5]])
    assert not sv.grid.is_hkl_axis_aligned()


def test_w_matrix_event_lands_at_W_times_x():
    after = read_mantid_mdhisto(FIX / "w_matrix_after.nxs")
    before = read_mantid_mdhisto(FIX / "w_matrix_before.nxs")
    diff = after.intensity - before.intensity
    changed = np.argwhere(diff != 0)
    exp = EXPECTED["files"]["w_matrix_after.nxs"]
    assert len(changed) == 1 and diff[tuple(changed[0])] == exp["weight"]
    assert np.allclose(after.grid.hkl()[tuple(changed[0])], exp["event_hkl"])
    assert np.allclose(after.provenance["W_MATRIX"], [[1, 0, 1], [1, 0, -1], [0, 1, 0]])
    assert any("units length 1.67" in c for c in after.provenance["consistency_checks"])


def test_histo_formula_axis_order_errors_mask_nan():
    sv = read_mantid_mdhisto(FIX / "histo_formula.nxs")
    i0, i1, i2 = np.meshgrid(*[np.arange(n) for n in sv.grid.shape], indexing="ij")
    f = 100.0 * i0 + 10 * i1 + i2 + 1
    ok = ~sv.mask
    assert np.array_equal(sv.intensity[ok], f[ok])
    assert np.allclose(sv.variance[ok], (0.1 * f[ok]) ** 2)
    assert sv.mask[0].all() and sv.mask[2, 3, 1] and sv.mask.sum() == 21
    assert sv.mask_reasons == ("Mantid mask", "NaN in Mantid signal (no coverage)")
    assert np.allclose(sv.grid.corner, (-1.25, -1.0, 0.25)) and np.allclose(sv.grid.axis_vectors, np.diag([0.5] * 3))


def test_lattice_and_ub():
    sv = read_mantid_mdhisto(FIX / "aligned_hkl.nxs")
    assert np.allclose(sv.grid.lattice.parameters, (5, 6, 7, 90, 105, 90))
    UB = np.array(sv.provenance["UB"])
    assert np.allclose(UB.T @ UB, sv.grid.lattice.reciprocal_metric)
    U = np.array(sv.provenance["U"])
    assert np.allclose(U @ U.T, np.eye(3))


# --------------------------------------------------------------------------- refusals
def edited(tmp_path, name, fn):
    q = tmp_path / f"edited_{name}"
    shutil.copyfile(FIX / name, q)
    with h5py.File(q, "r+") as f:
        fn(f["MDHistoWorkspace"])
    return q


def set_attr(ds, key, value):
    ds.attrs[key] = np.bytes_(value.encode())


@pytest.mark.parametrize("name, fn, match", [
    ("projected_hh0.nxs", lambda g: set_attr(g["data/D0"], "long_name", "[H,0,0]"), "says"),
    ("w_matrix_after.nxs", lambda g: set_attr(g["data/D0"], "units", "in 9.99 A^-1"), "units say"),
    ("aligned_hkl.nxs", lambda g: set_attr(g["data/D1"], "units", "Angstrom^-1"), "only 'r.l.u.'"),
    ("aligned_hkl.nxs", lambda g: set_attr(g["data/D2"], "frame", "QSample"), "QSample"),
    ("aligned_hkl.nxs", lambda g: g.attrs.__setitem__("SaveMDVersion", np.int32(1)), "SaveMDVersion"),
    ("aligned_hkl.nxs", lambda g: g["data/D0"].__setitem__(3, g["data/D0"][3] + 0.01), "non-uniform"),
    ("aligned_hkl.nxs", lambda g: g.__delitem__("experiment0/sample/oriented_lattice"), "oriented lattice"),
    ("aligned_hkl.nxs", lambda g: set_attr(g["data/D2"], "frame", "General Frame"), "HKL axes"),
])
def test_refusals(tmp_path, name, fn, match):
    with pytest.raises(MantidFormatError, match=match):
        read_mantid_mdhisto(edited(tmp_path, name, fn))


def test_integrated_extra_axis_is_dropped_and_recorded(tmp_path):
    def add_energy(g):
        d = g["data"]
        for k in ("signal", "errors_squared", "mask", "num_events"):
            a = d[k][()]
            del d[k]
            d.create_dataset(k, data=a[None, ...])
        d["signal"].attrs["axes"] = np.bytes_(b"D3:D2:D1:D0")
        d["signal"].attrs["signal"] = np.int32(1)
        e = d.create_dataset("D3", data=np.array([-1.0, 1.0]))
        set_attr(e, "frame", "General Frame")
        set_attr(e, "long_name", "DeltaE")
        set_attr(e, "units", "meV")
        t = np.eye(5)
        t[:3, :3] = g["transform_to_orig"][()].reshape(4, 4)[:3, :3]
        del g["transform_to_orig"]
        ds = g.create_dataset("transform_to_orig", data=t.ravel().astype(np.float32))
        ds.attrs["rows"], ds.attrs["columns"] = np.int32(5), np.int32(5)
    sv = read_mantid_mdhisto(edited(tmp_path, "aligned_hkl.nxs", add_energy))
    assert sv.grid.shape == (24, 24, 12)
    assert sv.provenance["integrated_axes"] == [{"name": "DeltaE", "frame": "General Frame", "units": "meV",
                                                 "range": [-1.0, 1.0]}]


@pytest.mark.parametrize("name, ok", [("[H,H,0]", [1, 1, 0]), ("[-H,H,0]", [-1, 1, 0]), ("[0,0,L]", [0, 0, 1]),
                                      ("[2H,-H,0]", [2, -1, 0]), ("[0.5H,0,0]", [0.5, 0, 0]), ("[H,1,0]", None),
                                      ("DeltaE", None), ("[H,K]", None)])
def test_parse_axis_name(name, ok):
    got = parse_axis_name(name)
    assert (got is None) if ok is None else np.allclose(got, ok)


# --------------------------------------------------------------------------- writing (checked by Mantid)
def test_mantid_accepted_our_files():
    report = json.loads((RT / "mantid_check.json").read_text())
    assert report["mantid_version"] == "7.0.0" and report["failures"] == []
    assert set(report["files"]) == {"aligned_hkl", "projected_hh0", "w_matrix_after", "histo_formula", "skewed"}
    assert all(all(v for v in e["checks"].values()) for e in report["files"].values())


@pytest.mark.parametrize("name", ["aligned_hkl", "projected_hh0", "w_matrix_after", "histo_formula"])
def test_round_trip_through_mantid(name, tmp_path):
    source = read_mantid_mdhisto(FIX / f"{name}.nxs")
    ours = read_mantid_mdhisto(write_mantid_mdhisto(source, tmp_path / "x.nxs"))
    resaved = read_mantid_mdhisto(RT / f"{name}.resaved.nxs")       # our file, saved again by Mantid
    for other in (ours, resaved):
        assert other.grid.isclose(source.grid, rtol=0, atol=1e-6)
        m = source.mask if source.mask is not None else np.zeros(source.grid.shape, bool)
        assert np.array_equal(other.intensity[~m], source.intensity[~m])
        assert np.array_equal(other.variance[~m], source.variance[~m])
        assert np.array_equal(other.mask if other.mask is not None else np.zeros_like(m), m)


def test_skewed_monoclinic_volume_through_mantid(tmp_path):
    sv = skewed_volume()
    exp = json.loads((RT / "skewed.expected.json").read_text())
    for key, hkl in exp["hkl_of_index"].items():
        assert np.allclose(sv.grid.index_to_hkl([int(x) for x in key.split(",")]), hkl)
    for back in (read_mantid_mdhisto(write_mantid_mdhisto(sv, tmp_path / "s.nxs")),
                 read_mantid_mdhisto(RT / "skewed.resaved.nxs")):
        assert back.grid.isclose(sv.grid, rtol=0, atol=1e-6)          # float32 transform in the file
        assert np.array_equal(back.mask, sv.mask)
        assert np.array_equal(back.intensity[~sv.mask], sv.intensity[~sv.mask])
        assert np.allclose(back.grid.lattice.parameters, sv.grid.lattice.parameters)


def test_writer_layout_and_refusals(tmp_path):
    sv = skewed_volume()
    p = write_mantid_mdhisto(sv, tmp_path / "s.nxs")
    with h5py.File(p, "r") as f:
        d = f["MDHistoWorkspace/data"]
        assert d["signal"].shape == (4, 5, 6) and d["signal"].attrs["axes"] == b"D2:D1:D0"
        assert d["mask"][3, 2, 1] == 1 and d["signal"][3, 2, 1] == 0      # voxel (1,2,3) is masked
        assert d["signal"][1, 2, 3] == sv.intensity[3, 2, 1]
        assert f["MDHistoWorkspace"].attrs["SaveMDVersion"] == 2
        UB = f["MDHistoWorkspace/experiment0/sample/oriented_lattice/orientation_matrix"][()]
        assert np.allclose(UB.T @ UB, sv.grid.lattice.reciprocal_metric)
    with pytest.raises(FileExistsError):
        write_mantid_mdhisto(sv, p)
    plane = ScatteringVolume(ReciprocalGrid.from_ranges(Lattice.from_parameters(4, 4, 4, 90, 90, 90),
                                                        (-1, 1, 5), (-1, 1, 5), (0, 0, 1)), np.ones((5, 5, 1)), "a")
    with pytest.raises(ValueError, match="3D"):
        write_mantid_mdhisto(plane, tmp_path / "p.nxs")


def test_decreasing_axis_written_with_increasing_edges(tmp_path):
    lat = Lattice.from_parameters(4, 4, 4, 90, 90, 90)
    g = ReciprocalGrid(lat, (1, 0, 0), [[-0.25, 0, 0], [0, 0.25, 0], [0, 0, 0.5]], (5, 4, 3))
    sv = ScatteringVolume(g, np.arange(60.0).reshape(5, 4, 3), "a")
    back = read_mantid_mdhisto(write_mantid_mdhisto(sv, tmp_path / "d.nxs"))
    assert back.grid.isclose(g, rtol=0, atol=1e-6) and np.array_equal(back.intensity, sv.intensity)
