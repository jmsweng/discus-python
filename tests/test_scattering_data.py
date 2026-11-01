"""Tests for ScatteringVolume, its xarray adapter and HDF5 persistence (Package B, Features 7 and 19).

Independent references: the existing delta_pdf_fractional / delta_pdf_from_diffuse transforms on the
real DISCUS C3 fixture and on synthetic data (the adapter must reproduce them exactly), and a
brute-force enumeration of reflections for the Bragg mask.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from discuspdf import (Lattice, ReciprocalGrid, ScatteringVolume, delta_pdf_fractional,
                       delta_pdf_from_diffuse, mask_near_bragg)

h5py = pytest.importorskip("h5py")
FIX = Path(__file__).parent / "fixtures" / "discus_6.23.01"


def ortho_grid(shape=(8, 9, 10), steps=(0.25, 0.25, 0.25)):
    return ReciprocalGrid.centered(Lattice.from_parameters(4.0, 5.0, 6.0, 90, 90, 90), shape, steps)


def mono_grid():
    return ReciprocalGrid.centered(Lattice.from_parameters(5, 6, 7, 90, 105, 90), (9, 8, 7), (0.25, 0.25, 0.25))


def sym_data(shape, seed=0):
    """Real, Friedel-symmetric data on a centred grid (so the transform is real)."""
    rng = np.random.default_rng(seed)
    a = rng.normal(size=shape)
    un = np.fft.ifftshift(a)
    opp = un[np.ix_(*[(-np.arange(n)) % n for n in shape])]
    return np.fft.fftshift(0.5 * (un + opp))


def c3():
    from discuspdf import read_discus_volume
    return read_discus_volume(FIX / "C3_total.h5")


# --------------------------------------------------------------------------- validation rules
def test_missing_data_must_be_masked():
    g = ortho_grid()
    I = np.ones(g.shape)
    I[1, 2, 3] = np.nan
    with pytest.raises(ValueError, match="unmasked"):
        ScatteringVolume(g, I, "arbitrary")
    m = np.zeros(g.shape, bool)
    m[1, 2, 3] = True
    v = ScatteringVolume(g, I, "arbitrary", mask=m, mask_reasons=("detector gap",))
    assert np.isnan(v.intensity[1, 2, 3])


@pytest.mark.parametrize("kwargs, match", [
    ({"units": ""}, "units"),
    ({"bragg": "subtracted"}, "bragg"),
    ({"bragg": "bragg_masked"}, "bragg_masked"),
    ({"mask": np.zeros((8, 9, 10), dtype=np.uint8)}, "boolean"),
    ({"mask": np.zeros((8, 9, 9), dtype=bool)}, "shape"),
    ({"variance": -np.ones((8, 9, 10))}, "variance"),
    ({"variance": np.ones((8, 9))}, "variance"),
    ({"mask_reasons": ("x",)}, "without a mask"),
    ({"provenance": {"x": float("nan")}}, "JSON"),
    ({"provenance": {"x": object()}}, "JSON"),
    ({"corrections": ("",)}, "non-empty"),
])
def test_invalid_metadata(kwargs, match):
    args = {"units": "arbitrary", **kwargs}
    with pytest.raises(ValueError, match=match):
        ScatteringVolume(ortho_grid(), np.zeros((8, 9, 10)), **args)


def test_shape_and_dtype():
    with pytest.raises(ValueError, match="shape"):
        ScatteringVolume(ortho_grid(), np.zeros((8, 9, 9)), "a")
    with pytest.raises(ValueError, match="numeric"):
        ScatteringVolume(ortho_grid(), np.zeros((8, 9, 10), dtype=complex), "a")
    with pytest.raises(TypeError):
        ScatteringVolume("grid", np.zeros((8, 9, 10)), "a")


def test_immutable_arrays_and_negative_values_kept():
    I = -np.arange(720, dtype=float).reshape(8, 9, 10)
    v = ScatteringVolume(ortho_grid(), I, "a", bragg="bragg_subtracted")
    assert v.intensity.min() == -719.0
    with pytest.raises(ValueError):
        v.intensity[0, 0, 0] = 1
    I[0, 0, 1] = 5
    assert v.intensity[0, 0, 1] == -1.0          # copied, not aliased


# --------------------------------------------------------------------------- Bragg mask
def brute_bragg_mask(grid, radius):
    B = 2 * np.pi * grid.lattice.reciprocal_matrix
    q = grid.q().reshape(-1, 3)
    lo = np.floor(grid.hkl().reshape(-1, 3).min(0)) - 1
    hi = np.ceil(grid.hkl().reshape(-1, 3).max(0)) + 1
    refl = np.array(np.meshgrid(*[np.arange(a, b + 1) for a, b in zip(lo, hi)], indexing="ij")).reshape(3, -1).T
    qr = refl @ B
    d2 = ((q[:, None, :] - qr[None, :, :]) ** 2).sum(-1).min(1)
    return (d2 <= radius ** 2).reshape(grid.shape)


@pytest.mark.parametrize("grid", [ortho_grid((9, 9, 9), (0.125, 0.125, 0.125)),
                                  ReciprocalGrid(Lattice.from_parameters(5, 6, 7, 80, 105, 70), (-1.3, -1, -1.1),
                                                 [[0.1, 0.02, 0], [0, 0.11, 0.01], [0.015, 0, 0.12]], (21, 19, 17))])
def test_bragg_mask_matches_brute_force(grid):
    r = 0.25
    m = mask_near_bragg(grid, r)
    assert m.any() and not m.all()
    assert np.array_equal(m, brute_bragg_mask(grid, r))


def test_bragg_mask_radius_limits_and_labels():
    g = ortho_grid((9, 9, 9), (0.125, 0.125, 0.125))
    with pytest.raises(ValueError, match="radius"):
        mask_near_bragg(g, 2 * np.pi / 6)            # half of shortest |c*| = pi/6 ... too large
    v = ScatteringVolume(g, np.ones(g.shape), "a", bragg="bragg_present").mask_bragg(0.2)
    assert v.bragg == "bragg_masked" and v.mask_reasons[0].startswith("bragg:")
    assert v.mask[4, 4, 4]                           # hkl = 0
    sub = ScatteringVolume(g, np.ones(g.shape), "a", bragg="bragg_subtracted")
    with pytest.raises(ValueError, match="already"):
        sub.mask_bragg(0.2)


def test_with_mask_union_and_correction_record():
    g = ortho_grid()
    a = np.zeros(g.shape, bool)
    a[0] = True
    b = np.zeros(g.shape, bool)
    b[:, 0] = True
    v = ScatteringVolume(g, np.ones(g.shape), "a").with_mask(a, "edge").with_mask(b, "gap")
    assert np.array_equal(v.mask, a | b) and v.mask_reasons == ("edge", "gap")
    w = v.with_correction(v.intensity - 0.5, "background: constant 0.5 subtracted")
    assert w.corrections == ("background: constant 0.5 subtracted",)
    assert np.allclose(w.intensity, 0.5) and np.array_equal(w.mask, v.mask)


# --------------------------------------------------------------------------- transforms (parity)
def test_delta_pdf_reproduces_delta_pdf_fractional_on_real_c3():
    vol = c3()
    sv = ScatteringVolume.from_discus_volume(vol)
    assert sv.bragg == "bragg_present" and sv.provenance["sha256"] == vol.provenance["sha256"]
    with pytest.raises(ValueError, match="Patterson"):
        sv.delta_pdf()
    ours = sv.delta_pdf(require_diffuse=False)
    ref = delta_pdf_fractional(vol, require_diffuse=False)
    # Equal to rounding: beta comes back from the lattice matrix as 105.00000000000001.
    assert np.allclose(ours.complex_volume, ref.complex_volume, rtol=0, atol=1e-13 * np.abs(ref.complex_volume).max())
    for a, b in zip(ours.u_axes, ref.u_axes):
        assert np.array_equal(a, b)
    assert np.allclose(ours.direct_basis, ref.direct_basis)


def test_cartesian_and_fractional_paths_agree_on_orthogonal_cell():
    g = ortho_grid((8, 9, 10), (0.25, 0.2, 0.125))
    v = ScatteringVolume(g, sym_data(g.shape, 3), "a", bragg="bragg_subtracted")
    frac = v.delta_pdf()
    cart = delta_pdf_from_diffuse(v.to_reciprocal_volume())
    assert np.allclose(frac.complex_volume, cart.complex_volume, rtol=1e-12, atol=1e-15)
    for u, r, a in zip(frac.u_axes, cart.r_axes, (4.0, 5.0, 6.0)):
        assert np.allclose(u * a, r)


def test_monoclinic_has_no_cartesian_path():
    v = ScatteringVolume(mono_grid(), sym_data((9, 8, 7)), "a", bragg="bragg_subtracted")
    with pytest.raises(ValueError, match="Cartesian"):
        v.to_reciprocal_volume()
    assert np.isfinite(v.delta_pdf().complex_volume).all()


def test_masked_voxels_need_explicit_fill():
    g = ortho_grid((8, 9, 10))
    I = sym_data(g.shape, 5)
    m = mask_near_bragg(g, 0.3)
    I_nan = np.where(m, np.nan, I)
    v = ScatteringVolume(g, I_nan, "a", bragg="bragg_masked", mask=m, mask_reasons=("bragg",))
    with pytest.raises(ValueError, match="fill_masked"):
        v.delta_pdf()
    with pytest.raises(ValueError, match="fill_masked"):
        v.to_reciprocal_volume()
    res = v.delta_pdf(fill_masked=0.0)
    filled = np.where(m, 0.0, I)
    ref = delta_pdf_from_diffuse(filled, ortho_grid((8, 9, 10)).q_axes())
    assert np.allclose(res.complex_volume, ref.complex_volume)
    assert "masked voxels filled with 0" in res.provenance["corrections"]
    rv = v.to_reciprocal_volume(fill_masked=0.0)
    assert rv.preprocessing == ("bragg_excluded",)
    with pytest.raises(ValueError, match="finite"):
        v.delta_pdf(fill_masked=np.nan)


def test_skewed_grid_refuses_transforms():
    pytest.importorskip("h5py")
    from discuspdf import read_discus_volume
    sv = ScatteringVolume.from_discus_volume(read_discus_volume(FIX / "C4_total.h5"))
    with pytest.raises(ValueError, match="skewed"):
        sv.delta_pdf(require_diffuse=False)


# --------------------------------------------------------------------------- xarray
xr = pytest.importorskip("xarray")


def full_volume():
    g = ortho_grid((8, 9, 10))
    rng = np.random.default_rng(1)
    m = np.zeros(g.shape, bool)
    m[0, :, 0] = True
    I = np.where(m, np.nan, rng.normal(size=g.shape))
    return ScatteringVolume(g, I, "e^2 per crystal", bragg="bragg_subtracted",
                            variance=np.where(m, np.nan, rng.random(g.shape)), mask=m,
                            mask_reasons=("edge",), corrections=("scaled x2",),
                            provenance={"engine": "DISCUS 6.23.01", "sha256": "abc"})


def same(a: ScatteringVolume, b: ScatteringVolume) -> bool:
    eq = lambda x, y: (x is None and y is None) or np.array_equal(x, y, equal_nan=True)
    return (a.grid.isclose(b.grid) and eq(a.intensity, b.intensity) and eq(a.variance, b.variance)
            and eq(a.mask, b.mask) and a.units == b.units and a.bragg == b.bragg
            and a.mask_reasons == b.mask_reasons and a.corrections == b.corrections
            and a.provenance == b.provenance)


def test_xarray_round_trip_and_coords():
    v = full_volume()
    ds = v.to_xarray()
    assert ds["intensity"].dims == ("i0", "i1", "i2") and ds["mask"].dtype == bool
    assert np.allclose(ds["h"].values, v.grid.hkl_axes()[0])
    assert np.allclose(ds["qy"].values, v.grid.q_axes()[1])
    assert ds["intensity"].attrs["units"] == "e^2 per crystal"
    assert same(ScatteringVolume.from_xarray(ds), v)
    # label-based selection works after swapping dims
    val = ds.swap_dims(i0="h", i1="k", i2="l")["intensity"].sel(h=0.25, k=-0.5, l=0.0).item()
    i = tuple(int(x) for x in np.round(v.grid.hkl_to_index((0.25, -0.5, 0.0))))
    assert val == v.intensity[i]


@pytest.mark.parametrize("sel", [dict(i0=slice(2, 7, 2)), dict(i0=slice(None, None, -1), i2=slice(1, 9, 3)),
                                 dict(i2=[4]), dict(i1=slice(3, 4))])
def test_xarray_slicing_rebuilds_the_grid(sel):
    v = full_volume()
    sub = ScatteringVolume.from_xarray(v.to_xarray().isel(**sel))
    idx = tuple(sel.get(d, slice(None)) for d in ("i0", "i1", "i2"))
    hkl = v.grid.hkl()
    expected = hkl[idx[0]][:, idx[1]][:, :, idx[2]]
    assert np.allclose(sub.grid.hkl(), expected)
    assert np.array_equal(sub.intensity, v.intensity[idx[0]][:, idx[1]][:, :, idx[2]], equal_nan=True)


def test_xarray_refuses_stale_or_reordered_data():
    ds = full_volume().to_xarray()
    with pytest.raises(ValueError, match="h does not match"):
        ScatteringVolume.from_xarray(ds.assign_coords(h=("i0", ds["h"].values + 0.1)))
    with pytest.raises(ValueError, match="transpose"):
        ScatteringVolume.from_xarray(ds.transpose("i2", "i1", "i0"))
    with pytest.raises(ValueError, match="keep_attrs"):
        bare = ds.copy()
        bare.attrs = {}
        ScatteringVolume.from_xarray(bare)
    with pytest.raises(ValueError, match="boolean"):          # arithmetic turns the mask into numbers
        ScatteringVolume.from_xarray(ds * 1)
    with pytest.raises(ValueError, match="evenly"):
        ScatteringVolume.from_xarray(ds.isel(i0=[0, 1, 3]))
    with pytest.raises(ValueError, match="version"):
        ScatteringVolume.from_xarray(ds.assign_attrs(discuspdf_schema_version=2))


def test_xarray_skewed_grid_round_trip():
    from discuspdf import read_discus_volume
    sv = ScatteringVolume.from_discus_volume(read_discus_volume(FIX / "C4_total.h5"))
    ds = sv.to_xarray()
    assert "h" not in ds.coords and "qx" not in ds.coords
    assert same(ScatteringVolume.from_xarray(ds), sv)


# --------------------------------------------------------------------------- HDF5 persistence
@pytest.mark.parametrize("compression", [None, "gzip"])
def test_hdf5_round_trip(tmp_path, compression):
    v = full_volume()
    p = v.save_hdf5(tmp_path / "v.h5", compression=compression)
    assert same(ScatteringVolume.load_hdf5(p), v)
    with h5py.File(p, "r") as f:
        assert f.attrs["schema"] == "discuspdf.ScatteringVolume" and f.attrs["schema_version"] == 1
        assert f["mask"].dtype == np.uint8


def test_hdf5_real_discus_round_trip_keeps_provenance(tmp_path):
    sv = ScatteringVolume.from_discus_volume(c3())
    back = ScatteringVolume.load_hdf5(sv.save_hdf5(tmp_path / "c3.h5"))
    assert same(back, sv)
    assert back.provenance["discus_metadata"]["audit_creation_method"] == "6.23.01"


def test_hdf5_refusals(tmp_path):
    v = full_volume()
    p = v.save_hdf5(tmp_path / "v.h5")
    with pytest.raises(FileExistsError):
        v.save_hdf5(p)
    v.save_hdf5(p, overwrite=True)
    with pytest.raises(ValueError, match="compression"):
        v.save_hdf5(tmp_path / "w.h5", compression="lzf")

    def edited(name, fn):
        q = tmp_path / name
        q.write_bytes(p.read_bytes())
        with h5py.File(q, "r+") as f:
            fn(f)
        return q

    with pytest.raises(ValueError, match="version 2"):
        ScatteringVolume.load_hdf5(edited("v2.h5", lambda f: f.attrs.__setitem__("schema_version", 2)))
    with pytest.raises(ValueError, match="schema"):
        ScatteringVolume.load_hdf5(edited("other.h5", lambda f: f.attrs.__setitem__("schema", "x")))
    with pytest.raises(ValueError, match="SHA-256"):
        ScatteringVolume.load_hdf5(edited("tamper.h5", lambda f: f["intensity"].__setitem__((1, 1, 1), 99.0)))
    with pytest.raises(ValueError, match="variance does not match"):
        ScatteringVolume.load_hdf5(edited("tvar.h5", lambda f: f["variance"].__setitem__((1, 1, 1), 9.0)))
    with pytest.raises(ValueError, match="mask does not match"):
        ScatteringVolume.load_hdf5(edited("tmask.h5", lambda f: f["mask"].__setitem__((5, 5, 5), 1)))
    with pytest.raises(ValueError, match="0 and 1"):
        ScatteringVolume.load_hdf5(edited("mask.h5", lambda f: f["mask"].__setitem__((1, 1, 1), 2)))
    with pytest.raises(ValueError, match="missing"):
        ScatteringVolume.load_hdf5(edited("noattr.h5", lambda f: f.attrs.__delitem__("bragg")))

    def strings(f):
        del f["variance"]
        f.create_dataset("variance", data=np.array([b"x"] * 3))
    with pytest.raises(ValueError, match="non-numeric"):
        ScatteringVolume.load_hdf5(edited("str.h5", strings))
    (tmp_path / "text.h5").write_text("not hdf5")
    with pytest.raises(ValueError, match="not an HDF5"):
        ScatteringVolume.load_hdf5(tmp_path / "text.h5")
