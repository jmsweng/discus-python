import numpy as np
import pytest

from discuspdf import ReciprocalVolume, Diffuse3DRunner


def test_volume_metadata_and_shape():
    q = (np.arange(8) - 4) * .2
    v = ReciprocalVolume(np.ones((8, 8, 8)), (q, q, q),
                         "cartesian_q", "1/angstrom",
                         preprocessing=("bragg_subtracted",))
    v.require_bragg_subtracted()
    with pytest.raises(ValueError, match="Bragg"):
        ReciprocalVolume(v.intensity, v.q_axes, "cartesian_q", "1/angstrom").require_bragg_subtracted()


@pytest.mark.parametrize("system,unit", [("hkl", "rlu"), ("cartesian_q", "rlu")])
def test_reject_unverified_units(system, unit):
    q = (np.arange(8) - 4) * .2
    with pytest.raises(ValueError):
        ReciprocalVolume(np.ones((8, 8, 8)), (q, q, q), system, unit)


def test_reject_offset_axes():
    q = np.arange(8) * .2
    with pytest.raises(ValueError, match="centered"):
        ReciprocalVolume(np.ones((8, 8, 8)), (q, q, q), "cartesian_q", "1/angstrom")


def test_prepare_isolated_3d_macro(tmp_path):
    script = tmp_path / "input.mac"
    script.write_text("# independently validated macro goes here\n")
    asset = tmp_path / "model.cll"
    asset.write_text("input")
    runner = Diffuse3DRunner(runs_directory=tmp_path / "runs")
    prepared = runner.prepare(script, (asset,))
    assert (prepared.directory / "run.mac").read_text() == script.read_text()
    (prepared.directory / "model.cll").write_text("tampered")
    with pytest.raises(ValueError, match="modified"):
        runner.execute(prepared)
    # The failed preflight is recorded; the directory does not stay "prepared".
    import json
    assert json.loads((prepared.directory / "status.json").read_text())["status"] == "failed"


def test_read_hdf5_volume_validation(tmp_path):
    import h5py
    from discuspdf import read_hdf5_volume
    q = (np.arange(8) - 4) * .2
    with h5py.File(tmp_path / "v.h5", "w") as f:
        f.create_dataset("intensity", data=np.ones((8, 8, 8)))
        for key in ("qx", "qy", "qz"):
            f.create_dataset(key, data=q)
    intensity, axes = read_hdf5_volume(tmp_path / "v.h5", "intensity", ("qx", "qy", "qz"))
    assert intensity.shape == (8, 8, 8)
    assert all(np.allclose(a, q) for a in axes)
    # Non-HDF5 bytes are rejected with a clear message.
    (tmp_path / "bad.h5").write_bytes(b"garbage")
    with pytest.raises(ValueError, match="not a valid HDF5"):
        read_hdf5_volume(tmp_path / "bad.h5", "intensity", ("qx", "qy", "qz"))
    # Missing dataset is rejected.
    with pytest.raises(ValueError, match="not found"):
        read_hdf5_volume(tmp_path / "v.h5", "missing", ("qx", "qy", "qz"))
    # Mismatched axis length is rejected.
    with h5py.File(tmp_path / "bad2.h5", "w") as f:
        f.create_dataset("intensity", data=np.ones((8, 8, 8)))
        f.create_dataset("qx", data=q)
        f.create_dataset("qy", data=q)
        f.create_dataset("qz", data=np.arange(4) * .2)
    with pytest.raises(ValueError, match="does not match"):
        read_hdf5_volume(tmp_path / "bad2.h5", "intensity", ("qx", "qy", "qz"))
