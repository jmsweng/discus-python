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
