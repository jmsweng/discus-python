"""E1/E3 software tests for Diffuse3DRunner (fake executable; no live DISCUS)."""
import json
import sys
from pathlib import Path

import numpy as np
import pytest

from discuspdf import Diffuse3DRunner, OutputSpec


def make_hdf5(path: Path, *, n=8, dq=.2, nan=False, intensity_shape=None,
              q_lengths=None) -> None:
    import h5py
    shape = intensity_shape or (n, n, n)
    q = (np.arange(n) - n // 2) * dq
    with h5py.File(path, "w") as f:
        f.create_dataset("intensity", data=np.ones(shape) if not nan
                         else np.full(shape, np.nan))
        lengths = q_lengths or (n, n, n)
        for i, key in enumerate(("qx", "qy", "qz")):
            f.create_dataset(key, data=(np.arange(lengths[i]) - lengths[i] // 2) * dq)


def fake_runner(tmp_path: Path, script_body: str) -> tuple[Diffuse3DRunner, Path]:
    script = tmp_path / "fake_discus.py"
    script.write_text(script_body)
    macro = tmp_path / "input.mac"
    macro.write_text("# validated 3D macro\n")
    runner = Diffuse3DRunner(runs_directory=tmp_path / "runs",
                             launcher=[sys.executable, str(script)])
    return runner, macro


def spec() -> OutputSpec:
    return OutputSpec(intensity_dataset="intensity",
                      q_datasets=("qx", "qy", "qz"),
                      expected_shape=(8, 8, 8))


def status_of(run) -> dict:
    return json.loads((run.directory / "status.json").read_text())


# --- E1 acceptance cases -------------------------------------------------

def test_valid_hdf5_with_spec_is_complete(tmp_path):
    body = ("import sys\n"
            "sys.path.insert(0, sys.argv[1])\n"
            "import h5py, numpy as np\n"
            "with h5py.File('diffuse.h5', 'w') as f:\n"
            "    f.create_dataset('intensity', data=np.ones((8,8,8)))\n"
            "    for k in ('qx','qy','qz'):\n"
            "        f.create_dataset(k, data=(np.arange(8)-4)*.2)\n")
    runner, macro = fake_runner(tmp_path, body)
    run = runner.prepare(macro)
    result = runner.execute(run, output_spec=spec())
    st = status_of(result)
    assert st["status"] == "complete"
    assert st["returncode"] == 0
    assert st["output_validation"]["valid_hdf5"] is True
    assert st["output_validation"]["intensity_shape"] == [8, 8, 8]
    # The exact staged macro is the exact executed macro.
    assert (result.directory / "run.mac").read_bytes() == macro.read_bytes()


def test_non_hdf5_bytes_are_failed(tmp_path):
    body = "from pathlib import Path\nPath('diffuse.h5').write_bytes(b'not-hdf5-bytes')\n"
    runner, macro = fake_runner(tmp_path, body)
    run = runner.prepare(macro)
    with pytest.raises(RuntimeError):
        runner.execute(run, output_spec=spec())
    st = status_of(run)
    assert st["status"] == "failed"
    assert st["returncode"] == 0
    assert "not a valid HDF5" in st["reason"]


def test_absent_hdf5_is_failed(tmp_path):
    body = "print('engine ran but wrote nothing')\n"
    runner, macro = fake_runner(tmp_path, body)
    run = runner.prepare(macro)
    with pytest.raises(RuntimeError, match="no output"):
        runner.execute(run)
    st = status_of(run)
    assert st["status"] == "failed"
    assert st["returncode"] == 0


def test_valid_hdf5_without_spec_is_produced_unverified(tmp_path):
    body = ("import sys\n"
            "sys.path.insert(0, sys.argv[1])\n"
            "import h5py, numpy as np\n"
            "with h5py.File('diffuse.h5', 'w') as f:\n"
            "    f.create_dataset('intensity', data=np.ones((8,8,8)))\n"
            "    for k in ('qx','qy','qz'):\n"
            "        f.create_dataset(k, data=(np.arange(8)-4)*.2)\n")
    runner, macro = fake_runner(tmp_path, body)
    run = runner.prepare(macro)
    result = runner.execute(run)
    st = status_of(result)
    assert st["status"] == "produced_unverified"
    assert st["returncode"] == 0


def test_malformed_datasets_are_failed(tmp_path):
    # 2D intensity + Q axis length mismatch.
    body = ("import sys\n"
            "sys.path.insert(0, sys.argv[1])\n"
            "import h5py, numpy as np\n"
            "with h5py.File('diffuse.h5', 'w') as f:\n"
            "    f.create_dataset('intensity', data=np.ones((8,8)))\n"
            "    for k in ('qx','qy','qz'):\n"
            "        f.create_dataset(k, data=(np.arange(8)-4)*.2)\n")
    runner, macro = fake_runner(tmp_path, body)
    run = runner.prepare(macro)
    with pytest.raises(RuntimeError):
        runner.execute(run, output_spec=spec())
    st = status_of(run)
    assert st["status"] == "failed"
    assert "3D" in st["reason"]

    # Right dimensionality, wrong expected shape.
    body2 = ("import sys\n"
             "sys.path.insert(0, sys.argv[1])\n"
             "import h5py, numpy as np\n"
             "with h5py.File('diffuse.h5', 'w') as f:\n"
             "    f.create_dataset('intensity', data=np.ones((8,8,4)))\n"
             "    for k in ('qx','qy','qz'):\n"
             "        f.create_dataset(k, data=(np.arange(8)-4)*.2)\n")
    runner2, _ = fake_runner(tmp_path, body2)
    run2 = runner2.prepare(macro)
    with pytest.raises(RuntimeError):
        runner2.execute(run2, output_spec=spec())
    st2 = status_of(run2)
    assert st2["status"] == "failed"
    assert "shape" in st2["reason"]


def test_nan_output_is_failed(tmp_path):
    body = ("import sys\n"
            "sys.path.insert(0, sys.argv[1])\n"
            "import h5py, numpy as np\n"
            "with h5py.File('diffuse.h5', 'w') as f:\n"
            "    f.create_dataset('intensity', data=np.full((8,8,8), np.nan))\n"
            "    for k in ('qx','qy','qz'):\n"
            "        f.create_dataset(k, data=(np.arange(8)-4)*.2)\n")
    runner, macro = fake_runner(tmp_path, body)
    run = runner.prepare(macro)
    with pytest.raises(RuntimeError):
        runner.execute(run, output_spec=spec())
    st = status_of(run)
    assert st["status"] == "failed"
    assert "NaN/Inf" in st["reason"]


def test_nonzero_exit_despite_output_is_failed(tmp_path):
    body = ("import sys\n"
            "sys.path.insert(0, sys.argv[1])\n"
            "import h5py, numpy as np\n"
            "with h5py.File('diffuse.h5', 'w') as f:\n"
            "    f.create_dataset('intensity', data=np.ones((8,8,8)))\n"
            "    for k in ('qx','qy','qz'):\n"
            "        f.create_dataset(k, data=(np.arange(8)-4)*.2)\n"
            "sys.exit(3)\n")
    runner, macro = fake_runner(tmp_path, body)
    run = runner.prepare(macro)
    with pytest.raises(RuntimeError, match="failed"):
        runner.execute(run, output_spec=spec())
    st = status_of(run)
    assert st["status"] == "failed"
    assert st["returncode"] == 3


def test_engine_error_log_despite_zero_exit_is_failed(tmp_path):
    body = ("import sys\n"
            "sys.path.insert(0, sys.argv[1])\n"
            "import h5py, numpy as np\n"
            "with h5py.File('diffuse.h5', 'w') as f:\n"
            "    f.create_dataset('intensity', data=np.ones((8,8,8)))\n"
            "    for k in ('qx','qy','qz'):\n"
            "        f.create_dataset(k, data=(np.arange(8)-4)*.2)\n"
            "print('FATAL ERROR: dataset truncated')\n")
    runner, macro = fake_runner(tmp_path, body)
    run = runner.prepare(macro)
    with pytest.raises(RuntimeError, match="failed"):
        runner.execute(run, output_spec=spec())
    st = status_of(run)
    assert st["status"] == "failed"
    assert st["returncode"] == 0
    assert st["engine_error"] == "FATAL ERROR"
    # The full log is preserved for inspection.
    log = (run.directory / "run.log").read_text()
    assert "FATAL ERROR: dataset truncated" in log


def test_pgplot_x11_noise_is_not_fatal(tmp_path):
    body = ("import sys\n"
            "sys.path.insert(0, sys.argv[1])\n"
            "import h5py, numpy as np\n"
            "with h5py.File('diffuse.h5', 'w') as f:\n"
            "    f.create_dataset('intensity', data=np.ones((8,8,8)))\n"
            "    for k in ('qx','qy','qz'):\n"
            "        f.create_dataset(k, data=(np.arange(8)-4)*.2)\n"
            "print('xdpyinfo: unable to open display \"\"')\n"
            "print('pgxwin_server: cannot connect to X server []')\n")
    runner, macro = fake_runner(tmp_path, body)
    run = runner.prepare(macro)
    result = runner.execute(run, output_spec=spec())
    assert status_of(result)["status"] == "complete"


def test_missing_binary_is_failed_without_engine(tmp_path):
    runner = Diffuse3DRunner(runs_directory=tmp_path / "runs",
                             launcher=["definitely_not_a_binary_xyz"])
    macro = tmp_path / "input.mac"
    macro.write_text("# macro\n")
    run = runner.prepare(macro)
    with pytest.raises(RuntimeError):
        runner.execute(run)
    st = status_of(run)
    assert st["status"] == "failed"
    assert st["engine_started"] is False
    assert "returncode" not in st


def test_timeout_is_failed_without_engine(tmp_path):
    body = "import time\ntime.sleep(30)\n"
    runner, macro = fake_runner(tmp_path, body)
    runner.timeout_seconds = 1
    run = runner.prepare(macro)
    with pytest.raises(RuntimeError):
        runner.execute(run)
    st = status_of(run)
    assert st["status"] == "failed"
    assert st["engine_started"] is False


# --- E3: preflight bookkeeping --------------------------------------------

def test_tampered_staged_input_fails_preflight(tmp_path):
    runner, macro = fake_runner(tmp_path, "print('ok')\n")
    asset = tmp_path / "model.cll"
    asset.write_text("input")
    run = runner.prepare(macro, (asset,))
    (run.directory / "model.cll").write_text("tampered")
    with pytest.raises(ValueError, match="modified"):
        runner.execute(run)
    st = status_of(run)
    assert st["status"] == "failed"
    assert st["engine_started"] is False
    assert "returncode" not in st


def test_tampered_macro_fails_preflight(tmp_path):
    runner, macro = fake_runner(tmp_path, "print('ok')\n")
    run = runner.prepare(macro)
    (run.directory / "run.mac").write_text("tampered macro")
    with pytest.raises(ValueError, match="modified"):
        runner.execute(run)
    st = status_of(run)
    assert st["status"] == "failed"
    assert st["engine_started"] is False


def test_rerun_rejected_with_failed_status(tmp_path):
    body = ("import sys\n"
            "sys.path.insert(0, sys.argv[1])\n"
            "import h5py, numpy as np\n"
            "with h5py.File('diffuse.h5', 'w') as f:\n"
            "    f.create_dataset('intensity', data=np.ones((8,8,8)))\n"
            "    for k in ('qx','qy','qz'):\n"
            "        f.create_dataset(k, data=(np.arange(8)-4)*.2)\n")
    runner, macro = fake_runner(tmp_path, body)
    run = runner.prepare(macro)
    runner.execute(run, output_spec=spec())
    with pytest.raises(FileExistsError):
        runner.execute(run)
    st = status_of(run)
    assert st["status"] == "failed"
    assert st["engine_started"] is False


def test_missing_inputs_json_fails_preflight(tmp_path):
    runner, macro = fake_runner(tmp_path, "print('ok')\n")
    run = runner.prepare(macro)
    (run.directory / "inputs.json").unlink()
    with pytest.raises(ValueError, match="inputs.json"):
        runner.execute(run)
    assert status_of(run)["status"] == "failed"


def test_output_spec_validation():
    with pytest.raises(ValueError):
        OutputSpec(intensity_dataset="i", q_datasets=("qx", "qy"))
    with pytest.raises(ValueError):
        OutputSpec(intensity_dataset="i", q_datasets=("qx", "", "qz"))
    with pytest.raises(ValueError):
        OutputSpec(intensity_dataset="i", q_datasets=("qx", "qy", "qz"),
                   expected_shape=(8, 8))
    s = OutputSpec(intensity_dataset="i", q_datasets=("qx", "qy", "qz"))
    assert s.expected_shape is None
