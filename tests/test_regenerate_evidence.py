"""The evidence content hash must ignore DISCUS's run-date stamp and detect any change of the data."""
import sys
from pathlib import Path

import numpy as np
import pytest

h5py = pytest.importorskip("h5py")
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks" / "discus_6_23_01"))
import regenerate_evidence as R  # noqa: E402


def _write(p, values, date):
    with h5py.File(p, "w") as f:
        g = f.create_group("entry").create_group("data")
        g["data_values"] = values
        g["audit_creation_date"] = np.bytes_(date)
        for k in R.META_KEYS:
            g[k] = np.bytes_("x") if k.startswith(("data_type", "data_rad")) else np.arange(3.0)


def test_content_hash_ignores_date_and_detects_changes(tmp_path):
    v = np.random.default_rng(0).random((4, 4, 4))
    _write(tmp_path / "a.h5", v, "2026-10-28"); _write(tmp_path / "b.h5", v, "2031-01-01")
    w = v.copy(); w[1, 2, 3] *= 1 + 1e-15; _write(tmp_path / "c.h5", w, "2026-10-28")
    a, b, c = (R.content_hash(tmp_path / n) for n in ("a.h5", "b.h5", "c.h5"))
    assert a == b and a != c


def test_registry_covers_every_evidence_directory():
    reg = R.load_registry()["evidence"]
    assert len(reg) == 5 and all(len(e["volumes"]) == 4 for e in reg.values())
    for rel, e in reg.items():
        d = R.REPO / rel
        man = R.read_manifest(d)
        assert all(R.sha256_file(d / n) == man[n] for n in e["inputs"])   # inputs to rerun are intact
        assert set(e["original_file_sha256"].items()) <= set(man.items())
