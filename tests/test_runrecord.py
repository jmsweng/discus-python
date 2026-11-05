"""The shared run-directory helpers keep each runner's status.json layout byte for byte."""
from __future__ import annotations

import json

from discuspdf._runrecord import changed_inputs, sha256_bytes, sha256_file, write_status


def test_status_layouts(tmp_path):
    write_status(tmp_path, "prepared", indent=None)
    assert (tmp_path / "status.json").read_text() == '{"status": "prepared"}\n'
    write_status(tmp_path, "failed", {"reason": "x"}, newline=False)            # pdf, diffuse3d
    assert (tmp_path / "status.json").read_text() == json.dumps({"status": "failed", "reason": "x"}, indent=2)
    write_status(tmp_path, "failed", {"reason": "x"}, indent=1)                 # mc
    assert (tmp_path / "status.json").read_text() == json.dumps({"status": "failed", "reason": "x"}, indent=1) + "\n"
    write_status(tmp_path, "complete", {"a": 1})                                # suite runner
    assert (tmp_path / "status.json").read_text() == json.dumps({"status": "complete", "a": 1}, indent=2) + "\n"


def test_changed_inputs(tmp_path):
    (tmp_path / "a").write_bytes(b"one")
    (tmp_path / "b").write_bytes(b"two")
    hashes = {"a": sha256_file(tmp_path / "a"), "b": sha256_bytes(b"two"), "c": sha256_bytes(b"")}
    assert changed_inputs(tmp_path, hashes) == [("c", "missing")]
    (tmp_path / "b").write_bytes(b"TWO")
    assert changed_inputs(tmp_path, hashes) == [("b", "modified"), ("c", "missing")]
