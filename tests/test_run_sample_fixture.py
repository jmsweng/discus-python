"""Consistency of the committed GaAs 1D reference run (``examples/run_sample``). No DISCUS needed.

The opt-in engine comparison (``test_integration_discus.py``) can only be trusted if the reference
directory is one coherent run. In October 2026 ``run.mac`` was edited by hand after the run
(``outf result.pdf`` -> ``result.csv``) and nobody noticed, because nothing checked it. These tests
catch that kind of drift in standard CI:

* every file matches ``MANIFEST.sha256``;
* the macro is exactly what the current ``render_pdf_macro`` produces for the recorded parameters;
* the macro's commands are exactly the commands DISCUS echoed in ``run.log``;
* the structure in ``examples/gaas.cll`` and the staged ``model.cll`` have the recorded hash;
* ``result.pdf`` has the r grid the recorded parameters ask for.
"""
import hashlib
import json
import re
import warnings
from pathlib import Path

import numpy as np

from discuspdf import PDFConfig, QmaxLimitWarning, read_pdf, render_pdf_macro

REPO_ROOT = Path(__file__).resolve().parent.parent
RUN = REPO_ROOT / "examples" / "run_sample"
STRUCTURE = REPO_ROOT / "examples" / "gaas.cll"

# Parameters PDFConfig accepts; config.json also holds bookkeeping keys (hashes, launcher, time).
CONFIG_KEYS = ("supercell", "rmax", "dr", "qmax", "qsig", "radiation", "qmin", "dq", "wavelength")


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def recorded_config() -> PDFConfig:
    """PDFConfig of the archived run, with the structure taken from the repository.

    config.json keeps the absolute structure path of the machine the run was made on; the file
    itself is identified by ``input_sha256`` (checked in ``test_structure_hash_matches_record``).
    """
    params = json.loads((RUN / "config.json").read_text(encoding="utf-8"))
    kwargs = {k: params[k] for k in CONFIG_KEYS if k in params}
    kwargs["supercell"] = tuple(kwargs["supercell"])
    with warnings.catch_warnings():  # qmax 30 is above the DISCUS limit at 1.5406 Å; known and documented
        warnings.simplefilter("ignore", QmaxLimitWarning)
        return PDFConfig(structure=STRUCTURE, **kwargs)


def test_files_match_manifest():
    lines = (RUN / "MANIFEST.sha256").read_text(encoding="ascii").splitlines()
    listed = {}
    for line in lines:
        digest, name = line.split(maxsplit=1)
        listed[name] = digest
    assert set(listed) == {"config.json", "model.cll", "result.pdf", "run.log", "run.mac", "status.json"}
    for name, digest in listed.items():
        assert _sha256(RUN / name) == digest, f"{name} differs from MANIFEST.sha256"


def test_structure_hash_matches_record():
    recorded = json.loads((RUN / "config.json").read_text(encoding="utf-8"))["input_sha256"]
    assert _sha256(RUN / "model.cll") == recorded
    assert _sha256(STRUCTURE) == recorded, "examples/gaas.cll is no longer the structure of the reference run"


def test_macro_is_what_the_generator_writes():
    assert (RUN / "run.mac").read_text(encoding="utf-8") == render_pdf_macro(recorded_config())


def test_macro_commands_are_the_commands_in_the_log():
    """DISCUS echoes each command after its prompt (``suite> discus``, ``discus/powd> set qmax,30``, ...).

    The echoed commands must be run.mac line by line. The runner's closing ``exit`` on stdin is echoed
    as ``suite >exit`` (other spacing) and is not part of the macro.
    """
    log = re.sub(r"\x1b\[[0-9;]*m", "", (RUN / "run.log").read_text(encoding="utf-8"))
    echoed = [m.group(1).strip()
              for m in re.finditer(r"^\s*(?:suite|discus(?:/\w+)?)> (\S.*)$", log, flags=re.M)]
    macro = [ln.strip() for ln in (RUN / "run.mac").read_text(encoding="utf-8").splitlines()
             if ln.strip() and not ln.lstrip().startswith("#")]
    assert echoed == macro


def test_result_grid_matches_record():
    cfg = recorded_config()
    r, gr = read_pdf(RUN / "result.pdf")
    n = round(cfg.rmax / cfg.dr) + 1
    assert len(r) == n == 601
    assert np.allclose(r, np.arange(n) * cfg.dr, rtol=0, atol=1e-8)
    assert np.all(np.isfinite(gr))


def test_status_complete():
    assert json.loads((RUN / "status.json").read_text(encoding="utf-8"))["status"] == "complete"
