"""Atomic text-file writing (shared by ``discus_io`` and ``cif_io``).

The logic is moved here from ``discus_io.py`` unchanged: refuse an existing
target unless ``overwrite=True``; write to a temporary file in the same
directory; ``os.replace``; on failure delete the temp file, leave the target
untouched.  The written text is ASCII with ``\\n`` line endings.
"""
from __future__ import annotations

import os
import tempfile
from pathlib import Path


def atomic_write_text(
    path: Path | str,
    text: str,
    *,
    overwrite: bool = False,
) -> Path:
    """Write ``text`` to ``path`` atomically.

    Refuses an existing path unless ``overwrite=True``.  Writes to a temporary
    file in the same directory, then ``os.replace``.  The written text is ASCII
    with ``\\n`` line endings.  On any failure the temporary file is deleted and
    the target is left untouched.
    """
    target = Path(path)
    if target.exists() and not overwrite:
        raise FileExistsError(f"refusing to overwrite existing file: {target}")

    try:
        text.encode("ascii")
    except UnicodeEncodeError as e:
        raise ValueError(f"non-ASCII text: {e}") from e

    directory = target.parent if str(target.parent) else Path(".")
    fd, tmp_name = tempfile.mkstemp(dir=str(directory), prefix=".tmp_", suffix=".tmp")
    os.close(fd)
    tmp = Path(tmp_name)
    try:
        with open(tmp, "w", newline="\n", encoding="ascii") as f:
            f.write(text)
        os.replace(tmp, target)
    except Exception:
        if tmp.exists():
            os.unlink(tmp)
        raise
    return target
