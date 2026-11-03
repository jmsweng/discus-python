"""Back end of the discus-disorder-builder web page: build a disordered crystal model and run DISCUS on it.

The page exports a model (JSON, ``schema: "discuspdf-disorder-builder"``) and a short script that embeds it
and calls :func:`main`. Everything that script needs lives here:

* :mod:`.model` — reading the model, occupant choices (atoms, molecule orientations, vacancies);
* :mod:`.realize` — the configuration :class:`~discuspdf.Structure` with VOID padding, rigid-body and atomic
  disorder, and a clash report;
* :mod:`.ordering` — Monte Carlo ordering of occupants (neighbour energies, forbidden contacts);
* :mod:`.runs` — DISCUS 3D Fourier, 1D powder PDF and mmc runs through discuspdf's runners;
* :mod:`.outputs` — result files in the selected formats and the numpy preview;
* :mod:`.workflow` — :func:`run` and the command line (``python -m discuspdf.builder run model.json``).
"""
from .model import EXPORT_DEFAULTS, SCHEMA, ModelError, load_model, export_formats
from .realize import build_structure, clash_report
from .workflow import ENGINE_VERSION, cli, main, run

__all__ = ["EXPORT_DEFAULTS", "SCHEMA", "ModelError", "load_model", "export_formats", "build_structure",
           "clash_report", "ENGINE_VERSION", "cli", "main", "run"]
