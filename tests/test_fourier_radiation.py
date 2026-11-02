"""Radiation (X-ray / neutron) and wavelength in FourierConfig.

Independent references:
* the D reference probe macros DISCUS 6.23.01 ran (benchmarks/discus_6_23_01/d_reference/inputs, written
  by its own make_inputs.py): the generated ``fourier`` and ``output`` blocks must equal them byte for
  byte, for X-ray and for neutrons;
* the NeXus files those runs wrote (…/evidence/20261101T221751Z_LAPTOP-NGU3BT0N/work): they record the
  radiation and wavelength, and ScatteringVolume.from_discus_file(config=…) must accept the matching
  configuration and refuse a mismatching one.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from discuspdf import Lattice, ReciprocalGrid, ScatteringVolume
from discuspdf.fourier_macro import WAVELENGTHS, FourierConfig, StructureSource

ROOT = Path(__file__).resolve().parents[1]
DREF = ROOT / "benchmarks" / "discus_6_23_01" / "d_reference"
WORK = DREF / "evidence" / "20261101T221751Z_LAPTOP-NGU3BT0N" / "work"
CUAU = Lattice.from_parameters(3.75, 3.75, 3.75, 90, 90, 90)
C4 = ReciprocalGrid.from_ranges(CUAU, (-2, 2, 17), (-2, 2, 17), (0, 2, 9))   # the probe's 'c4' grid


def c4_config(run_id: str, radiation: str, aver: float = 0.0, **kw) -> FourierConfig:
    return FourierConfig(C4, (("int", f"{run_id}.h5"),), aver_percent=aver, temperature="use",
                         radiation=radiation, **kw)


def probe_blocks(run_id: str) -> list[str]:
    lines = (DREF / "inputs" / f"{run_id}.mac").read_text().splitlines()
    return lines[lines.index("fourier"):-2]          # fourier … output … exit (drop the final 'exit exit')


@pytest.mark.parametrize("run_id, radiation, aver", [
    ("cuau_n_tot_c4", "neutron", 0.0), ("cuau_n_dif_c4", "neutron", 100.0),
    ("cuau_x_tot_c4", "xray", 0.0), ("cuau_x_dif_c4", "xray", 100.0),
])
def test_blocks_equal_the_probe_macros(run_id, radiation, aver):
    cfg = c4_config(run_id, radiation, aver)
    assert cfg.fourier_lines() + cfg.output_lines() == probe_blocks(run_id)


def test_default_is_xray_moa1_and_unchanged():
    cfg = FourierConfig(C4, (("int", "a.h5"),))
    assert (cfg.radiation, cfg.wavelength) == ("xray", "moa1")
    assert cfg.fourier_lines()[:3] == ["fourier", "  xray", "  wvle moa1"]


def test_radiation_line_comes_before_wvle():
    # 'neut' resets the wavelength symbol in DISCUS 6.23.01 (fourier.f90), so it must come first.
    lines = c4_config("a", "neutron", wavelength=1.54).fourier_lines()
    assert lines[1:3] == ["  neut", "  wvle 1.54"]


@pytest.mark.parametrize("given, written", [
    (1.5406, "1.5406"), (2, "2"), (0.70926, "0.70926"), (1 / 3, "0.333333333333"),
    ("CuA1", "cua1"), (" mob ", "mob"), ("wa1", "wa1"),
])
def test_wavelength_forms(given, written):
    cfg = c4_config("a", "xray", wavelength=given)
    assert cfg.fourier_lines()[2] == f"  wvle {written}"


def test_symbol_table_matches_discus():
    # per_symwl in lib_f90/element_data_mod.f90 (DISCUS 6.23.01): K-alpha1, K-alpha2, K-beta per anode.
    assert len(WAVELENGTHS) == 30 and {"moa1", "cua1", "cua2", "cub", "aga1", "wa1"} <= set(WAVELENGTHS)


def test_round_trip_and_v1_compatibility():
    cfg = c4_config("a", "neutron", wavelength=1.2)
    d = json.loads(json.dumps(cfg.to_dict()))
    assert d["version"] == 2 and d["radiation"] == "neutron" and d["wavelength"] == 1.2
    assert FourierConfig.from_dict(d) == cfg
    v1 = {k: v for k, v in c4_config("a", "xray").to_dict().items() if k != "radiation"}
    v1["version"] = 1
    assert FourierConfig.from_dict(v1) == c4_config("a", "xray")


@pytest.mark.parametrize("run_id, radiation", [("cuau_n_tot_c4", "neutron"), ("cuau_x_tot_c4", "xray")])
def test_engine_file_metadata_cross_check(run_id, radiation):
    path = WORK / f"{run_id}.h5"
    vol = ScatteringVolume.from_discus_file(path, config=c4_config(run_id, radiation))
    assert vol.provenance["fourier_config"]["radiation"] == radiation
    other = "xray" if radiation == "neutron" else "neutron"
    with pytest.raises(ValueError, match="radiation"):
        ScatteringVolume.from_discus_file(path, config=c4_config(run_id, other))
    with pytest.raises(ValueError, match="symbol"):
        ScatteringVolume.from_discus_file(path, config=c4_config(run_id, radiation, wavelength="cua1"))
    # The file records 0.70926 Å for MOA1: a numeric config is compared with that length.
    ScatteringVolume.from_discus_file(path, config=c4_config(run_id, radiation, wavelength=0.70926))
    with pytest.raises(ValueError, match="wavelength"):
        ScatteringVolume.from_discus_file(path, config=c4_config(run_id, radiation, wavelength=1.5406))


def test_render_full_macro_neutron():
    text = c4_config("job", "neutron", wavelength=1.0).render(StructureSource.stru("model.stru"))
    assert "\n  neut\n  wvle 1\n" in text and "xray" not in text
