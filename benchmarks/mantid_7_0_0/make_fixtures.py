"""Generate Mantid MDHisto fixtures for discuspdf's reader (run with Mantid's Python, not discuspdf's).

    <mantid-python> benchmarks/mantid_7_0_0/make_fixtures.py tests/fixtures/mantid_7_0_0

Every number the tests check is decided here, before Mantid writes anything, and stored in
expected.json: event positions (always at bin centres), event weights, the CreateMDHistoWorkspace
signal formula and the masked box. Mantid version is recorded.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import mantid
from mantid.simpleapi import (AddSampleLog, BinMD, ConvertToMD, CreateMDHistoWorkspace, CreateMDWorkspace,
                              CreateSampleWorkspace, FakeMDEventData, MaskMD, SaveMD, SetUB, mtd)

CELL = dict(a=5.0, b=6.0, c=7.0, alpha=90.0, beta=105.0, gamma=90.0)


def md_event_ws(name):
    CreateMDWorkspace(Dimensions=3, Extents="-3,3,-3,3,-3,3", Names="[H,0,0],[0,K,0],[0,0,L]",
                      Units="r.l.u.,r.l.u.,r.l.u.", Frames="HKL,HKL,HKL", SplitInto=2, MaxRecursionDepth=1,
                      OutputWorkspace=name)
    AddSampleLog(name, LogName="discuspdf_fixture", LogText=name)   # creates the experiment info
    SetUB(name, **CELL)


def add_events(name, events):
    for n, x, y, z in events:
        FakeMDEventData(name, PeakParams=f"{n},{x},{y},{z},1e-7")


def main(out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    expected = {"mantid_version": mantid.__version__, "cell": CELL, "files": {}}

    # 1. axis-aligned BinMD: bins of 0.25 (H, K) and 0.5 (L) from -3; centres at -2.875 + 0.25 i etc.
    md_event_ws("ev")
    ev_aligned = [(100, 1.125, -0.375, 0.25), (300, -1.875, 1.125, -0.75)]
    add_events("ev", ev_aligned)
    BinMD("ev", AlignedDim0="[H,0,0],-3,3,24", AlignedDim1="[0,K,0],-3,3,24", AlignedDim2="[0,0,L],-3,3,12",
          OutputWorkspace="aligned")
    SaveMD("aligned", Filename=str(out / "aligned_hkl.nxs"))
    expected["files"]["aligned_hkl.nxs"] = {"events_hkl": [[x, y, z, n] for n, x, y, z in ev_aligned],
                                            "shape": [24, 24, 12]}

    # 2. projected BinMD: basis [H,H,0], [-H,H,0], [0,0,L]; events at projected bin centres.
    md_event_ws("ev2")
    proj_pts = [(200, 0.625, -1.125, 0.25), (50, -1.875, 0.375, -2.25)]   # projected coordinates (x, y, z)
    U, V, Wv = np.array([1, 1, 0.]), np.array([-1, 1, 0.]), np.array([0, 0, 1.])
    hkl_events = [(n, *(x * U + y * V + z * Wv)) for n, x, y, z in proj_pts]
    add_events("ev2", hkl_events)
    BinMD("ev2", AxisAligned=False, BasisVector0="[H,H,0],r.l.u.,1,1,0", BasisVector1="[-H,H,0],r.l.u.,-1,1,0",
          BasisVector2="[0,0,L],r.l.u.,0,0,1", OutputExtents="-3,3,-3,3,-3,3", OutputBins="24,24,12",
          NormalizeBasisVectors=False, OutputWorkspace="projected")
    SaveMD("projected", Filename=str(out / "projected_hh0.nxs"))
    expected["files"]["projected_hh0.nxs"] = {"events_hkl": [[h, k, l, n] for n, h, k, l in hkl_events],
                                              "events_projected": [[x, y, z] for _, x, y, z in proj_pts],
                                              "shape": [24, 24, 12]}

    # 3. ConvertToMD with Uproj/Vproj/Wproj (writes W_MATRIX), then an aligned BinMD; one fake event
    #    added between two otherwise identical binnings, so the difference is that event alone.
    m = CreateSampleWorkspace(OutputWorkspace="m", WorkspaceType="Event", NumBanks=1, BankPixelWidth=4,
                              XUnit="TOF", XMin=1000, XMax=20000)
    SetUB("m", **CELL)
    ConvertToMD("m", QDimensions="Q3D", Q3DFrames="HKL", dEAnalysisMode="Elastic", QConversionScales="HKL",
                Uproj="1,1,0", Vproj="0,0,1", Wproj="1,-1,0", MinValues="-5,-5,-5", MaxValues="5,5,5",
                OutputWorkspace="cmd")
    names = [mtd["cmd"].getDimension(i).name for i in range(3)]
    binning = dict(AlignedDim0=f"{names[0]},-5,5,20", AlignedDim1=f"{names[1]},-5,5,20",
                   AlignedDim2=f"{names[2]},-5,5,20")
    BinMD("cmd", OutputWorkspace="cmd_before", **binning)
    xyz = (1.25, -0.75, 2.25)                       # bin centres of 0.5-wide bins from -5
    FakeMDEventData("cmd", PeakParams=f"500,{xyz[0]},{xyz[1]},{xyz[2]},1e-7")
    BinMD("cmd", OutputWorkspace="cmd_after", **binning)
    SaveMD("cmd_before", Filename=str(out / "w_matrix_before.nxs"))
    SaveMD("cmd_after", Filename=str(out / "w_matrix_after.nxs"))
    W = np.array([[1, 0, 1], [1, 0, -1], [0, 1, 0]], float)      # columns u, v, w as requested above
    expected["files"]["w_matrix_after.nxs"] = {"event_projected": list(xyz), "event_hkl": (W @ np.array(xyz)).tolist(),
                                               "weight": 500, "names": names, "difference_with": "w_matrix_before.nxs",
                                               "shape": [20, 20, 20]}

    # 4. CreateMDHistoWorkspace with a known signal formula (Mantid's input order: D0 fastest),
    #    a masked box and one NaN.
    n0, n1, n2 = 6, 5, 4
    i0, i1, i2 = np.meshgrid(np.arange(n0), np.arange(n1), np.arange(n2), indexing="ij")
    signal = 100.0 * i0 + 10.0 * i1 + i2 + 1.0
    error = 0.1 * signal
    flat = lambda a: a.ravel(order="F")                         # D0 fastest
    signal_in = flat(signal).copy()
    nan_index = (2, 3, 1)
    signal_in[np.ravel_multi_index(nan_index, (n0, n1, n2), order="F")] = np.nan
    CreateMDHistoWorkspace(SignalInput=",".join(repr(float(v)) for v in signal_in),
                           ErrorInput=",".join(repr(float(v)) for v in flat(error)), Dimensionality=3,
                           Extents="-1.5,1.5,-1.25,1.25,0,2", NumberOfBins=f"{n0},{n1},{n2}",
                           Names="[H,0,0],[0,K,0],[0,0,L]", Units="r.l.u.,r.l.u.,r.l.u.", Frames="HKL,HKL,HKL",
                           OutputWorkspace="histo")
    AddSampleLog("histo", LogName="discuspdf_fixture", LogText="histo")
    SetUB("histo", **CELL)
    MaskMD("histo", Dimensions="[H,0,0],[0,K,0],[0,0,L]", Extents="-1.5,-1.0,-1.25,1.25,0,2")   # first H bin
    SaveMD("histo", Filename=str(out / "histo_formula.nxs"))
    expected["files"]["histo_formula.nxs"] = {
        "signal_formula": "100*i0 + 10*i1 + i2 + 1 (i0 along [H,0,0])", "error_formula": "0.1*signal",
        "nan_index": list(nan_index), "masked": "i0 == 0", "shape": [n0, n1, n2],
        "extents": [[-1.5, 1.5], [-1.25, 1.25], [0, 2]]}

    (out / "expected.json").write_text(json.dumps(expected, indent=1) + "\n", encoding="utf-8")
    print("wrote", sorted(p.name for p in out.iterdir()))


if __name__ == "__main__":
    main(Path(sys.argv[1]).resolve())
