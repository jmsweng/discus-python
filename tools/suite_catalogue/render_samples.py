"""Sample DISCUS macros from the package's three macro generators, for the suite catalogue.

``sample_macros()`` returns ``(label, macro_text)`` pairs. Together they exercise every option of
:func:`discuspdf.pdf.render_pdf_macro`, :meth:`discuspdf.fourier_macro.FourierConfig.render` and
:meth:`discuspdf.mc.MCConfig.render` that changes the generated text. Nothing here runs DISCUS.

The generators read no input files. Structure names written into the macros are basenames of committed
fixtures, and the crystal sizes and atom counts come from those files, located from the repository root.
"""
from __future__ import annotations

from pathlib import Path

from discuspdf import Lattice, PDFConfig, ReciprocalGrid, read_stru, render_pdf_macro
from discuspdf.fourier_macro import SYMMETRY, FourierConfig, Lots, Seed, StructureSource
from discuspdf.mc import (BuckinghamTarget, Convergence, DisplacementTarget, LennardJonesTarget, MCConfig, Neighborhood,
                          OccupationalTarget, RepulsiveTarget, Shift, SpringTarget, Swap, SwapDisplacement, Vector)

REPO = Path(__file__).resolve().parents[2]
FIXTURES = REPO / "tests" / "fixtures" / "discus_6_23_01"
MC_INPUTS = REPO / "benchmarks" / "discus_6_23_01" / "d_mc_probes" / "inputs"
GAAS_CLL = FIXTURES / "gaas.cll"                  # 2-atom GaAs cell, read by PDF/Fourier only by name
CUAU_STRU = MC_INPUTS / "cuau.stru"               # Cu/Au, 16x16x16 cells, one site per cell
CU_DISP_STRU = MC_INPUTS / "cu_disp.stru"         # Cu, 12x12x12 cells, one site per cell
GAAS_222_STRU = FIXTURES / "gaas_222.stru"        # GaAs, 2x2x2 cells, eight sites per cell

SEED_A = Seed((11, 22, 33, 44, 55, 66, 77, 88))
SEED_B = Seed((81, 72, 63, 54, 45, 36, 27, 18))
PAIRINGS = ("all", "local", "slocal", "site")     # the values MCConfig accepts (mc._PAIRING)


def _pdf_samples() -> list[tuple[str, str]]:
    base = dict(structure=GAAS_CLL, qmax=7.5)
    cases = [
        ("pdf_defaults", {}),                                             # xray, default wavelength, qsig 0
        ("pdf_supercell_rgrid_qmin_dq", dict(supercell=(3, 4, 5), rmax=12.0, dr=0.03, qmin=0.8, dq=0.02)),
        ("pdf_neutron_default_wavelength", dict(radiation="neutron")),
        ("pdf_xray_explicit_wavelength", dict(qmax=16.0, wavelength=0.7107)),
        ("pdf_neutron_explicit_wavelength", dict(radiation="neutron", wavelength=1.2)),
        ("pdf_qsig_allowed", dict(qsig=0.03, allow_qsig=True)),          # the only path that writes a nonzero delta
    ]
    return [(label, render_pdf_macro(PDFConfig(**{**base, **kw}))) for label, kw in cases]


def _fourier_samples() -> list[tuple[str, str]]:
    cuau = read_stru(CUAU_STRU)
    cu_disp = read_stru(CU_DISP_STRU)
    lat: Lattice = cuau.structure.lattice
    cells_cuau = tuple(cuau.ncell[:3])
    cells_disp = tuple(cu_disp.ncell[:3])

    g_hkl = ReciprocalGrid.from_ranges(lat, (-2, 2, 17), (-2, 2, 17), (-2, 2, 17))          # volume
    g_hk0 = ReciprocalGrid.from_ranges(lat, (-2, 2, 41), (-2, 2, 41), (0, 0, 1))            # plane
    g_h00 = ReciprocalGrid.from_ranges(lat, (-2, 2, 41), (0, 0, 1), (0, 0, 1))              # line
    g_hhl = ReciprocalGrid.from_corners(lat, (-2, -2, -2), (2, 2, -2), (-2, -2, 2), (-2, -2, -2), (41, 41, 1))

    s_cell = StructureSource.cell(GAAS_CLL.name, (4, 4, 4))
    s_cuau = StructureSource.stru(CUAU_STRU.name, cells_cuau)
    s_cu_disp = StructureSource.stru(CU_DISP_STRU.name, cells_disp)
    s_cuau_no_cells = StructureSource.stru(CUAU_STRU.name)

    out: list[tuple[str, str]] = []

    def add(label, grid, outputs, structure, *, title=True, **options):
        cfg = FourierConfig(grid, tuple(outputs), **options)
        out.append((label, cfg.render(structure, title=label if title else None)))

    add("fourier_hkl_int_real_imag", g_hkl,
        [("int", "hkl_int.h5"), ("real", "hkl_real.h5"), ("imag", "hkl_imag.h5")], s_cuau)
    add("fourier_hk0_plane_cell_source_no_title", g_hk0, [("int", "hk0_int.h5")], s_cell, title=False)
    add("fourier_h00_line_stru_without_cells", g_h00, [("int", "h00_int.h5")], s_cuau_no_cells)
    add("fourier_hhl_plane_auto_labels", g_hhl, [("int", "hhl_int.h5")], s_cuau, aver_percent=100)
    add("fourier_hhl_plane_explicit_labels", g_hhl, [("int", "hhl_k_l_h.h5")], s_cuau,
        aver_percent=100, labels=("k", "l", "h"))
    add("fourier_3dpdf_output_aver100", g_hkl, [("3DPDF", "hkl_3dpdf.h5")], s_cuau, aver_percent=100)
    add("fourier_all_four_outputs_aver100", g_hkl,
        [("int", "all_int.h5"), ("real", "all_real.h5"), ("imag", "all_imag.h5"), ("3DPDF", "all_3dpdf.h5")],
        s_cuau, aver_percent=100)
    add("fourier_hdf5_suffix_plane", g_hk0, [("int", "out.hdf5")], s_cell, aver_percent=100)
    add("fourier_aver50_seed_a", g_hkl, [("int", "aver50_int.h5"), ("real", "aver50_real.h5")], s_cu_disp,
        aver_percent=50, seed=SEED_A)
    add("fourier_aver12p5_seed_b", g_hkl, [("int", "aver12_int.h5")], s_cu_disp,
        aver_percent=12.5, seed=SEED_B)
    add("fourier_lots_box_all_periodic", g_hkl, [("int", "lot_box_all_y.h5")], s_cuau,
        lots=Lots("box", (4, 4, 4), "all", True))
    add("fourier_lots_box_all_open_aver100", g_hkl, [("int", "lot_box_all_n.h5")], s_cuau,
        lots=Lots("box", (4, 4, 4), "all", False), aver_percent=100)
    add("fourier_lots_box_random_seeded", g_hkl, [("int", "lot_box_64.h5")], s_cu_disp,
        lots=Lots("box", (4, 4, 4), 64, True), seed=SEED_A)
    add("fourier_lots_ellipsoid_all_periodic", g_hkl, [("int", "lot_ell_all_y.h5")], s_cuau,
        lots=Lots("ellipsoid", (4, 4, 4), "all", True))
    add("fourier_lots_ellipsoid_random_open", g_hkl, [("int", "lot_ell_9_n.h5")], s_cuau,
        lots=Lots("ellipsoid", (4, 5, 6), 9, False), seed=SEED_B)
    add("fourier_neutron_default_wavelength", g_hkl, [("int", "neut_int.h5")], s_cuau, radiation="neutron")
    add("fourier_neutron_cua1_symbol", g_hkl, [("real", "neut_cua1_real.h5")], s_cuau,
        radiation="neutron", wavelength="cua1", aver_percent=100)
    add("fourier_xray_wavelength_angstrom", g_hkl, [("real", "xray_0p7107_real.h5")], s_cell,
        radiation="xray", wavelength=0.7107)
    add("fourier_temperature_use_apply_aver100", g_hkl,
        [("real", "t_use_real.h5"), ("imag", "t_use_imag.h5")], s_cuau,
        temperature="use", symmetry="apply", aver_percent=100)
    for sym in SYMMETRY:
        add(f"fourier_symmetry_{sym}", g_hkl, [("int", f"sym_{sym}.h5")], s_cuau,
            symmetry=sym, aver_percent=100)
    # phase 2c-A
    add("fourier_electron_energy", g_hkl, [("int", "elec_int.h5")], s_cuau, radiation="electron", energy=200.0)
    add("fourier_xray_energy_delf", g_hk0, [("int", "xe_int.h5")], s_cell, energy=17.48, dispersion=True,
        delf={"CU": (0.26, 1.6), "AU": (-2.0, 8.5)})
    add("fourier_xray_internal_anomalous", g_hkl, [("int", "xa_int.h5")], s_cuau, wavelength="cua1",
        dispersion=True)
    return out


def _mc_samples() -> list[tuple[str, str]]:
    cuau = read_stru(CUAU_STRU).structure
    cu_disp = read_stru(CU_DISP_STRU).structure
    gaas = read_stru(GAAS_222_STRU, supercell=(2, 2, 2)).structure

    nn = Neighborhood.from_translations([(1, 0, 0), (0, 1, 0), (0, 0, 1)])
    out: list[tuple[str, str]] = []

    def add(label, cfg, path, structure):
        cfg.check_structure(structure)  # the sample must be valid for the structure it names
        out.append((label, cfg.render(path.name, structure.n_atoms, title=label)))

    def mc(neighborhoods, targets, moves, **options):
        return MCConfig(neighborhoods=neighborhoods, targets=targets, moves=moves,
                        seed=options.pop("seed", SEED_A), **options)

    occ_corr = OccupationalTarget(0, "Cu", "Au", correlation=-0.3)
    occ_energy = OccupationalTarget(0, "Cu", "Au", energy=0.5)

    add("mc_occ_correlation_swap_all", mc((nn,), (occ_corr,), (Swap(),)), CUAU_STRU, cuau)
    add("mc_occ_energy_swap_local_hot", mc((nn,), (occ_energy,), (Swap(0.5, "local"),),
                                           moves_per_atom=2.5, feedback_per_atom=0.5, temperature=2.0),
        CUAU_STRU, cuau)
    add("mc_occ_energy_clustering_custom_output", mc((nn,), (OccupationalTarget(0, "Cu", "Au", energy=-0.4),),
                                                     (Swap(1.0, "slocal"),), seed=Seed((1, 2, 3, 4, 5, 6, 7, 8)),
                                                     output_file="cuau_clustered.stru"),
        CUAU_STRU, cuau)
    for move_cls in (Swap, SwapDisplacement):
        for pairing in PAIRINGS:
            add(f"mc_{move_cls.__name__.lower()}_{pairing}_p050",
                mc((nn,), (OccupationalTarget(0, "Cu", "Au", energy=1.0),), (move_cls(0.5, pairing),)),
                CUAU_STRU, cuau)
    add("mc_shift_two_species_p030", mc((nn,), (OccupationalTarget(0, "Cu", "Au", energy=0.2),),
                                        (Shift({"Cu": (0.01, 0.01, 0.02), "Au": (0.01, 0.01, 0.01)}, 0.3),)),
        CUAU_STRU, cuau)
    add("mc_spring_target_shift_single_species", mc((nn,), (SpringTarget(0, "Cu", "Cu", 3.0, 1.0),),
                                                    (Shift({"Cu": (0.01, 0.01, 0.01)}),)),
        CU_DISP_STRU, cu_disp)
    add("mc_displacement_correlation", mc((nn,), (DisplacementTarget(0, "Cu", "Cu", (1, 0, 0), correlation=0.3),),
                                          (Shift({"Cu": (0.01, 0.01, 0.01)}),)),
        CU_DISP_STRU, cu_disp)
    add("mc_displacement_energy_diagonal_dir", mc((nn,), (DisplacementTarget(0, "Cu", "Cu", (1, 1, 0), energy=1.0),),
                                                  (Shift({"Cu": (0.01, 0.01, 0.01)}),)),
        CU_DISP_STRU, cu_disp)
    add("mc_lennard_jones_default_exponents", mc((nn,), (LennardJonesTarget(0, "Cu", "Cu", 2.9, 1.0),),
                                                 (Shift({"Cu": (0.01, 0.01, 0.01)}),)),
        CU_DISP_STRU, cu_disp)
    add("mc_lennard_jones_custom_exponents", mc((nn,), (LennardJonesTarget(0, "Cu", "Cu", 2.5, 0.8, 9.0, 6.0),),
                                                (Shift({"Cu": (0.01, 0.01, 0.01)}),)),
        CU_DISP_STRU, cu_disp)

    # Two-site neighbourhoods (origin and target sites differ), an asymmetric neighbourhood, one of each
    # target kind on different neighbourhoods, and all three move types together.
    nb_sites = Neighborhood((Vector(1, 2, (0, 0, 0)), Vector(1, 2, (-1, 0, 0))))
    nb_asym = Neighborhood((Vector(1, 1, (1, 0, 0)), Vector(2, 2, (0, 1, 0))), symmetric=False)
    nb_diag = Neighborhood.from_translations([(1, 1, 0), (1, 0, 1), (0, 1, 1)], site=1)
    add("mc_gaas_two_sites_all_targets_all_moves",
        mc((nb_sites, nb_asym, nb_diag),
           (OccupationalTarget(0, "GA", "AS", correlation=-0.3),
            DisplacementTarget(1, "GA", "AS", (1, 0, 0), energy=0.5),
            LennardJonesTarget(2, "GA", "GA", 2.8, 0.5)),
           (Swap(0.5, "all"), Shift({"GA": (0.01, 0.01, 0.01)}, 0.3), SwapDisplacement(0.2, "site"))),
        GAAS_222_STRU, gaas)
    # phase 2b extensions
    add("mc_repulsive_buckingham_directional_shift",
        mc((nn,), (RepulsiveTarget(0, "Cu", "Cu", 1.0, b_scale=2.5, m=12.0),
                   BuckinghamTarget(0, "Cu", "Cu", 1000.0, 0.3, 10.0)),
           (Shift({"Cu": 0.01}, directions={"Cu": (1, 0, 0)}),)),
        CU_DISP_STRU, cu_disp)
    add("mc_occ_correlation_convergence", mc((nn,), (occ_corr,), (Swap(),),
                                             finish=Convergence(diff=0.02, rdiff=0.1, change=0.05, aver=0.01)),
        CUAU_STRU, cuau)
    return out


def _ops_samples() -> list[tuple[str, str]]:
    """discuspdf.discus_ops: every step and wave variant in StructureJobs."""
    from discuspdf.discus_ops import (Purge, Remove, ReplaceAtom, ReplaceRandom, StructureJob, Switch,
                                      Thermal, Wave)

    cell = StructureSource.cell(GAAS_CLL.name, (4, 4, 4))
    stru = StructureSource.stru(GAAS_222_STRU.name, (2, 2, 2))
    jobs = {
        "ops_edit_steps": StructureJob(cell, [ReplaceRandom("GA", "IN", 0.3), ReplaceRandom("AS", "P", 0.5, site=2),
                                              ReplaceAtom(1, "AL", b=0.6), Switch(1, 2), Remove(3, 4), Purge(types=True),
                                              Remove("last"), Purge()], seed=SEED_A),
        "ops_thermal": StructureJob(stru, [Thermal(), Thermal(two_d=True), Thermal(molecules=True)], seed=SEED_B),
        "ops_wave_long_sinus": StructureJob(cell, [Wave("longitudinal", (1, 0, 0), 20.0, amplitude=0.05)]),
        "ops_wave_trans_box": StructureJob(cell, [Wave("transverse", (1, 1, 0), 15.0, function="box", box_width=0.3,
                                                       amplitude=0.1, shift=0.01, phase="random", branch="optical",
                                                       oscillation=(0, 0, 1), select=("GA",))], seed=SEED_A),
        "ops_wave_density_triangle": StructureJob(cell, [Wave("density", (0, 0, 1), 30.0, function="triangle",
                                                              replace=("GA", "IN"), viceversa=True, plow=0.1,
                                                              phigh=0.9)], seed=SEED_B),
    }
    from discuspdf.discus_chem import (ChemAnalysis, ChemDisplacement, ChemOccupational, CorrelationVector,
                                       DistanceNeighbourhood, Neighbourhood, VectorNeighbourhood)
    from discuspdf.discus_domain import Domains
    from discuspdf.discus_stack import Stacking

    jobs["ops_stack_abab"] = StructureJob(Stacking(
        layers=("a.stru", "b.stru"), correlation=((0.1, 0.9), (0.9, 0.1)),
        translations=(((0, 0, 0.5), (0, 0, 0.5)), ((0, 0, 0.5), (0, 0, 0.5))), modulus=((1, 0, 0), (0, 1, 0)),
        average=(0, 0, 0.5), number=20, first="random", modulus_on=True,
        sigmas=(((0.01, 0, 0), (0, 0, 0)), ((0, 0, 0), (0, 0.01, 0))), random_probability=0.1,
        random_offset=(0.1, 0, 0), random_sigma=(0, 0.02, 0), atoms="on"), seed=SEED_A)
    jobs["ops_domain"] = StructureJob(cell, [Domains(content="dom.stru", inputfile="pseudo.stru", name="DOM1",
                                                     character="cube",
                                                     shape=((1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0)))],
                                      seed=SEED_A, readback="box")
    jobs["ops_chem"] = StructureJob(cell, [ChemAnalysis(
        vectors=(CorrelationVector(1, 1, 2, (0, 0, 0)),),
        neighbourhoods=(Neighbourhood(VectorNeighbourhood(1), direction=((1, 0, 0), (1, 0, 0))),
                        Neighbourhood(DistanceNeighbourhood((0.25, 0.25, 0.25), 0.1))),
        calculations=(ChemOccupational("GA", "AS"), ChemDisplacement("GA", "all")))], readback="none")
    from discuspdf.discus_rmc import FourierPlane, PlaneGrid, RmcFit, RmcLots, RmcPlane, RmcShift, RmcSwap

    grid = PlaneGrid((0, -1, 0), (0, 1, 0), (3, -1, 0), 21, 31)
    jobs["ops_fourier_plane"] = StructureJob(cell, [FourierPlane(grid, "plane.nipl")], readback="none")
    jobs["ops_rmc_swap_shift"] = StructureJob(cell, [RmcFit(
        planes=(RmcPlane("plane.nipl", grid, weights="sqrt"),),
        moves=(RmcSwap(0.5, "site"), RmcShift({"GA": (0.01, 0.01, 0.01), "AS": (0.01, 0.01, 0.01)}, 0.5)),
        select=("GA", "AS"), tries=1000, sigma=0.1, fit_background=True, min_distances={("GA", "AS"): 2.0},
        debye_waller=True, symmetry=True)], seed=SEED_A)
    jobs["ops_rmc_lots_aver_fork"] = StructureJob(cell, [RmcFit(
        planes=(RmcPlane("plane.nipl", grid, radiation="neutron", wavelength=1.5),), moves=(RmcSwap(),),
        select="all", tries=1000, fit_scale=False, scale=2.0, average_percent=10.0,
        lots=RmcLots((2, 2, 2), 10, "ellipsoid", False))], seed=SEED_A)
    return [(label, job.macro().render()) for label, job in jobs.items()]


def sample_macros() -> list[tuple[str, str]]:
    """Every sample macro as ``(label, macro_text)``, PDF, Fourier, MC, then discus_ops."""
    return _pdf_samples() + _fourier_samples() + _mc_samples() + _ops_samples()


if __name__ == "__main__":
    for label, text in sample_macros():
        print(f"{label}: {len(text.splitlines())} macro lines")
