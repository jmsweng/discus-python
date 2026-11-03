from .pdf import (PDFConfig, PDFRunner, PDFRun, read_pdf, read_pdf_csv, write_pdf_csv, PDF_CSV_NAME, render_pdf_macro, QmaxLimitWarning, default_wavelength,
                  SetDeltaError, SetDeltaWarning)
from .volume import (
    delta_pdf_from_diffuse,
    read_hdf5_volume,
    plot_delta_pdf_slice,
    FFTDiagnostics,
    FFTResult,
    FFTValidationError,
)

__all__ = [
    "PDFConfig", "PDFRunner", "PDFRun", "read_pdf", "read_pdf_csv", "write_pdf_csv", "PDF_CSV_NAME", "render_pdf_macro", "QmaxLimitWarning", "default_wavelength", "SetDeltaError", "SetDeltaWarning",
    "delta_pdf_from_diffuse", "read_hdf5_volume", "plot_delta_pdf_slice",
    "FFTDiagnostics", "FFTResult", "FFTValidationError",
]

from .reciprocal import ReciprocalVolume
from .diffuse3d import Diffuse3DRunner, Diffuse3DRun, OutputSpec
__all__ += ["ReciprocalVolume", "Diffuse3DRunner", "Diffuse3DRun", "OutputSpec"]

# Version-verified DISCUS 6.23.01 common/NeXus 3D output.
from .discus_output import DiscusVolume, FractionalFFTResult, read_discus_volume, delta_pdf_fractional, evaluate_real_space
__all__ += ["DiscusVolume", "FractionalFFTResult", "read_discus_volume", "delta_pdf_fractional", "evaluate_real_space"]

# General lattice core (Package A, Feature A5).
from .lattice import Lattice
__all__ += ["Lattice"]

# Regular reciprocal-space grids (Package B, Feature 6).
from .reciprocal_grid import ReciprocalGrid, FFTCompatibility
__all__ += ["ReciprocalGrid", "FFTCompatibility"]

# Labelled scattering data, xarray adapter and versioned HDF5 (Package B, Features 7 and 19).
from .scattering_data import ScatteringVolume, mask_near_bragg
__all__ += ["ScatteringVolume", "mask_near_bragg"]

# DISCUS 6.23.01 Fourier macro writer: Bragg subtraction and LOTS (Package B, Features 8 and 9).
from .fourier_macro import FourierConfig, Lots, Seed, StructureSource, prepare_fourier_run
__all__ += ["FourierConfig", "Lots", "Seed", "StructureSource", "prepare_fourier_run"]

# Physically labelled slices, isosurfaces and VTK export (Package B, Feature 20).
from .viz import SampledVolume, plot_slice, plot_isosurfaces, isosurface, write_vtk, read_vtk
__all__ += ["SampledVolume", "plot_slice", "plot_isosurfaces", "isosurface", "write_vtk", "read_vtk"]

# Mantid MDHistoWorkspace (SaveMD v2) import/export (Package B, Feature 18).
from .mantid_io import MantidFormatError, read_mantid_mdhisto, write_mantid_mdhisto
__all__ += ["MantidFormatError", "read_mantid_mdhisto", "write_mantid_mdhisto"]

# Editable Structure model (Package A, Feature A1).
from .structure import Structure
__all__ += ["Structure"]

# DISCUS 6.23.01 .cll / .stru interchange (Package A, Feature A3).
from .discus_io import (
    DiscusFormatError,
    DiscusCell,
    DiscusStru,
    detect_discus_kind,
    read_cll,
    parse_cll,
    format_cll,
    write_cll,
    read_stru,
    parse_stru,
    format_stru,
    write_stru,
)
__all__ += [
    "DiscusFormatError", "DiscusCell", "DiscusStru", "detect_discus_kind",
    "read_cll", "parse_cll", "format_cll", "write_cll",
    "read_stru", "parse_stru", "format_stru", "write_stru",
]

# Interop helpers (Package A, Feature A2).
from .interop import (
    InteropError,
    InteropWarning,
    B_FROM_U,
    u_to_b,
    b_to_u,
    element_from_type_symbol,
    LossReport,
)
__all__ += [
    "InteropError", "InteropWarning", "B_FROM_U",
    "u_to_b", "b_to_u", "element_from_type_symbol", "LossReport",
]

# ASE adapter (Package A, Feature A2).
from .ase_io import to_ase, from_ase
__all__ += ["to_ase", "from_ase"]

# DiffPy adapter (Package A, Feature A2).
from .diffpy_io import to_diffpy, from_diffpy
__all__ += ["to_diffpy", "from_diffpy"]

# Structure editing operations (Package A, Feature A4).
from .structure_ops import EditResult, repeat, substitute, displace, check_composition, sample_configuration
__all__ += ["EditResult", "repeat", "substitute", "displace", "check_composition", "sample_configuration"]

# Neighbor lists (Package C, Feature 12).
from .neighbors import NeighborList, Shell, label_sites_by_slot, neighbor_list
__all__ += ["NeighborList", "Shell", "label_sites_by_slot", "neighbor_list"]

# Occupational and displacement correlations (Package C, Features 13 and 14).
from .correlations import (
    OccupationalCorrelation,
    DisplacementCorrelation,
    occupational_correlation,
    occupational_correlation_table,
    displacement_correlation,
    displacement_correlation_table,
    displacements,
    concentrations,
)
__all__ += [
    "OccupationalCorrelation", "DisplacementCorrelation", "occupational_correlation",
    "occupational_correlation_table", "displacement_correlation", "displacement_correlation_table",
    "displacements", "concentrations",
]

# Independent structure-factor reference (Package D, Feature 10). The functions stay in the module
# (``from discuspdf import reference_scattering as ref``) so they are not mistaken for DISCUS results.
from . import reference_scattering
from .reference_scattering import ScatteringFactors
__all__ += ["reference_scattering", "ScatteringFactors"]

# DISCUS atom names and extended XYZ output (disorder builder back end). The builder itself is the
# subpackage discuspdf.builder (imported on use: python -m discuspdf.builder run model.json).
from . import atom_names
from .xyz_io import format_extxyz, write_extxyz
__all__ += ["atom_names", "format_extxyz", "write_extxyz"]

# The whole suite: one engine layer, a checked macro builder, a general runner (suite phase 1).
from .engine import Engine, EngineChoiceError, choose_engine, fork_reasons
from .macro import Macro, MacroError
from .suite_runner import SuiteRun, SuiteRunner
__all__ += ["Engine", "EngineChoiceError", "choose_engine", "fork_reasons", "Macro", "MacroError",
            "SuiteRun", "SuiteRunner"]

# DISCUS structure operations as steps of a StructureJob (suite phase 2).
from . import discus_ops
__all__ += ["discus_ops"]
