from .pdf import PDFConfig, PDFRunner, PDFRun, read_pdf, render_pdf_macro
from .volume import (
    delta_pdf_from_diffuse,
    read_hdf5_volume,
    plot_delta_pdf_slice,
    FFTDiagnostics,
    FFTResult,
    FFTValidationError,
)

__all__ = [
    "PDFConfig", "PDFRunner", "PDFRun", "read_pdf", "render_pdf_macro",
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
