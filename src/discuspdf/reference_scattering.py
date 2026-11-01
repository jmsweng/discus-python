"""Small, independent structure-factor reference (Package D, Feature 10).

This is a *reference* for analytical fixtures and for cross-checking DISCUS — a direct sum over atoms,
written for clarity, not a competing scattering engine.  It has no lots, no Fourier grids of its own,
no anomalous dispersion, no ADP tensors and no instrument effects.

Conventions (all explicit, see ``docs/REFERENCE_SCATTERING.md``):

* Scattering vector: ``Q = 2π · (hkl @ B)`` in Å⁻¹, ``B = inv(A).T`` (rows a*, b*, c*), with ``hkl`` in
  reciprocal-lattice units of the **parent** unit cell (``cell="parent"``, the DISCUS convention) or of
  the whole supercell (``cell="supercell"``); or give Cartesian ``q`` directly.
* ``s = sin θ / λ = |Q| / 4π`` (Å⁻¹).
* Amplitude ``F(Q) = Σ_j occ_j · f_j(s) · exp(−B_j s²) · exp(+i Q·r_j)`` with ``r_j`` Cartesian (Å).
  The ``+i`` sign is DISCUS's.  The Debye–Waller factor is applied only with ``debye_waller=True``.
* Units: X-ray ``f`` in electrons → ``|F|²`` in electrons²; neutron ``b`` in fm → ``|F|²`` in fm².
  ``normalize="none"`` is the whole sample; ``"per_atom"`` divides by Σ occ; ``"per_cell"`` by the
  number of parent cells.

Scattering factors are never guessed: they come from a :class:`ScatteringFactors` table that names its
source (user values, ``xraydb`` Waasmaier–Kirfel X-ray factors, or ``periodictable`` neutron lengths).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable, Mapping

import numpy as np

from ._elements import _ELEMENTS, _strip_suffix
from .structure import Structure

__all__ = [
    "ScatteringFactors",
    "scattering_vectors",
    "structure_factor",
    "intensity",
    "average_structure_factor",
    "diffuse_intensity",
]

_RADIATIONS = ("xray", "neutron")
_NORMALIZE = ("none", "per_atom", "per_cell")
_MAX_CHUNK_ELEMENTS = 4_000_000   # Q points × atoms held in one phase block


# --------------------------------------------------------------------------------------------- factors
def _element_symbol(label: str) -> str:
    """Canonical element symbol for a plain element label ('CU', 'cu', 'Cu1' -> 'Cu').

    Charged labels ('O2-', 'Fe3+') are refused: their factors must be given explicitly."""
    if not isinstance(label, str) or not label:
        raise ValueError(f"species label must be a non-empty string, got {label!r}")
    if label.rstrip().endswith(("+", "-")):
        raise ValueError(f"species {label!r} is an ion; give its scattering factor explicitly")
    stripped = _strip_suffix(label)
    canonical = stripped[0].upper() + stripped[1:].lower()
    if canonical not in _ELEMENTS:
        raise ValueError(f"{label!r} is not an element symbol; give its scattering factor explicitly")
    return canonical


def _wk_function(coef: Mapping[str, Any]) -> Callable[[np.ndarray], np.ndarray]:
    a = np.asarray(coef["a"], float)
    b = np.asarray(coef["b"], float)
    c = float(coef["c"])
    if a.shape != b.shape or a.ndim != 1 or a.size == 0:
        raise ValueError("Waasmaier-Kirfel coefficients need equal-length 'a' and 'b' lists")
    if not (np.isfinite(a).all() and np.isfinite(b).all() and np.isfinite(c)):
        raise ValueError("Waasmaier-Kirfel coefficients must be finite")

    def f(s: np.ndarray) -> np.ndarray:
        s2 = np.asarray(s, float) ** 2
        return c + np.sum(a[:, None] * np.exp(-b[:, None] * s2.ravel()[None, :]), axis=0).reshape(s2.shape)

    return f


@dataclass(frozen=True)
class ScatteringFactors:
    """Scattering factor per species label, with a stated radiation and source.

    Keys are matched to ``Structure.species`` exactly, then case-insensitively.  Use the constructors:
    :meth:`neutron`, :meth:`xray_wk`, :meth:`xray_functions`, :meth:`from_xraydb`,
    :meth:`from_periodictable`.  ``factor(species, s)`` returns f(s) (X-ray, electrons) or the constant
    coherent length b (neutron, fm) broadcast to the shape of ``s``.
    """
    radiation: str
    table: Mapping[str, Any]
    source: str
    _lookup: dict = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        if self.radiation not in _RADIATIONS:
            raise ValueError(f"radiation must be one of {_RADIATIONS}")
        if not isinstance(self.source, str) or not self.source:
            raise ValueError("state the source of the scattering factors")
        if not self.table:
            raise ValueError("no scattering factors given")
        lookup: dict[str, Any] = {}
        for k, v in self.table.items():
            if not isinstance(k, str) or not k:
                raise ValueError("species keys must be non-empty strings")
            if self.radiation == "neutron":
                if isinstance(v, (complex, np.complexfloating)) or not np.isfinite(float(v)):
                    raise ValueError(f"neutron length for {k!r} must be a finite real number (fm)")
                v = float(v)
            elif not callable(v):
                raise ValueError(f"X-ray factor for {k!r} must be a function of s")
            if k.upper() in lookup:
                raise ValueError(f"species {k!r} given twice (keys are case-insensitive)")
            lookup[k.upper()] = v
        object.__setattr__(self, "table", dict(self.table))
        object.__setattr__(self, "_lookup", lookup)

    # -- constructors ----------------------------------------------------------------------------
    @classmethod
    def neutron(cls, lengths: Mapping[str, float], source: str = "user") -> "ScatteringFactors":
        """Coherent neutron scattering lengths in fm, e.g. ``{"Cu": 7.718, "Au": 7.63}``."""
        return cls("neutron", dict(lengths), source)

    @classmethod
    def xray_wk(cls, coefficients: Mapping[str, Mapping[str, Any]], source: str = "user") -> "ScatteringFactors":
        """X-ray factors from Waasmaier–Kirfel style coefficients ``{"a": [...], "b": [...], "c": c}``:
        ``f(s) = c + Σ a_i exp(−b_i s²)``."""
        return cls("xray", {k: _wk_function(v) for k, v in coefficients.items()}, source)

    @classmethod
    def xray_functions(cls, functions: Mapping[str, Callable], source: str = "user") -> "ScatteringFactors":
        """X-ray factors as arbitrary functions ``f(s)`` (s in Å⁻¹, f in electrons)."""
        return cls("xray", dict(functions), source)

    @classmethod
    def from_xraydb(cls, species) -> "ScatteringFactors":
        """Neutral-atom Waasmaier–Kirfel factors from the ``xraydb`` package (``pip install xraydb``)."""
        try:
            import xraydb
        except ImportError as exc:                                     # pragma: no cover - env specific
            raise ImportError("from_xraydb needs the optional package 'xraydb'") from exc
        funcs = {}
        for label in sorted(set(species)):
            sym = _element_symbol(label)
            funcs[label] = (lambda s, _sym=sym: np.asarray(xraydb.f0(_sym, np.asarray(s, float).ravel()),
                                                           float).reshape(np.shape(s)))
        return cls("xray", funcs, f"xraydb {getattr(xraydb, '__version__', '?')} f0 (Waasmaier-Kirfel 1995)")

    @classmethod
    def from_periodictable(cls, species) -> "ScatteringFactors":
        """Bound coherent lengths ``b_c`` (fm) from ``periodictable`` (natural isotopic abundance).

        Elements whose length has an imaginary (absorption) part are refused: give them explicitly."""
        try:
            import periodictable
        except ImportError as exc:                                     # pragma: no cover - env specific
            raise ImportError("from_periodictable needs the optional package 'periodictable'") from exc
        lengths = {}
        for label in sorted(set(species)):
            el = getattr(periodictable, _element_symbol(label))
            b, bi = el.neutron.b_c, el.neutron.b_c_i
            if b is None:
                raise ValueError(f"periodictable has no coherent length for {label!r}")
            if bi not in (None, 0, 0.0):
                raise ValueError(f"{label!r} has a complex coherent length (b_c_i = {bi}); give it explicitly")
            lengths[label] = float(b)
        return cls("neutron", lengths, f"periodictable {periodictable.__version__} neutron.b_c")

    # -- evaluation -------------------------------------------------------------------------------
    def has(self, species: str) -> bool:
        return species in self.table or species.upper() in self._lookup

    def factor(self, species: str, s) -> np.ndarray:
        s = np.asarray(s, float)
        if species in self.table:
            v = self.table[species]
        elif species.upper() in self._lookup:
            v = self._lookup[species.upper()]
        else:
            raise KeyError(f"no {self.radiation} scattering factor for species {species!r} ({self.source})")
        if self.radiation == "neutron":
            return np.full(s.shape, float(v))
        out = np.asarray(v(s), float)
        if out.shape != s.shape or not np.isfinite(out).all():
            raise ValueError(f"X-ray factor for {species!r} returned a bad array (shape {out.shape})")
        return out

    def describe(self) -> dict:
        return {"radiation": self.radiation, "source": self.source, "species": sorted(self.table)}


# --------------------------------------------------------------------------------------------- geometry
def _parent_matrix(structure: Structure) -> np.ndarray:
    return structure.lattice.matrix / structure.supercell[:, None]


def scattering_vectors(structure: Structure, *, hkl=None, q=None, cell: str = "parent") -> np.ndarray:
    """Cartesian Q (Å⁻¹, 2π included) for ``hkl`` in r.l.u. of the parent cell or supercell, or ``q``."""
    if (hkl is None) == (q is None):
        raise ValueError("give exactly one of hkl= or q=")
    if q is not None:
        Q = np.asarray(q, float)
        if Q.shape[-1:] != (3,) or not np.isfinite(Q).all():
            raise ValueError("q must be finite with last dimension 3")
        return Q
    if cell not in ("parent", "supercell"):
        raise ValueError("cell must be 'parent' or 'supercell'")
    h = np.asarray(hkl, float)
    if h.shape[-1:] != (3,) or not np.isfinite(h).all():
        raise ValueError("hkl must be finite with last dimension 3")
    A = _parent_matrix(structure) if cell == "parent" else structure.lattice.matrix
    B = np.linalg.inv(A).T
    return 2.0 * np.pi * (h @ B)


def _check_inputs(structure: Structure, factors: ScatteringFactors, debye_waller: bool) -> None:
    if not isinstance(structure, Structure):
        raise TypeError("structure must be a Structure")
    if not isinstance(factors, ScatteringFactors):
        raise TypeError("factors must be ScatteringFactors")
    missing = sorted({sp for sp in structure.species if not factors.has(sp)})
    if missing:
        raise KeyError(f"no {factors.radiation} scattering factor for {missing} ({factors.source})")
    if debye_waller and structure.b_iso is None:
        raise ValueError("debye_waller=True but the structure has no b_iso values")


def _amplitude(Qflat: np.ndarray, r: np.ndarray, weight_fn, species: np.ndarray, chunk: int | None) -> np.ndarray:
    """Σ_j w_j(s) exp(i Q·r_j) for a flat (M, 3) Q array; w per species × atom-scalar."""
    M, N = len(Qflat), len(r)
    out = np.zeros(M, complex)
    if N == 0 or M == 0:
        return out
    step = chunk or max(1, _MAX_CHUNK_ELEMENTS // max(N, 1))
    groups = [(sp, np.flatnonzero(species == sp)) for sp in sorted(set(species))]
    for i0 in range(0, M, step):
        Qc = Qflat[i0:i0 + step]
        s = np.linalg.norm(Qc, axis=1) / (4.0 * np.pi)
        acc = np.zeros(len(Qc), complex)
        for sp, idx in groups:
            w = weight_fn(sp, s, idx)                        # (m, n_sp) real
            acc += np.einsum("mn,mn->m", w, np.exp(1j * (Qc @ r[idx].T)))
        out[i0:i0 + step] = acc
    return out


def _weights(structure: Structure, factors: ScatteringFactors, debye_waller: bool):
    occ = structure.occupancy
    biso = structure.b_iso if debye_waller else None

    def weight(sp, s, idx):
        f = factors.factor(sp, s)[:, None]                   # (m, 1)
        w = f * occ[idx][None, :]
        if biso is not None:
            w = w * np.exp(-np.outer(s ** 2, biso[idx]))
        return w

    return weight


def structure_factor(structure: Structure, factors: ScatteringFactors, *, hkl=None, q=None,
                     cell: str = "parent", debye_waller: bool = False, chunk: int | None = None) -> np.ndarray:
    """Complex amplitude ``F(Q) = Σ_j occ_j f_j(s) [exp(−B_j s²)] exp(i Q·r_j)``; shape of ``hkl[..., 0]``."""
    _check_inputs(structure, factors, debye_waller)
    Q = scattering_vectors(structure, hkl=hkl, q=q, cell=cell)
    shape = Q.shape[:-1]
    F = _amplitude(Q.reshape(-1, 3), structure.cartesian, _weights(structure, factors, debye_waller),
                   np.asarray(structure.species, dtype=object), chunk)
    return F.reshape(shape)


def _normalizer(structure: Structure, normalize: str) -> float:
    if normalize not in _NORMALIZE:
        raise ValueError(f"normalize must be one of {_NORMALIZE}")
    if normalize == "per_atom":
        n = float(structure.occupancy.sum())
        if n <= 0:
            raise ValueError("per_atom normalization needs a nonzero total occupancy")
        return n
    if normalize == "per_cell":
        return float(np.prod(structure.supercell))
    return 1.0


def intensity(structure: Structure, factors: ScatteringFactors, *, hkl=None, q=None, cell: str = "parent",
              debye_waller: bool = False, normalize: str = "none", chunk: int | None = None) -> np.ndarray:
    """``|F(Q)|²`` (electrons² or fm²) for the whole sample, or per atom / per parent cell."""
    norm = _normalizer(structure, normalize)
    F = structure_factor(structure, factors, hkl=hkl, q=q, cell=cell, debye_waller=debye_waller, chunk=chunk)
    return np.abs(F) ** 2 / norm


def average_structure_factor(structure: Structure, factors: ScatteringFactors, *, hkl=None, q=None,
                             cell: str = "parent", debye_waller: bool = False,
                             chunk: int | None = None) -> np.ndarray:
    """Amplitude of the average periodic structure, ``⟨F⟩ = L(Q) · Σ_sites ⟨f e^{−Bs²} e^{iQ·(r−R_c)}⟩_cells``.

    ``R_c = cell_index @ A_parent`` is the origin of the atom's parent cell and ``L(Q) = Σ_c e^{iQ·R_c}``
    runs over all ``prod(supercell)`` cells; a site missing in a cell (vacancy) contributes zero to the
    cell average.  This is the full (100 %) average of DISCUS ``set aver``; it needs a periodic
    configuration model."""
    _check_inputs(structure, factors, debye_waller)
    if structure.model != "configuration":
        raise ValueError("the average structure needs a configuration model (one realised crystal)")
    if not all(structure.pbc):
        raise ValueError("the average structure needs periodic boundaries in all three directions")
    Q = scattering_vectors(structure, hkl=hkl, q=q, cell=cell)
    shape = Q.shape[:-1]
    Qf = Q.reshape(-1, 3)
    n = structure.supercell
    ncell = int(np.prod(n))
    Ap = _parent_matrix(structure)
    R = structure.cell_index @ Ap
    local = structure.cartesian - R
    weight = _weights(structure, factors, debye_waller)
    # cell-averaged amplitude of one cell: (1/Ncell) Σ_j w_j exp(iQ·(r_j − R_cj))
    Fcell = _amplitude(Qf, local, weight, np.asarray(structure.species, dtype=object), chunk) / ncell
    # lattice sum over all cells, separable in the parent basis: Π_k Σ_{m=0}^{n_k−1} exp(i m Q·a_k)
    L = np.ones(len(Qf), complex)
    for k in range(3):
        phase = Qf @ Ap[k]
        m = np.arange(n[k])
        L *= np.exp(1j * np.outer(phase, m)).sum(axis=1)
    return (L * Fcell).reshape(shape)


def diffuse_intensity(structure: Structure, factors: ScatteringFactors, *, hkl=None, q=None,
                      cell: str = "parent", debye_waller: bool = False, normalize: str = "none",
                      chunk: int | None = None) -> np.ndarray:
    """``|F − ⟨F⟩|²`` — the amplitude subtraction that DISCUS ``set aver,100`` performs."""
    norm = _normalizer(structure, normalize)
    kw = dict(hkl=hkl, q=q, cell=cell, debye_waller=debye_waller, chunk=chunk)
    F = structure_factor(structure, factors, **kw)
    Fav = average_structure_factor(structure, factors, **kw)
    return np.abs(F - Fav) ** 2 / norm
