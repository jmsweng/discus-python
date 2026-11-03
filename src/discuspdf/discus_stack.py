"""DISCUS stacking-fault builder: the ``stack`` menu as a :class:`~discuspdf.discus_ops.CrystalBuilder`.

``Stacking`` builds a crystal whose layers are stacked by a Markov chain of layer types (the
correlation matrix) with translation vectors between the types, optional Gaussian spread and
random stacking faults, then decorates the layer origins with the layer atoms (``create`` then
``run``). It is the only :class:`CrystalBuilder` of the package::

    from discuspdf.discus_ops import StructureJob
    from discuspdf.discus_stack import Stacking
    from discuspdf.fourier_macro import Seed

    job = StructureJob(Stacking(layers=("cu_layer.stru", "ag_layer.stru"),
                                correlation=((0, 1), (1, 0)),
                                translations=(((0, 0, 0.5), (0, 0, 0.5)), ((0, 0, 0.5), (0, 0, 0.5))),
                                modulus=((1, 0, 0), (0, 1, 0)), number=8, first=1))
    job.macro().render()          # the DISCUS macro, checked against the command catalogue

The DISCUS facts behind each field are in ``docs/DISCUS_STACK_FACTS.md`` (file:line references
into ``discus/prog/stack.f90`` and the modules it uses). Left out on purpose: ``rotate`` (rotational
disorder), ``four`` (the Fourier shortcut), ``dist file|list`` (layer sequence read from a file),
``storage`` and ``average`` of the scalar ``set aver`` (it only changes the Fourier average
subtraction, not the crystal).

Units: translations and modulus vectors are fractional coordinates of the layer cell, which is the
cell of the saved crystal (DISCUS does not compare the cells of the layer files). The correlation
matrix row ``i`` holds the probabilities of the layer type that follows layer type ``i``.

Not verified by a run: this environment has no DISCUS build (no Fortran compiler), so the macros are
checked against the command catalogue and the source, not executed.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from pathlib import Path

from .discus_io import DiscusFormatError, read_stru
from .discus_ops import CrystalBuilder, _int, _prob, _real
from .fourier_macro import Seed
from .macro import Macro

ATOM_MODES = ("off", "on", "strict")
RANDOM_FIRST = "random"

_LAYER_FILE = re.compile(r"[A-Za-z0-9_][A-Za-z0-9_+.-]*\.stru", re.IGNORECASE)
_ROW_SUM_TOL = 1e-9  # a correlation row must sum to 1 within this
_DET_TOL = 1e-9      # |det[mod1 mod2 average]| must exceed this (stack.f90:1291-1294)
_CELL_TOL = 1e-6     # extents are read from 6-decimal text; this absorbs the rounding


def _layer_name(value) -> str:
    if not isinstance(value, str) or not _LAYER_FILE.fullmatch(value) or ".." in value:
        raise ValueError("layer must be a plain DISCUS .stru file name (no spaces, commas, quotes or "
                         f"paths), got {value!r}")
    if value.lower().startswith("internal"):
        raise ValueError(f"layer {value!r}: a name starting with 'internal' is DISCUS internal storage, "
                         "not a file (structur.f90:778)")
    return value


def _seq(name: str, value) -> tuple:
    if isinstance(value, (str, bytes)):
        raise ValueError(f"{name} must be a sequence, got {value!r}")
    try:
        return tuple(value)
    except TypeError:
        raise ValueError(f"{name} must be a sequence, got {value!r}") from None


def _vector(name: str, value, *, nonneg: bool = False, nonzero: bool = False) -> tuple:
    items = _seq(name, value)
    if len(items) != 3:
        raise ValueError(f"{name} must have three components, got {value!r}")
    out = tuple(_real(f"{name}[{k}]", x) for k, x in enumerate(items))
    if nonneg and any(x < 0.0 for x in out):
        raise ValueError(f"{name} must not have negative components, got {value!r}")
    if nonzero and all(x == 0.0 for x in out):
        raise ValueError(f"{name} must not be the zero vector")
    return out


def _square(name: str, value, n: int) -> tuple:
    rows = _seq(name, value)
    if len(rows) != n:
        raise ValueError(f"{name} must have {n} rows (one per layer type), got {len(rows)}")
    out = []
    for i, row in enumerate(rows):
        cells = _seq(f"{name}[{i}]", row)
        if len(cells) != n:
            raise ValueError(f"{name}[{i}] must have {n} entries (one per layer type), got {len(cells)}")
        out.append(cells)
    return tuple(out)


def _det3(a: tuple, b: tuple, c: tuple) -> float:
    """Determinant of the matrix whose columns are a, b, c."""
    return (a[0] * (b[1] * c[2] - b[2] * c[1])
            - b[0] * (a[1] * c[2] - a[2] * c[1])
            + c[0] * (a[1] * b[2] - a[2] * b[1]))


@dataclass(frozen=True)
class Stacking(CrystalBuilder):
    """A stacking-fault crystal from layer files (DISCUS ``stack`` menu, then ``save``).

    Fields (1-based layer types in DISCUS, 0-based here):

    * ``layers``: one ``.stru`` file per layer type, in type order (``layer <file>``). The file is
      read with the ``read stru`` path (``structur.f90:742``), so it must be a DISCUS structure
      file, typically written by ``save``; the layers must share one cell (not checked by DISCUS).
    * ``correlation``: ``correlation[i][j]`` is the probability that type ``j`` follows type ``i``.
      Rows must sum to 1 (``crow``). See the docs for the average and random-fault weights.
    * ``translations``: ``translations[i][j] = (tx, ty, tz)``, the fractional shift from the origin of
      a layer of type ``i`` to a layer of type ``j`` (``trans``). Every pair is written; a pair with
      zero probability is still used by random faults.
    * ``modulus``: the two modulo vectors ``(mod1, mod2)`` (``modulus``). Required even with
      ``modulus_on=False``: ``create`` inverts ``[mod1 mod2 average]`` (``stack.f90:1291-1294``).
    * ``number``: number of layers (``number``), at least 1.
    * ``first``: ``"random"`` (weighted by the correlation matrix, the DISCUS default) or the 1-based
      type of the first layer (``first``).
    * ``average``: ``None`` takes the average translation the DISCUS way, from the matrix
      (``set trans:average``, the equal-weight mean of the translations, ``stack.f90`` 1263-1268);
      a vector sets it with ``aver`` and ``set trans:fixed``. Must be non-zero. Modulo wrapping
      removes only whole multiples of the two modulus vectors from each origin
      (``stack.f90`` 1445-1451), so a sheared average (e.g. (1/3, 2/3, 1/3) for close packing)
      makes the in-plane origins drift with the layer number. For close packing give the
      straight stacking vector, ``average=(0, 0, 1/3)``.
    * ``modulus_on``: ``set modulus:on|off``. On, origins are wrapped into the modulo cell around the
      average, which keeps random in-plane offsets from accumulating (``stack.f90:1445-1451``).
    * ``sigmas``: ``None`` or ``sigmas[i][j] = (sx, sy, sz)``, Gaussian spread of each translation
      (``sigma``, non-negative). Pairs with all zeros are not written.
    * ``random_probability``: probability per layer of a *random-displacement* fault
      (``random prob``): the next type is drawn from the column sums of the matrix and the origin
      gets ``r1*mod1 + r2*mod2`` with r1, r2 uniform in [0, 1) (``stack.f90:1320-1336``,
      ``1429-1435``). This is not a close-packed growth or deformation fault; such sequences are
      written in the correlation matrix instead (``row i = P(next | i)``).
    * ``random_offset`` / ``random_sigma``: extra vector added at a random fault and its Gaussian
      spread (``random offset`` / ``random sigma``), both as in the DISCUS help. Written only if non-zero.
    * ``atoms``: ``None`` (DISCUS default ``off``) or ``"on"`` / ``"strict"`` (``set atoms``), shifts
      single atoms at the crystal edge into periodic positions (``appl_dis.hlp``, topic ``set``).

    Not supported: rotational disorder (``rotate``), the Fourier shortcut (``four``), and
    ``dist file|list`` (sequences read from a file).

    ``build`` writes ``seed`` (when given) before the ``stack`` menu, then ``layer``, ``number``,
    ``crow``, ``first``, ``trans``, ``sigma``, ``modulus``, ``aver``, ``set``, ``random``, ``create``
    and ``run``. ``create`` builds the origin list and ``run`` decorates it with atoms (``run``
    before ``create`` stops with -55, ``stack.f90:1788-1790``).
    """

    layers: tuple
    correlation: tuple
    translations: tuple
    modulus: tuple
    number: int
    first: object = RANDOM_FIRST
    average: object = None
    modulus_on: bool = False
    sigmas: object = None
    random_probability: float = 0.0
    random_offset: tuple = (0.0, 0.0, 0.0)
    random_sigma: tuple = (0.0, 0.0, 0.0)
    atoms: object = None

    def __post_init__(self):
        layers = tuple(_layer_name(x) for x in _seq("layers", self.layers))
        n = len(layers)
        if n < 1:
            raise ValueError("layers needs at least one layer type")
        corr = tuple(tuple(_prob(f"correlation[{i}][{j}]", x) for j, x in enumerate(row))
                     for i, row in enumerate(_square("correlation", self.correlation, n)))
        for i, row in enumerate(corr):
            if abs(sum(row) - 1.0) > _ROW_SUM_TOL:
                raise ValueError(f"correlation row {i} sums to {sum(row)!r}; a row is P(next type | "
                                 "previous type) and must sum to 1")
        raw_trans = _square("translations", self.translations, n)
        trans = tuple(tuple(_vector(f"translations[{i}][{j}]", raw_trans[i][j]) for j in range(n))
                      for i in range(n))

        mods = _seq("modulus", self.modulus)
        if len(mods) != 2:
            raise ValueError(f"modulus must be two vectors (mod1, mod2), got {len(mods)}")
        mod1 = _vector("modulus[0]", mods[0], nonzero=True)
        mod2 = _vector("modulus[1]", mods[1], nonzero=True)

        number = _int("number", self.number, 1)

        if isinstance(self.first, str):
            if self.first != RANDOM_FIRST:
                raise ValueError(f"first must be 'random' or a layer type number, got {self.first!r}")
            first = RANDOM_FIRST
        else:
            first = _int("first", self.first, 1)
            if first > n:
                raise ValueError(f"first {first} is larger than the number of layer types ({n})")

        # the average DISCUS uses (stack.f90:1263-1268): sum of correlation * translation over the
        # matrix, divided by the sum of all correlation elements
        total = sum(sum(row) for row in corr)
        computed = tuple(sum(corr[i][j] * trans[i][j][k] for i in range(n) for j in range(n)) / total
                         for k in range(3))
        if self.average is None:
            if all(x == 0.0 for x in computed):
                raise ValueError("the average translation from the matrix is zero; DISCUS stops with "
                                 "-57 (stack.f90:1270-1273). Give average=(u, v, w)")
            average = None
            used = computed
        else:
            average = _vector("average", self.average, nonzero=True)
            used = average
        if abs(_det3(mod1, mod2, used)) <= _DET_TOL:
            raise ValueError("modulus vectors and the average translation are coplanar (det 0); DISCUS "
                             "inverts [mod1 mod2 average] at create even with modulus off "
                             "(stack.f90:1291-1294, matrix_mod.f90:262-263)")

        if not isinstance(self.modulus_on, bool):
            raise ValueError(f"modulus_on must be True or False, got {self.modulus_on!r}")

        sigmas = None
        if self.sigmas is not None:
            raw = _square("sigmas", self.sigmas, n)
            sigmas = tuple(tuple(_vector(f"sigmas[{i}][{j}]", raw[i][j], nonneg=True) for j in range(n))
                           for i in range(n))

        probability = _prob("random_probability", self.random_probability)
        offset = _vector("random_offset", self.random_offset)
        r_sigma = _vector("random_sigma", self.random_sigma, nonneg=True)

        if self.atoms is not None and self.atoms not in ATOM_MODES:
            raise ValueError(f"atoms must be None or one of {ATOM_MODES}, got {self.atoms!r}")

        object.__setattr__(self, "layers", layers)
        object.__setattr__(self, "correlation", corr)
        object.__setattr__(self, "translations", trans)
        object.__setattr__(self, "modulus", (mod1, mod2))
        object.__setattr__(self, "number", number)
        object.__setattr__(self, "first", first)
        object.__setattr__(self, "average", average)
        object.__setattr__(self, "sigmas", sigmas)
        object.__setattr__(self, "random_probability", probability)
        object.__setattr__(self, "random_offset", offset)
        object.__setattr__(self, "random_sigma", r_sigma)

    @property
    def uses_random_numbers(self) -> bool:
        """True when the crystal depends on the seed: random first layer, any matrix entry strictly
        between 0 and 1, a random fault probability, a translation spread, or random offsets.
        DISCUS still draws numbers for a deterministic chain; those do not change the result."""
        if self.first == RANDOM_FIRST or self.random_probability > 0.0:
            return True
        if any(0.0 < p < 1.0 for row in self.correlation for p in row):
            return True
        if any(x != 0.0 for x in self.random_offset + self.random_sigma):
            return True
        return bool(self.sigmas) and any(x != 0.0 for row in self.sigmas for v in row for x in v)

    def input_names(self) -> tuple[str, ...]:
        return tuple(dict.fromkeys(self.layers))

    def build(self, m: Macro, seed: Seed | None) -> None:
        if seed is not None:
            m.command("seed", *seed.values)
        n = len(self.layers)
        with m.menu("stack"):
            for name in self.layers:
                m.command("layer", name)
            m.command("number", self.number)
            for i, row in enumerate(self.correlation, start=1):
                m.command("crow", i, *row)
            m.command("first", self.first)
            for i in range(1, n + 1):
                for j in range(1, n + 1):
                    m.command("trans", i, j, *self.translations[i - 1][j - 1])
                    if self.sigmas is not None and any(self.sigmas[i - 1][j - 1]):
                        m.command("sigma", i, j, *self.sigmas[i - 1][j - 1])
            m.command("modulus", *self.modulus[0], *self.modulus[1])
            if self.average is not None:
                m.command("aver", *self.average)
            m.command("set", f"trans:{'fixed' if self.average is not None else 'average'}")
            m.command("set", f"modulus:{'on' if self.modulus_on else 'off'}")
            if self.atoms is not None:
                # last on purpose: every other 'set' line resets st_mod_atom to off unless it
                # carries atoms: itself (stack.f90 921-934)
                m.command("set", f"atoms:{self.atoms}")
            if self.random_probability > 0.0:
                m.command("random", "prob", self.random_probability)
            if any(self.random_offset):
                m.command("random", "offset", *self.random_offset)
            if any(self.random_sigma):
                m.command("random", "sigma", *self.random_sigma)
            m.command("create")
            m.command("run")

    def read_result(self, path: Path, check_order: bool):
        """Read the saved stack as one periodic box (:func:`discuspdf.discus_ops.read_box`).

        ``run`` appends the atoms layer by layer (``stack.f90`` 2186-2191), so they are not in
        DISCUS's cell-major order and cannot be labelled by unit cell. DISCUS writes no ``ncell``
        line; the box is DISCUS's own estimate from the atom extents, ``floor(span) + 1`` per axis
        (``estimate_ncells``, ``stack.f90`` 2245-2247, ``discus_estimate.f90`` 57-63). With a
        sheared average stacking vector the layer origins drift in plane, and the estimated box
        grows accordingly (see ``docs/DISCUS_STACK_FACTS.md``)."""
        from .discus_ops import read_box

        return read_box(path)


__all__ = ["Stacking", "ATOM_MODES", "RANDOM_FIRST"]
