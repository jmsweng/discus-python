"""Structure editing basics: repeat, substitute, displace, check_composition, write .stru.

Run:  python examples/structure_ops_basics.py
"""
import numpy as np

from discuspdf import (
    Lattice,
    Structure,
    check_composition,
    displace,
    repeat,
    substitute,
    write_stru,
)
from discuspdf.cif_io import read_cif

CIF = "tests/fixtures/interop/gaas_cod_9008845.cif"


def main() -> None:
    # Read the COD GaAs CIF (average model, 8 atoms).
    cif = read_cif(CIF)
    s_avg = cif.structure
    print(f"Input: {s_avg!r}")
    assert s_avg.model == "average"
    assert s_avg.n_atoms == 8

    # Sample a configuration (seed=42).
    from discuspdf import sample_configuration
    s_conf = sample_configuration(s_avg, seed=42).structure
    print(f"Sampled: {s_conf!r}")
    assert s_conf.model == "configuration"
    assert s_conf.n_atoms == 8

    # Repeat 4x4x4.
    s_rep = repeat(s_conf, (4, 4, 4)).structure
    print(f"Repeated: {s_rep!r}")
    assert s_rep.n_atoms == 8 * 64  # 512
    assert tuple(s_rep.supercell) == (4, 4, 4)

    # Substitute 10% Ga -> In (seed=42).
    s_sub = substitute(s_rep, old="Ga", new="In", seed=42, fraction=0.10).structure
    n_in = int(np.sum(s_sub.species == "In"))
    n_ga = int(np.sum(s_sub.species == "Ga"))
    print(f"Substituted: {n_ga} Ga, {n_in} In")
    # 10% of 256 Ga = 25.6 -> half-up = 26
    assert n_in == 26, f"expected 26 In, got {n_in}"
    assert n_ga == 256 - 26

    # check_composition.
    check_composition(s_sub, {"Ga": 230, "In": 26, "As": 256})
    print("check_composition: OK")

    # Displace sigma=0.05 A (seed=42).
    s_disp = displace(s_sub, seed=42, sigma=0.05).structure
    print(f"Displaced: {s_disp!r}")
    assert s_disp.n_atoms == 512

    # Write .stru.
    out = "structure_ops_basics_out.stru"
    write_stru(out, s_disp)
    print(f"Wrote {out}")

    print("All asserts passed.")


if __name__ == "__main__":
    main()
