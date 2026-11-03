# DISCUS stacking-fault menu: facts from the source

Source: DISCUS 6.23.01 (`/home/claude/dc_6.23.01`). File references are `discus/prog/<file>` unless
another path is given. These facts were read from the source; no DISCUS build was run (this
environment has no Fortran compiler), so behaviour not stated here is unverified.

Implemented in `src/discuspdf/discus_stack.py` (`Stacking`).

## Command order and what each command must satisfy

- Menu dispatch: `stack.f90:5` (`SUBROUTINE stack`). Branches: `average` 147-165, `layer` 360-429,
  `modulus` 430-449, `number` 450-466, `random` 467-520, `rotate` 521-711, `run` 712-716, `set`
  717-719, `show`, `sigma` 727-755, `translation` (`trans`) 756-785, `first` 320-346, `crow` 246-272,
  `ccol` 186-209, `cele` 210-233, `create` 234-245, `distribution` 273-319.
- `run` needs a prior `create`: `do_stack_fill` stops with -55 when no layers were created
  (`stack.f90:1788-1790`). `create` builds the origin list and restores the previous crystal
  (`stack.f90:1521-1590`); `run` builds the atoms (`do_stack_fill`, `stack.f90:1703-2265`).
- `create` inverts the matrix `[mod1 mod2 average]` and stops with -56 when it is singular
  (`stack.f90:1291-1295`, `lib_f90/matrix_mod.f90:262-263`). This happens with `set modulus off`
  too, so the two modulo vectors and a non-coplanar average are required.
- With `set trans` on average, a zero average translation stops `create` with -57
  (`stack.f90:1270-1273`).
- `first <k>` must satisfy `1 <= k <= number of layer types` (`stack.f90:332`, `1167`).
- `crow`/`ccol` accept a row only if its length exceeds the number of layer types read so far
  (`stack.f90:189`, `249`); the wrapper emits all layers before the matrix.

## Correlation matrix

- `crow nr,c1..cn` sets row `nr`: `st_corr(nr, c)` (`stack.f90:255`). `ccol` sets column `nr`
  (`stack.f90:195`); `cele i,j,v` sets one element (`stack.f90:220`). All three write `st_corr(row, col)`.
- The next-layer probabilities come from the row of the previous layer:
  `stack_prob` with `lprev` uses `st_corr(st_type(l), m)` (`stack.f90:1653`), called with `l = i-1`
  (`stack.f90:1333`). So `correlation[i][j]` is P(type `j` follows type `i`).
- `stack_type` draws `ran1 * ptot` with `ptot` the sum of the row (`stack.f90:1688`, `1665-1699`), so
  DISCUS normalises each row itself. The wrapper still requires rows to sum to 1 (the help text,
  `appl_dis.hlp:9639`, asks only for equal row sums, which the chain does not need).
- Random-fault weights are column sums: `prob_n(i) = sum_j st_corr(j,i)` (`stack.f90:1620-1621`),
  normalised by the total (`stack.f90:1657`).
- Average translation (used with `set trans average`): `st_mod(:,3) = sum_ij st_corr(i,j)*st_trans(i,j,:) / sum_ij st_corr(i,j)`
  (`stack.f90:1263`, `1268`). This is the per-layer mean only when every row carries equal weight;
  for a non-uniform stationary distribution it is not the true mean.

## Translations, sigma, modulus

- `trans i,j,tx,ty,tz` sets `st_trans(i,j,:)` (`stack.f90:767`). The origin of a layer of type `j`
  after type `i` adds `st_trans(i,j,:)` (`stack.f90:1346-1351`). Units: fractional coordinates
  of the layer cell (atom positions are fractional: `cr_pos = st_pos + st_origin`,
  `stack.f90:1942-1945`, `2186`).
- `sigma i,j,sx,sy,sz` sets `st_sigma(i,j,:)` (`stack.f90:738`); `gasdev(sigma)` is added to every
  component of every layer (`stack.f90:1346-1352`).
- Origins accumulate: `origin(i) = origin(i-1) + trans` (`stack.f90:1413-1415`).
- `set modulus on`: `u = st_inv * (origin - origin(1))`, then subtract `int(u1+EPS)*mod1 + int(u2+EPS)*mod2`
  (`stack.f90:1445-1451`). The origins stay in a modulo cell around the average.
- `set modulus off` (the default, `stack_mod.f90:67`): no wrapping, so random offsets accumulate.
- `aver u,v,w` sets `st_t_aver` (`stack.f90:147-165`, `155`), used as the average when
  `set trans fixed` (`st_tra_aver` false, `stack.f90:1279`).
- The scalar `set aver:<v>` sets `st_aver` (`stack.f90:903`), used only for the Fourier average
  subtraction (`stack.f90:957`, `3027`). It does not change the crystal.

## Keywords of `set trans`

- The source accepts `average` and `fixed` only (`stack.f90:889-891`, `916-918`). The help text
  (`appl_dis.hlp:9853`, `9856`, `9884`, `9885`) and the catalogue show `calculate`. `calculate` is
  not accepted: `str_comp` compares the typed prefix with the keyword (`lib_f90/str_comp_mod.f90:33-34`,
  `49-52`), so `calculate` fails the `average` test, and the `else` branch of the optional
  parameters sets `ier_num = -6` (`stack.f90:920-922`). The wrapper writes `set trans:average`.
- Optional style `set trans:average` and `set modulus:on` are parsed by `get_optional`
  (`stack.f90:850`); the positional style is `set trans,average` (`stack.f90:889`).

## Random stacking faults

- `random prob,p` sets `st_prob` (`stack.f90:479`, default `0`, `stack_mod.f90:81`).
- Per layer: `if ran1 < st_prob` the type is drawn from the column sums (`stack.f90:1320-1327`),
  otherwise from the previous row (`stack.f90:1332-1334`). The random fault then adds
  `r1*mod1 + r2*mod2 + offset + gasdev(sigma_off)` to the origin (`stack.f90:1429-1435`); `r1`, `r2`
  come from `ran1`. The offset is `random offset` (`stack.f90:492`), its spread `random sigma`
  (`stack.f90:506`).
- With matrix distribution (the default, `stack_mod.f90:60`), `ran1` is drawn for every layer after
  the first, even when `p = 0`; `first random` draws one more (`stack.f90:1152-1154`). Results depend
  on the seed only when a choice is random.

## Layer files

- `layer <file>` (`stack.f90:360-429`): the name goes through `do_build_name` (363) and
  `read_to_internal` (397; `structur.f90:8555`). A regular file goes through `import_test`
  (`structur.f90:809`) and `do_readstru_disk` (`structur.f90:815`), i.e. the `read stru` path. So a
  layer must be a DISCUS structure file (`format` line, 12 atom columns, as `save` writes). A `.cll`
  was not verified and is refused by the wrapper.
- A name starting with `internal` is DISCUS internal storage, not a file (`structur.f90:778`;
  `stack.f90:376`).
- The crystal takes the lattice of the layer read (`stack.f90:1911`, `1934`). DISCUS does not check
  that the layers share a space group or cell (`appl_dis.hlp:9714`, topic `four`).

## Crystal size and periodicity after `run`

- Atoms per layer are appended; `cr_natoms` grows by the atoms of each layer (`stack.f90:2191`).
- `cr_icc` is `estimate_ncells`: `n_i = int(max_i - min_i) + 1` over the atom coordinates
  (`stack.f90:2245-2247`; `discus_estimate.f90:57-63`). `cr_ncatoms = natoms / prod(n)`
  (`stack.f90:2248`). This is an estimate of the extent, not a count of layers.
- The crystal is not periodic: `chem_purge` true, `chem_period` false (`stack.f90:2239-2243`).

## Saving

- `save` writes an `ncell` line only with the `ncell` option (`discus_save.f90:910`;
  `sav_w_ncell` is false by default, `discus_save_mod.f90:29`). The wrapper's `save result.stru`
  therefore writes none, and `read_result` takes the cell count from the atom extents
  (see `Stacking.read_result`).
- `read_stru` needs `N` divisible by `prod(n)` and atoms in cell-major order, unless
  `check_order=False` (`discus_io.py`, `parse_stru`). A stack whose layers do not wrap into whole
  cells therefore fails; the wrapper raises `DiscusFormatError` and `StructureJob.run` raises
  `ResultNotAStructure` with the run.

## Not used by the wrapper

- Rotational disorder (`rotate`, `stack.f90:521-711`), the Fourier shortcut (`four`,
  `stack.f90:347-359`), `dist file|list` (`stack.f90:273-319`), `storage`.
