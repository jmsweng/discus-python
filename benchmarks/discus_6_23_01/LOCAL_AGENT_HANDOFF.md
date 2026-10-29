# Handoff to the local agent: run the pre-registered DISCUS experiments (Gate 2 and Gate 1)

**From:** Claude (independent reviewer), 2026-10-28, on behalf of the project owner (jmsweng).
**Repo / branch:** `jmsweng/discus-python`, branch `science/discus-3d-dpdf-benchmarks` (PR #8, **draft**).
**Why you:** the DISCUS_SUITE 6.23.01 engine exists only on the owner's machine: WSL Ubuntu-22.04, `/usr/local/bin/discus_suite`.

## Your job in one paragraph

Run six pre-registered DISCUS macros exactly as committed, in fresh evidence directories. Run the two committed analysis scripts. Commit the raw evidence (structures, macros, logs, HDF5, analysis output, hash manifests) to the PR branch, push, and report the results verbatim. **You are an operator, not an analyst.**
- Do not change any registered file.
- Do not tune anything.
- Do not interpret a failure as a convention problem.
- Do not edit gate statuses or reports.

If anything deviates from this document, stop and report.

## Hard rules (owner's standing instructions)

1. **No merge** of PR #8. Do not mark it ready for review. Do not force-push or rewrite history.
2. **No deletion** of any file or directory, anywhere. Never use `spike_discus_3d/run_cases.sh` (its `--overwrite` has deleted directories before). Never `rm -rf`. Never `git clean`, `git reset --hard` or `git checkout -- <file>` on files you did not create.
3. **Do not modify** anything under:
   - `benchmarks/discus_6_23_01/gate2_prereg/` except creating new `evidence/<RUN_ID>/` directories;
   - `benchmarks/discus_6_23_01/gate1_c1_dyadic/` except creating new `evidence/<RUN_ID>/` directories;
   - any document or test anywhere.

   `analyze.py` checks SHA-256 hashes and will declare the run INVALID if registered files change.
4. **Do not reinstall or rebuild DISCUS**, and do not touch other WSL distributions. Ask the owner before any installation. The only exception is a missing Python package in the project venv: `pip install h5py numpy pytest`, nothing else.
5. **Each attempt gets a new evidence directory.** Never write into an existing one. If a run fails, keep everything, start a new directory, and record why.

## Step 0: make a fresh, LF-only clone inside the WSL filesystem

**Revised 2026-10-28 by the reviewer, after the first attempt correctly stopped here.** The existing clone at `/mnt/c/Users/jmswe/Documents/Programming projects/discus-python` shows ~53 files as modified, and they differ only in line endings. git stores every text file with LF; the Windows-side checkout wrote CRLF. If the registered inputs arrived the same way, their SHA-256 hashes would not match and DISCUS might misread the macros.

So: **do not use, change, clean, pull or commit in the `/mnt/c` clone.** Leave its modified and untracked files exactly as they are; the owner will deal with them separately. Make a new clone on the Linux filesystem, with no line-ending conversion:

```bash
WORK=~/gate2_prereg_runs/discus-python
test ! -e "$WORK" || { echo "STOP: $WORK already exists"; exit 1; }
mkdir -p ~/gate2_prereg_runs
git -c core.autocrlf=false clone --branch science/discus-3d-dpdf-benchmarks \
    https://github.com/jmsweng/discus-python "$WORK"     # if authentication fails: STOP and report
cd "$WORK"
git config core.autocrlf false        # this clone only
git log -1 --format='%H %s'           # record (REGISTERED_COMMIT); must include the gate2-prereg-v2 commit (see 'Registration v2' below)
git status --short                    # must print nothing

# Line-ending check: every registered file must be LF in git AND in the working tree.
git ls-files --eol benchmarks/discus_6_23_01/gate2_prereg benchmarks/discus_6_23_01/gate1_c1_dyadic tests/test_gate2_prereg.py \
  | awk '$1 != "i/-text" && ($1 != "i/lf" || $2 != "w/lf")'   # must print nothing
ls benchmarks/discus_6_23_01/gate2_prereg/inputs benchmarks/discus_6_23_01/gate1_c1_dyadic/inputs
```

Use `"$WORK"` wherever this document says `<REPO>`.

**Python environment for this clone.** Create a fresh venv here; installing only these three packages is approved:

```bash
python3 -m venv .venv && source .venv/bin/activate && pip install numpy h5py pytest
```

Nothing needs to be installed into the package itself, because the scripts add `src/` to the path.

Expected inputs:
- **Gate 2:** `gate2_block.stru`, `gate2_random.stru`, `gate2_{block,random}_{total,diffuse}.mac`, `gate2_probe.mac`.
- **Gate 1:** `C1_frozen.stru`, `c1_dyadic.mac`, `c1_repro.mac`.

## Gate 3 run (2026-10-28): READ THIS FIRST for the Gate 3 request

Gate 2 is finished. Do **not** rerun Gate 1 or Gate 2. Gate 3 is a short pre-registered run: two macros, about 1 minute in total.

1. `git pull --ff-only`. The head must include the commit "Gate 3: pre-register native 3DPDF experiment". Then:
   - `python -m pytest -q tests/test_gate2_prereg.py tests/test_gate3_prereg.py tests/test_regenerate_evidence.py` → `15 passed`.
   - `python benchmarks/discus_6_23_01/gate3/gate3_selftest.py` → 8 lines starting `ok`.
   - `sha256sum benchmarks/discus_6_23_01/gate3/predictions.json` → `e6368b229db439b40b2e81829875c4cf1a8fe34c6ba6f9b3d53848bb004493e6`.
2. Engine identity: same as Step 2 (binary hash `d4274f13…45706f3`). Write ENGINE_ENV.txt into the new evidence directory.
3. Run:
   ```sh
   E3="benchmarks/discus_6_23_01/gate3/evidence/$RUN_ID"; mkdir -p "$E3"
   cp -n benchmarks/discus_6_23_01/gate3/inputs/* "$E3/"; cd "$E3"
   for m in gate3_block_total gate3_block_diffuse; do
     timeout 3600 /usr/local/bin/discus_suite < "$m.mac" > "$m.log" 2>&1; echo "$m exit=$?" | tee -a RUN_STATUS.txt
   done
   cd "<REPO>"
   python benchmarks/discus_6_23_01/gate3/gate3_analyze.py "$E3"; echo "analyze exit=$?" | tee -a "$E3/RUN_STATUS.txt"
   python benchmarks/discus_6_23_01/regenerate_evidence.py hash "$E3"
   ( cd "$E3" && sha256sum * > MANIFEST.sha256 )
   ```
   Stop conditions:
   - any `***APPL***` line in a log;
   - a missing `.h5` file;
   - `Temp. factors` not `ignored`.

   A FAIL from `analyze.py` is not a stop condition.
4. Commit only `$E3`. The `.h5` files are git-ignored; that is intended. Push, then report:
   - `GATE3_RESULTS_PREREG.md` in full;
   - the `hash` output;
   - RUN_STATUS.txt;
   - any deviations.

## Evidence volumes (2026-10-28)

Engine `.h5` outputs in `evidence/` directories are git-ignored and are no longer committed. After a run, add `python benchmarks/discus_6_23_01/regenerate_evidence.py hash <evidence_dir>` output to your report so the reviewer can record the content hashes. Everything else in the evidence directory is committed as before. See `EVIDENCE_VOLUMES.md`.

## Registration v4 (2026-10-28): READ THIS FIRST. It overrides everything below where it conflicts, including the v3 section.

v3 failed one criterion because of rounding ties (`gate2_prereg/RESULT_v3.md`). The owner approved **gate2-prereg-v4**, which makes two changes: tie-aware criteria (new `R` rows) and **fresh structure seeds**. The macros are unchanged, but they now read new `.stru` files.

- `git pull --ff-only`; the head must include the commit titled "Gate 2: register v4". Expected: `pytest tests/test_gate2_prereg.py` → `9 passed`; `python benchmarks/discus_6_23_01/gate2_prereg/selftest.py` → 12 lines starting `ok`; `gate2_prereg/predictions.json` SHA-256 = `bebdc86a80817439091834bb1252f5bb63f2f19644930c3bd298c4c87acba311`.
- **Do NOT rerun Gate 1.** Only Gate 2 runs, with a new `RUN_ID`.
- Order and stop conditions: same as v3. Run Step 1, then Step 2 (the binary hash must be `d4274f13…45706f3`), then the probe, which must give 4096 / 2 / 1. Then run `gate2_block_total.mac` first: its log must contain `Temp. factors      : ignored` and no `***APPL***` line. Then the other three macros, then `analyze.py`. A FAIL is not a stop condition.
- Commit only the new v4 evidence directory. Do not edit any earlier evidence. Report as before; also paste the `tie_rounding` block from `gate2_results.json` for each model.

## Registration v3 (2026-10-28): superseded by v4 above

The v2 Gate 2 run (your report) FAILED because DISCUS applied an unrequested Debye–Waller factor (B = 1 Å²); see `gate2_prereg/DIAGNOSIS_v2_RESULT.md`. The owner approved **gate2-prereg-v3**: macros use `temp ignore`, and `analyze.py` has a corrected error check plus a `T` criterion. Predictions and thresholds are unchanged.

- `git pull --ff-only`; the head must include the commit titled "Gate 2: register v3". Expected: `pytest tests/test_gate2_prereg.py` → `8 passed`; `gate2_prereg/predictions.json` SHA-256 = `40dd2b140e6380713d275e5d1f59d15ec561da4c20803631bf6b308c31cb488d`.
- **Gate 1 is already done and passed. Do NOT rerun it.** Skip every Gate 1 step. Only Gate 2 is rerun, with a new `RUN_ID`.
- Run order: Step 1 (self-checks), Step 2 (engine identity; the binary hash must still be `d4274f13…45706f3`), then the probe, then `gate2_block_total.mac` **first**. Stop condition (new): its log must contain a line `Temp. factors      : ...` that does NOT end in `used`, and DISCUS must not reject `temp ignore` (any `***APPL***` line). If it does, stop and report the exact log lines; do not try alternatives.
- Then the other three macros and `analyze.py`. A FAIL is still not a stop condition. Report as before.
- Commit as before: the v3 evidence directory, plus the untracked v1 Gate 1 `ENGINE_ENV.txt` directory (`gate1_c1_dyadic/evidence/20261028T232749Z_LAPTOP-NGU3BT0N`). Do not edit the v2 evidence.

## Registration v2 (2026-10-28): if you already have the fresh clone

The first probe (`evidence/20261028T232749Z_LAPTOP-NGU3BT0N`) correctly stopped the run. DISCUS rejected the v1 structure files at `read stru` with error −49: the atom lines were shifted one column against DISCUS's fixed-column format. The reviewer has re-registered as **gate2-prereg-v2**. Only the two `.stru` files and the scripts that write them changed; every predicted number and threshold is identical. The deviation is recorded in `predictions.json → deviation_from_v1`.

In the existing fresh clone `~/gate2_prereg_runs/discus-python`:

```bash
cd ~/gate2_prereg_runs/discus-python
git status --short     # only the untracked v1 probe directory may appear; anything else: STOP
git pull --ff-only     # must fast-forward
git log -1 --format='%H %s'
```

**Keep the v1 probe directory exactly as it is.** It is the evidence for the deviation, and it is committed together with the v2 evidence in Step 6. It holds v1 inputs, so never run `analyze.py` on it. Then continue from Step 1 with **new** `RUN_ID` directories.

## Step 1: Python environment and registration self-checks

Use the venv created in Step 0 (same WSL context as DISCUS):

```bash
cd "$WORK" && source .venv/bin/activate
python -c "import numpy, h5py; print(numpy.__version__, h5py.__version__)"
python -m pytest -q tests/test_gate2_prereg.py      # expect: 8 passed
python benchmarks/discus_6_23_01/gate2_prereg/selftest.py "$(mktemp -d)"   # expect 8 lines starting "ok"
sha256sum benchmarks/discus_6_23_01/gate2_prereg/predictions.json benchmarks/discus_6_23_01/gate1_c1_dyadic/predictions.json
```

Registered values:

| File | SHA-256 |
|---|---|
| `gate2_prereg/predictions.json` | `8c1a66f89ded847bfe2360ce9f757c30e6d329fc975dca968cec115ba6f82760` |
| `gate1_c1_dyadic/predictions.json` | `77d7972effa5eb04e261e530f411cc3113aae025112e3b6bbc3f0c4148480a73` |

If either test step fails or a hash differs: **STOP and report.** Do not "fix" it. A hash difference on a fresh LF clone is not expected; report `git ls-files --eol` for the affected file.

## Step 2: record the engine identity

```bash
RUN_ID=$(date -u +%Y%m%dT%H%M%SZ)_$(hostname)
E2="benchmarks/discus_6_23_01/gate2_prereg/evidence/$RUN_ID"
E1="benchmarks/discus_6_23_01/gate1_c1_dyadic/evidence/$RUN_ID"
test ! -e "$E2" && test ! -e "$E1" && mkdir -p "$E2" "$E1"     # must be new
{
  echo "registered_commit: $(git rev-parse HEAD)"
  echo "discus_path: $(command -v discus_suite)"
  sha256sum "$(command -v discus_suite)"
  ls -l "$(readlink -f "$(command -v discus_suite)")"
  echo "nproc: $(nproc)   OMP_NUM_THREADS=${OMP_NUM_THREADS:-unset}"
  uname -a; cat /etc/os-release | head -2
  python --version; python -c "import numpy, h5py; print('numpy', numpy.__version__, 'h5py', h5py.__version__)"
} > "$E2/ENGINE_ENV.txt" 2>&1
cp "$E2/ENGINE_ENV.txt" "$E1/ENGINE_ENV.txt"
cat "$E2/ENGINE_ENV.txt"
```

Report the **full** SHA-256 of the binary. Earlier work recorded prefix `d4274f13…` and suffix `…45706f3`; say whether it matches. Do not use `/usr/bin/time`: if it is absent, earlier tooling falsely reported exit 127.

## Step 3: Gate 2 runs

DISCUS writes its outputs to the current directory, so always `cd` into the evidence directory first. Copy the inputs with `cp -n` (never overwrite):

```bash
cp -n benchmarks/discus_6_23_01/gate2_prereg/inputs/* "$E2/"
cd "$E2"
run() { local m="$1"; local t0=$(date -u +%s)
        timeout 3600 /usr/local/bin/discus_suite < "$m.mac" > "$m.log" 2>&1; local rc=$?
        echo "$m exit=$rc seconds=$(( $(date -u +%s) - t0 ))" | tee -a RUN_STATUS.txt; }
run gate2_probe
grep -n -i -E "n\[1\]|n\[2\]|n\[3\]|error|\*\*\*" gate2_probe.log   # report these lines verbatim
```

**Probe check.** The three `eval` lines should show n[1] = 4096 atoms, n[2] = 2 atom types, and n[3] = 1 atom per unit cell. If DISCUS reports anything else, or a read error, **STOP and report**: the structure is not being read as 16×16×16 cells, and the experiment would be meaningless.

```bash
run gate2_block_total
run gate2_block_diffuse
run gate2_random_total
run gate2_random_diffuse
ls -l *.h5        # expect four files of ~16.8 MB: gate2_{block,random}_{total,diffuse}.h5
grep -H -E "Atom typ|Used .* of the crystal|Friedel|Fourier technique" gate2_*_{total,diffuse}.log
```

Each macro computes 128³ points × 4096 atoms. Expect seconds to a few minutes each. If one exits non-zero, times out or produces no HDF5: do not rerun into this directory. Keep it, write the reason in `RUN_STATUS.txt`, report, and wait for instructions.

## Step 4: Gate 2 analysis

```bash
cd "<REPO>"
python benchmarks/discus_6_23_01/gate2_prereg/analyze.py "$E2"; echo "analyze exit=$?" | tee -a "$E2/RUN_STATUS.txt"
cat "$E2/GATE2_RESULTS_PREREG.md"
```

Exit codes: 0 = PASS (every criterion, both models), 1 = FAIL, 2 = INVALID. Any of these is a valid outcome to report. **Do not edit thresholds, rerun with changes, or add a "corrected" analysis.**

## Step 5: Gate 1 C1 runs and analysis

```bash
cp -n benchmarks/discus_6_23_01/gate1_c1_dyadic/inputs/* "$E1/"
cd "$E1"
run() { local m="$1"; local t0=$(date -u +%s)
        timeout 3600 /usr/local/bin/discus_suite < "$m.mac" > "$m.log" 2>&1; local rc=$?
        echo "$m exit=$rc seconds=$(( $(date -u +%s) - t0 ))" | tee -a RUN_STATUS.txt; }
run c1_dyadic
run c1_repro
ls -l *.h5        # expect c1_dyadic_{total,diffuse}.h5 and c1_repro_{total,diffuse}.h5
cd "<REPO>"
python benchmarks/discus_6_23_01/gate1_c1_dyadic/c1_dyadic.py analyze "$E1"; echo "analyze exit=$?" | tee -a "$E1/RUN_STATUS.txt"
```

## Step 6: hash manifests and commit the evidence

```bash
( cd "$E2" && sha256sum * > MANIFEST.sha256 )
( cd "$E1" && sha256sum * > MANIFEST.sha256 )
du -sh "$E2" "$E1"; find "$E2" "$E1" -size +90M     # must print nothing (GitHub limit 100 MB/file)
git status --short   # only the new evidence directories (and the v1 probe directory) may appear
git add "$E2" "$E1" benchmarks/discus_6_23_01/gate2_prereg/evidence/20261028T232749Z_LAPTOP-NGU3BT0N
git commit -m "Gate 2 / Gate 1 pre-registered engine runs ($RUN_ID): evidence only

Gate 2 analyze.py status: <PASS|FAIL|INVALID>; Gate 1 c1_dyadic status: <...>.
Engine: discus_suite 6.23.01 sha256 <full hash>. No registered file modified."
git push origin science/discus-3d-dpdf-benchmarks     # plain push; never --force
```

The evidence directories (`evidence/`, not `runs/`) are deliberately outside the repo-wide `runs/` gitignore. If `git status` shows anything other than the two evidence directories, **do not commit it**; report it.

## Step 7: report back to the owner

Paste verbatim, without summarizing away numbers:
1. `REGISTERED_COMMIT`, the commit you pushed, and `ENGINE_ENV.txt`;
2. the probe-log `eval` lines;
3. `RUN_STATUS.txt` from both directories;
4. the full `GATE2_RESULTS_PREREG.md` and `C1_RESULTS_PREREG.md`;
5. any warnings or errors in the logs, and anything that differed from this document.

State the outcome as the scripts printed it. Do **not** write "Gate 2 passed" or "failed" into any project document, and do not update `docs/SCIENTIFIC_VALIDATION_REPORT.md`, `GATE2_RESULTS.md` or `REVIEW_2026-10-28.md`. The owner and the reviewer will do that after reading the evidence.

## Stop conditions (report and wait)

- dirty working tree, non-fast-forward pull, failing self-checks, or hash mismatch (Steps 0–1);
- binary not 6.23.01, or not at `/usr/local/bin/discus_suite`;
- probe reports something other than 4096 / 2 / 1, or any read error;
- engine exit ≠ 0, timeout, missing or unexpectedly sized HDF5;
- `analyze.py` exit 2 (INVALID);
- anything requiring deletion, installation, or editing a registered file.

A FAIL (exit 1) is **not** a stop condition for committing evidence: commit and report it like a PASS.
