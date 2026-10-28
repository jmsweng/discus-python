# Instructions for the local agent: commit and push Claude's DISCUS 3D work

Written 2026-10-22. Claude could not push (no access to the private repo from its workspace), so the work
sits **untracked** in this working tree. Your job: put it on two branches, push, open PRs, then STOP.

## Rules
- Do **not** merge into `main`, force-push, rewrite history, rebase `main`, or reopen/re-merge PRs #1–#3.
- Stage **explicit paths only**. Never `git add -A` / `git add .` / `git commit -a`: 15 tracked files show as
  "modified" only because of CRLF-vs-LF; they are not Claude's changes. Verify: `git diff --ignore-cr-at-eol --stat` prints nothing.
- Do not commit `spike_discus_3d/data_run*/`, `*.tgz`, `*.pdf`, `suite.hlp`, `DiscusSuite.txt`, `share_listing.txt`,
  `PR2_REVIEW.md`, `PR3_REVIEW.md`, or this file. (Most are in `.git/info/exclude`; still check `git status` before every commit.)
- No `sudo`, no installs, no credentials in commits. Use this machine's existing SSH key; do not copy it anywhere.

## 0. Check the base
```bash
cd "/mnt/c/Users/jmswe/Documents/Programming projects/discus-python"
git fetch origin
git log --oneline -1 origin/main        # expected ef06490 (PR #4 merge) or newer; record the SHA you use
git status -sb
```
If `origin/main` is newer than `ef06490`, use it, and re-run the tests (step 1) on that base.

## 1. Test before committing (WSL venv)
```bash
source ~/projects/discuspdf_starter/.venv/bin/activate
python -m pytest -q          # Claude saw: 69 passed, 1 skipped (Linux VM, Python 3.10). Report any difference.
```
Optional (real engine): `DISCUS_INTEGRATION=1 python -m pytest -q tests/test_integration_discus.py`.
It is **expected to be unreliable**: the committed 1D reference was made with qmax=30, the test uses qmax=25.0. Report the result; do not "fix" it here.

## 2. PR A — `spike/verify-real-discus-3d` (evidence + probe kit, no library code)
```bash
git switch -c spike/verify-real-discus-3d origin/main
git add spike_discus_3d/AUDIT_BASELINE.md spike_discus_3d/FINDINGS.md spike_discus_3d/PROBE_PROTOCOL.md \
        spike_discus_3d/probe_environment.sh spike_discus_3d/run_cases.sh spike_discus_3d/seed_test.sh \
        spike_discus_3d/inspect_output.py spike_discus_3d/macros
git status -s                      # only the paths above may be staged
git commit -F - <<'MSG'
Add verified DISCUS 6.23.01 3D probe kit and findings

Real runs of DISCUS_SUITE 6.23.01 (WSL2): 3D hkl-grid single-crystal Fourier,
nexus output layout, Bragg subtraction via 'set aver', independent kinematic
cross-check, skewed/anisotropic grid tests. See FINDINGS.md.

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_015dwUbYBUwUi8AvGeYXKGNX
MSG
git push -u origin spike/verify-real-discus-3d
```

## 3. PR B — `feat/discus-3d-native-adapter` (library code + real fixtures + tests)
Branch from the **same** base, not from PR A:
```bash
git switch -c feat/discus-3d-native-adapter origin/main
git add src/discuspdf/discus_output.py tests/test_discus_output.py tests/fixtures/discus_6.23.01
git status -s
python -m pytest -q                # must still pass on this branch
git commit -F - <<'MSG'
Add DISCUS 6.23.01 common/nexus reader and fractional-coordinate 3D-ΔPDF

Reads real DISCUS output with source-metadata validation, keeps the cell metric
(no relabelling hkl as Q), transforms in fractional coordinates, refuses skewed,
uncentered or Bragg-present input by default. Real fixtures: two 0.5 MB DISCUS files.

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_015dwUbYBUwUi8AvGeYXKGNX
MSG
git push -u origin feat/discus-3d-native-adapter
```
This adds new files only. It deliberately does **not** change `diffuse3d.py`, `volume.py`, `reciprocal.py`, `__init__.py`, or docs.

## 4. Open the PRs (base `main`)
Use `gh pr create --base main --head <branch> --title ... --body-file ...` if `gh` is installed, otherwise open
`https://github.com/jmsweng/discus-python/pull/new/<branch>`. Body must contain these sections, copied from `spike_discus_3d/FINDINGS.md`:
**Observed / Independently verified / Assumed / Not yet tested**, plus the exact base SHA, open limitations (unreproducible `seed`, unexplained ~1e-3 residual vs
independent calculation, DISCUS 3DPDF header step vs FFT spacing, qmax=30 vs 25 mismatch), and the CI run link once it appears.
PR B should say: "Depends on no other PR; PR A holds the evidence." End each PR body with:
```
🤖 Generated with [Claude Code](https://claude.com/claude-code)

https://claude.ai/code/session_015dwUbYBUwUi8AvGeYXKGNX
```

## 5. Then STOP and report to the owner
Return: both PR URLs, the base SHA, branch tip SHAs, `pytest` summary lines, CI status, and anything that differed from the above.
**Do not merge.** The owner approves merges.

## Optional, only if the owner asks
- `bash spike_discus_3d/seed_test.sh` (finds out why `seed 1` is not reproducible) and send back `spike_discus_3d/seed_test.tgz`.
