# Integration history — PRs #1–#3

This note records how the October 22, 2026 PDF/3D infrastructure changes
were integrated. **The Git graph is intact and was not rewritten.** GitHub's
PR status is separate from Git's commit ancestry.

| Original PR | Feature branch tip incorporated | Integration commit on `main` | GitHub PR status |
|---|---|---|---|
| [#1](https://github.com/jmsweng/discus-python/pull/1) — 1D PDF runner | `407a7443a929467a732f4799e9610cd57efae793` | [`f8b2ba1`](https://github.com/jmsweng/discus-python/commit/f8b2ba1feef23ac50bcf4cc581514fa6de376fe3) | Merged |
| [#2](https://github.com/jmsweng/discus-python/pull/2) — 3D reciprocal infrastructure | `ca16d93f09cc477e532bafc5fcba24cbd126270e` | [`b6bc243`](https://github.com/jmsweng/discus-python/commit/b6bc243e633b6ba70046542294273bec081fa386) | Closed (GitHub does not mark it merged) |
| [#3](https://github.com/jmsweng/discus-python/pull/3) — FFT diagnostics | `d97626f52a7596efd4003843d832bd7cf0692f6f` | [`31cfa33`](https://github.com/jmsweng/discus-python/commit/31cfa335c3baa3cf9386ac9ba9ba08c31415aaf4) | Closed (GitHub does not mark it merged) |

## Verified ancestry

The integration commits above **all have two parents**. In order:

- `f8b2ba1`: parents `2cb450839751` and `407a7443a929`.
- `b6bc243`: parents `f8b2ba1feef2` and `ca16d93f09cc`.
- `31cfa33`: parents `b6bc243e633b` and `d97626f52a75`.

Thus, PR #2 and PR #3 were **integrated by ordinary Git merge commits**,
rather than by using GitHub's pull-request merge operation. Their code
and original feature-branch ancestry are present in `main`; GitHub's
`merged=false` flags refer to the PR records, not the contents or ancestry
of the repository.

## Notes for future agents

- Do **not** reopen or re-merge #2 or #3 just to change GitHub's status.
- Do **not** force-push or rebase `main` to cosmetically alter these old commits.
- Keep the original PRs closed and link the integration commit when referring
  to the implementation.
- Subsequent cleanup should be proposed in a new PR against `main`.
- Scientific status is unchanged: 1D PDF integration has local DISCUS evidence;
  3D-ΔPDF does not yet have a verified end-to-end DISCUS benchmark.
