# DISCUS + discuspdf setup plan

- [x] 1. Install apt prerequisites in Ubuntu (git, unzip, python3-venv, python3-pip, curl) — done with sudo password
- [x] 2. Run official `bbb_install_script.sh` in Ubuntu to install DISCUS suite — done (compiled from source, v6.23.01, `/usr/local/bin/discus_suite`)
- [x] 3. Locate the DISCUS executable and verify it runs (banner + clean exit) — done (banner v6.23.01, exit 0)
- [x] 4. Copy starter into `~/projects/discuspdf_starter`, create `.venv`, `pip install -e '.[notebook,dev]'` — done (Python 3.10.12, pip 26.2.2, all extras)
- [x] 5. Run `python -m pytest -q` — **7 passed**
- [x] 6. Generate + inspect the 1D PDF macro; validate release-compatible CLI flags — done (CLI `-noautorun -macro` confirmed in manual; macro syntax corrected to v6 powder→PDF workflow)
- [x] 7. Manually run the generated macro with the installed suite — done (real `result.pdf`, 601 lines)
- [x] 8. Run a real 1D PDF job through `PDFRunner` from Python; verify finite numeric `result.pdf` — done (601 samples, finite, increasing r)
- [x] 9. Verify notebook 1 (real, 2 plots) and notebook 2 (synthetic, 1 plot) with the `.venv` kernel — done
- [x] 10. Write `SETUP_REPORT.md` and `docs/VERIFIED_DISCUS_MACRO.md` — done (both written + synced to WSL)
