"""Profile the selftest components."""
import sys
import time
import tempfile
import json
import hashlib
from pathlib import Path

sys.path.insert(0, str(Path(r"C:\Users\jmswe\.local\share\opencode\worktree\9d7e5c\calm-sailor\benchmarks\discus_6_23_01\gate0_prereg")))
import gate0_selftest as ST
import gate0_analyze as AN
import oracle as O
import numpy as np

work = Path(tempfile.mkdtemp(prefix="gate0_profile_"))

# Time _make_clean_run (includes oracle computation)
t0 = time.time()
ST._make_clean_run(work)
t1 = time.time()
print(f"_make_clean_run (128 oracle calls + h5 write): {t1 - t0:.1f}s")

# Time _make_gaussian_run (reuses oracle cache)
t0 = time.time()
ST._make_gaussian_run(work)
t1 = time.time()
print(f"_make_gaussian_run (reuses cache + h5 write): {t1 - t0:.1f}s")

# Time a single analysis
pred = json.loads((ST.HERE / "predictions.json").read_text())
for name in pred["code_sha256"]:
    p = ST.HERE / name
    if p.exists():
        pred["code_sha256"][name] = hashlib.sha256(p.read_bytes()).hexdigest()

from pathlib import Path as _Path
import json as _json

orig_read_text = _Path.read_text

def patched_read_text(self, *args, **kwargs):
    if self.name == "predictions.json" and str(self.parent) == str(ST.HERE):
        return _json.dumps(pred, indent=1) + "\n"
    return orig_read_text(self, *args, **kwargs)

_Path.read_text = patched_read_text

t0 = time.time()
res = AN.analyze(work / "clean")
t1 = time.time()
print(f"AN.analyze (clean, 128 h5 reads): {t1 - t0:.1f}s")
print(f"Status: {res['status']}")
for arm, arm_data in res.get("arms", {}).items():
    for c in arm_data["criteria"]:
        status = "PASS" if c["pass"] else "FAIL"
        print(f"  {c['criterion']}: {c['value']:.4e} [{status}]")

_Path.read_text = orig_read_text
