"""Inline the sources into one page: webui/index.html (served by GitHub Pages).
Also writes dist/artifact.html (body only, for a claude.ai artifact preview; not committed)."""
from pathlib import Path
import html
root = Path(__file__).parent
src = root / "src"
engine = (root / "engine.py").read_text(encoding="utf-8")
assert "</script" not in engine.lower()
body = (src / "shell.html").read_text(encoding="utf-8")
parts = {"/*STYLE*/": src / "style.css", "/*CORE*/": src / "core.js", "/*EXAMPLES*/": src / "examples.js",
         "/*VIEWER*/": src / "viewer.js", "/*APP*/": src / "app.js"}
for k, p in parts.items():
    t = p.read_text(encoding="utf-8")
    assert "</script" not in t.lower(), p
    body = body.replace(k, t)
body = body.replace("/*ENGINE*/", engine)
out = root / "dist"
out.mkdir(exist_ok=True)
(out / "artifact.html").write_text(body, encoding="utf-8")
head_end = body.index("<header")
full = ("<!doctype html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n"
        "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1, viewport-fit=cover\">\n"
        "<style>:root{padding-top:env(safe-area-inset-top,0px);padding-bottom:env(safe-area-inset-bottom,0px)}"
        "body{margin:0}[hidden]{display:none!important}img{max-width:100%}</style>\n"
        + body[:head_end] + "</head>\n<body>\n" + body[head_end:] + "\n</body>\n</html>\n")
(root / "index.html").write_text(full, encoding="utf-8")
print("built", len(full) // 1024, "KB")
