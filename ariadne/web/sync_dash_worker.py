"""Rewrite the dashboard's embedded WORKER_SRC from engine.mjs + the worker protocol tail.

dash.html carries the whole simulation engine inline, because the page must work offline from a
laptop with no server. That copy drifted from engine.mjs once already (an engine change was made and
the page kept running the old physics), so it is now GENERATED rather than maintained: engine.mjs is
the only source of the model, and this script restamps the page from it.
"""
import os, re

HERE = os.path.dirname(os.path.abspath(__file__))
DASH = "/tmp/claude-0/-home-claude/4dd36763-85b0-510b-8425-34461c4e8d03/scratchpad/story/dash.html"

engine = open(f"{HERE}/engine.worker.js", encoding="utf-8").read()
tail = re.search(r'WORKER_TAIL = r"""(.*?)"""', open(f"{HERE}/inject_upload.py", encoding="utf-8").read(),
                 re.S).group(1)

# Backticks and ${ would end or interpolate the template literal the page wraps this in.
body = engine + tail
assert "`" not in body and "${" not in body, "worker source would break the template literal"

d = open(DASH, encoding="utf-8").read()
m = re.search(r"(const WORKER_SRC = `)(.*?)(\n`;)", d, re.S)
assert m, "WORKER_SRC block not found in dash.html"
old = m.group(2)
d = d[:m.start(2)] + "\n" + body + d[m.end(2):]
open(DASH, "w", encoding="utf-8").write(d)
print(f"WORKER_SRC: {len(old)/1024:.1f} KiB -> {len(body)/1024:.1f} KiB")
for k in ("LENGTH_MS = ", "DT = "):
    for line in body.split("\n"):
        if line.strip().startswith(("const " + k, "let " + k)):
            print("  ", line.strip())
print("  probe command:", "cmd === \"probe\"" in body)
