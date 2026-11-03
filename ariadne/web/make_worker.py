"""engine.mjs -> engine.worker.js (classic script for a Blob Worker).

The worker cannot be a module in a Blob URL across every browser, so the module's `export ` keywords
are stripped and the file is used as a classic script. Nothing else is changed, so the worker runs
byte-identical physics to the module the Node tests exercise.
"""
import os, re
HERE = os.path.dirname(os.path.abspath(__file__))
s = open(f"{HERE}/engine.mjs", encoding="utf-8").read()
s = re.sub(r"^export\s+", "", s, flags=re.M)
assert "export " not in s, "an export survived the strip"
open(f"{HERE}/engine.worker.js", "w", encoding="utf-8").write(s)
print(f"wrote engine.worker.js  {len(s)/1024:.1f} KiB")
