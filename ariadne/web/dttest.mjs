import fs from "fs";
import { build, calibrate, thresholdFor, boundaryX0, setDT, X0_UNIT, K_DEFAULT } from "./engine.mjs";
const raw = fs.readFileSync("/tmp/W_stanford_raw.csv","utf8").trim().split("\n").map(l=>l.split(",").map(Number));
const { sp, n } = build(raw); const ONSET = 38;
for (const dt of [0.1, 0.2, 0.25, 0.4]) {
  setDT(dt);
  let t=Date.now();
  const cal = calibrate(sp, n, ONSET, K_DEFAULT, null);
  const B = cal.ok ? boundaryX0(sp, n, ONSET, cal.x0Healthy, K_DEFAULT, 7, null) : null;
  const ltg = cal.ok ? thresholdFor(sp, n, ONSET, cal.x0Healthy, K_DEFAULT, 0.941, 0.3, 6, null) : null;
  const cnb = B===null?null:(B/(2.541*X0_UNIT));
  console.log(`dt=${dt}  x0=${cal.ok?cal.x0Healthy:"FAIL"}  boundary=${B?B.toFixed(5):"-"}  CNB=${cnb?cnb.toFixed(3):"-"} (py 0.344)  LTG=${ltg} (py 0.719)   ${((Date.now()-t)/1000).toFixed(0)}s`);
}
