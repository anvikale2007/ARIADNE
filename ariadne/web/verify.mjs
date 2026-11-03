import fs from "fs";
import { build, calibrate, thresholdFor, boundaryX0, X0_UNIT, K_DEFAULT } from "./engine.mjs";

const raw = fs.readFileSync("/tmp/W_stanford_raw.csv", "utf8").trim()
  .split("\n").map(l => l.split(",").map(Number));
const { sp, n } = build(raw);
const ONSET = 38;
console.log(`n=${n}  onset=${ONSET}`);

let sims = 0; const tick = () => { sims++; };
let t = Date.now();
const cal = calibrate(sp, n, ONSET, K_DEFAULT, tick);
console.log(`calibrate -> x0_healthy ${cal.x0Healthy}   (python: -2.12)   healthy recruited ${cal.healthyRecruited}   ${sims} sims  ${((Date.now()-t)/1000).toFixed(1)}s`);

t = Date.now(); sims = 0;
const B = boundaryX0(sp, n, ONSET, cal.x0Healthy, K_DEFAULT, 7, tick);
console.log(`boundary x0 reduction = ${B.toFixed(5)}   ${sims} sims  ${((Date.now()-t)/1000).toFixed(1)}s`);

// pure-excitability drugs via the ray identity; python values for comparison
const PURE = { cenobamate:[2.541,0.344], carbamazepine:[0.534,1.594], lacosamide:[0.530,1.594], phenytoin:[0.369,2.281] };
for (const [k,[xe,py]] of Object.entries(PURE)) {
  const th = B / (xe * X0_UNIT);
  console.log(`  ${k.padEnd(14)} ray ${th.toFixed(3).padStart(6)}   python ${String(py).padStart(6)}   diff ${(th-py).toFixed(3)}`);
}

// mixed / coupling-only drugs need their own bisection
const MIX = { lamotrigine:[0.941,0.3,0.719], valproate:[0.450,0.2,1.531], pregabalin:[0,0.6,2.156], levetiracetam:[0,0.403,3.219] };
for (const [k,[xe,ce,py]] of Object.entries(MIX)) {
  t = Date.now(); sims = 0;
  const th = thresholdFor(sp, n, ONSET, cal.x0Healthy, K_DEFAULT, xe, ce, 6, tick);
  console.log(`  ${k.padEnd(14)} sim ${String(th).padStart(6)}   python ${String(py).padStart(6)}   diff ${th===null?"n/a":(th-py).toFixed(3)}   ${sims} sims ${((Date.now()-t)/1000).toFixed(1)}s`);
}
