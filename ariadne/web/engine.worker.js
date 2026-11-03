/* ARIADNE in the browser.
 *
 * A direct port of ariadne/model.py (network Epileptor-2D), ariadne/metrics.py and the calibration
 * and bisection in ariadne/engine.py, so an uploaded connectome runs the SAME equations as the
 * Python pipeline rather than a lookalike.
 *
 * Two deliberate differences from Python, both stated in the UI:
 *   1. Noise is a seeded JS PRNG, not numpy's Mersenne Twister, so the noise realisation differs.
 *      Amplitude is 1e-4, far below the ~11% bisection resolution, but a threshold can land one
 *      bisection step away from the Python value.
 *   2. Thresholds for drugs with no coupling effect are taken from the containment boundary by the
 *      identity threshold = boundary / (x0_effect * X0_UNIT), which the project verified holds
 *      exactly. That turns 6 bisections into 1 and is what makes this fast enough to run live.
 */

// ---- model constants, from ariadne/model.py ----
const X0_EPILEPTOGENIC = -1.6;
const R_SLOW = 0.00035;
const IEXT = 3.1;
// 0.2 ms verified to give results identical to 0.1 ms on the Stanford case (calibration, boundary and
// every drug threshold match exactly) at twice the speed. 0.25 ms already shifts lamotrigine one step.
let DT = 0.2;                // ms
const NOISE = 1e-4;
const RECORD_EVERY_MS = 10.0;

// ---- engine constants, from ariadne/engine.py ----
const X0_UNIT = 0.004;
const COUPLING_UNIT = 0.06;
const SPREAD_FRACTION = 0.5;
// 45 s verified to give results identical to 60 s on the Stanford case - calibration, boundary and all
// seven drug thresholds match exactly - at 0.82x the time. 35 s already shifts pregabalin and
// levetiracetam one bisection step each, so 45 s is the edge, not a guess.
const LENGTH_MS = 45000;
const S_MAX = 4.0;
const K_DEFAULT = -5.0;

// ---- metrics constants, from ariadne/metrics.py ----
const BURN_IN_MS = 15000;
const RECRUIT_AMPLITUDE = 1.0;

function mulberry32(a) {
  return function () {
    a |= 0; a = (a + 0x6D2B79F5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}
/** Box-Muller, cached, so each call costs about one uniform draw. */
function gaussian(rng) {
  let spare = null;
  return function () {
    if (spare !== null) { const v = spare; spare = null; return v; }
    let u, v, s;
    do { u = rng() * 2 - 1; v = rng() * 2 - 1; s = u * u + v * v; } while (s >= 1 || s === 0);
    const m = Math.sqrt(-2 * Math.log(s) / s);
    spare = v * m; return u * m;
  };
}

/** connectome.prepare(): log1p, symmetrise, zero diagonal, scale so the largest row sum is 1. */
function prepare(rows) {
  const n = rows.length;
  const W = new Float64Array(n * n);
  for (let i = 0; i < n; i++)
    for (let j = 0; j < n; j++)
      W[i * n + j] = Math.log1p(Math.max(0, rows[i][j]));
  for (let i = 0; i < n; i++)
    for (let j = i + 1; j < n; j++) {
      const m = (W[i * n + j] + W[j * n + i]) / 2;
      W[i * n + j] = m; W[j * n + i] = m;
    }
  for (let i = 0; i < n; i++) W[i * n + i] = 0;
  let mx = 0;
  for (let i = 0; i < n; i++) {
    let s = 0;
    for (let j = 0; j < n; j++) s += W[i * n + j];
    if (s > mx) mx = s;
  }
  if (mx > 0) for (let k = 0; k < W.length; k++) W[k] /= mx;
  return { W, n };
}

/** Compressed sparse rows — the coupling term is the only O(n^2) part of the inner loop. */
function csr(W, n) {
  const idx = [], val = [], ptr = new Int32Array(n + 1);
  const rsum = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      const w = W[i * n + j];
      if (w !== 0) { idx.push(j); val.push(w); rsum[i] += w; }
    }
    ptr[i + 1] = idx.length;
  }
  return { idx: Int32Array.from(idx), val: Float64Array.from(val), ptr, rsum };
}

/** Prepare raw rows straight into the sparse form the simulator wants. */
function setDT(v){ DT = v; }

function build(rows) {
  const { W, n } = prepare(rows);
  return { sp: csr(W, n), n, W };
}

/**
 * One simulation. Returns the recruited mask and count.
 * Only min/max per region after burn-in are tracked, because that is all metrics.analyse needs.
 */
function simulate(sp, n, onset, x0Shift, couplingScale, x0Healthy, K, seed, stopAt) {
  const { idx, val, ptr, rsum } = sp;
  const x0 = new Float64Array(n).fill(x0Healthy + x0Shift);
  if (onset !== null && onset !== undefined) x0[onset] = X0_EPILEPTOGENIC + x0Shift;
  const Keff = K * couplingScale;

  const x1 = new Float64Array(n).fill(-1.8);
  const z = new Float64Array(n).fill(3.0);
  const a1 = new Float64Array(n), a2 = new Float64Array(n);
  const b1 = new Float64Array(n), b2 = new Float64Array(n);
  const xp = new Float64Array(n), zp = new Float64Array(n);
  const lo = new Float64Array(n).fill(Infinity), hi = new Float64Array(n).fill(-Infinity);
  const got = new Uint8Array(n);      // already past the recruitment amplitude
  let count = 0;
  // Amplitude over the measured window can only grow as the window extends, so once more than half
  // the regions are recruited the verdict is fixed and the rest of the run cannot change it.
  const limit = (stopAt === undefined || stopAt === null) ? Infinity : stopAt;

  const randn = gaussian(mulberry32(seed * 7919 + 13));
  const nsteps = Math.round(LENGTH_MS / DT);
  const recEvery = Math.round(RECORD_EVERY_MS / DT);
  const burnSteps = Math.round(BURN_IN_MS / DT);
  const sq = Math.sqrt(DT);

  const dfun = (X, Z, dX, dZ) => {
    for (let i = 0; i < n; i++) {
      let c = -rsum[i] * X[i];
      for (let p = ptr[i]; p < ptr[i + 1]; p++) c += val[p] * X[idx[p]];
      const xi = X[i], zi = Z[i];
      const f = xi < 0 ? xi * xi + 2.0 * xi : -0.6 * (zi - 4.0) * (zi - 4.0) + 5.0 * xi;
      dX[i] = 1.0 - zi + IEXT - f * xi;
      let e = 0.0;
      if (zi < 0) { const z2 = zi * zi, z3 = z2 * zi; e = -0.1 * z3 * z3 * zi; }
      dZ[i] = R_SLOW * (4.0 * (xi - x0[i]) + e - zi + Keff * c);
    }
  };

  for (let k = 0; k < nsteps; k++) {
    dfun(x1, z, a1, a2);
    for (let i = 0; i < n; i++) { xp[i] = x1[i] + DT * a1[i]; zp[i] = z[i] + DT * a2[i]; }
    dfun(xp, zp, b1, b2);
    for (let i = 0; i < n; i++) {
      x1[i] += 0.5 * DT * (a1[i] + b1[i]) + NOISE * sq * randn();
      z[i] += 0.5 * DT * (a2[i] + b2[i]);
    }
    if (k >= burnSteps && k % recEvery === 0) {
      for (let i = 0; i < n; i++) {
        const v = x1[i];
        if (v < lo[i]) lo[i] = v;
        if (v > hi[i]) hi[i] = v;
        if (!got[i] && hi[i] - lo[i] > RECRUIT_AMPLITUDE) { got[i] = 1; count++; }
      }
      if (count > limit) break;        // exact: the verdict can no longer change
    }
  }
  return { mask: got, count };
}

const spreads = (count, n) => count > SPREAD_FRACTION * n;
/** Run only as far as the spread/no-spread verdict needs. */
const verdict = (sp, n, onset, shift, scale, x0h, K, seed) =>
  spreads(simulate(sp, n, onset, shift, scale, x0h, K, seed, SPREAD_FRACTION * n).count, n);

/**
 * ONE spread/no-spread verdict, exported so a worker pool can evaluate many points of a search at
 * once instead of walking a bisection one step at a time.
 *
 * Calibration, the containment boundary and a drug's threshold are all the same shape of problem:
 * find the smallest value of one parameter at which the seizure stops spreading, where spreading is
 * monotone in that parameter. A bisection is the fewest TOTAL simulations, but they are strictly
 * sequential, so on an 8-core machine seven of eight cores sit idle through the whole search. Fanning
 * a grid of points out across the pool costs more simulations and far less wall-clock time, which is
 * the trade that matters when someone is standing at a booth waiting for the answer.
 */
function probeSpread(sp, n, onset, x0Shift, couplingScale, x0Healthy, K) {
  return verdict(sp, n, onset, x0Shift, couplingScale, x0Healthy, K, 0);
}

/**
 * engine.calibrate(): the lowest background excitability at which one seeded region recruits most of
 * the network. Python scans the grid from the bottom and takes the first hit; this bisects the same
 * grid, which is identical because spreading is monotonic in x0 and costs ~4 sims instead of ~15.
 */
function calibrate(sp, n, onset, K, onStep) {
  const grid = [];
  for (let v = -2.20; v <= -2.0585; v += 0.01) grid.push(Math.round(v * 1000) / 1000);
  let lo = 0, hi = grid.length - 1, found = null;
  if (verdict(sp, n, onset, 0, 1, grid[0], K, 0)) found = grid[0];
  else {
    onStep && onStep();
    while (lo < hi) {
      const mid = (lo + hi) >> 1;
      const sp2 = verdict(sp, n, onset, 0, 1, grid[mid], K, 0);
      onStep && onStep();
      if (sp2) hi = mid; else lo = mid + 1;
    }
    if (verdict(sp, n, onset, 0, 1, grid[lo], K, 0)) found = grid[lo];
    onStep && onStep();
  }
  if (found === null) return { ok: false };
  const healthy = simulate(sp, n, null, 0, 1, found, K, 0, SPREAD_FRACTION * n);
  onStep && onStep();
  return { ok: true, x0Healthy: found, K, healthyRecruited: healthy.count };
}

/** Bisection on drug strength, matching engine._drug_job. */
function thresholdFor(sp, n, onset, x0Healthy, K, x0Effect, couplingEffect, iters, onStep) {
  const at = (s) => {
    const shift = -s * x0Effect * X0_UNIT;
    const scale = Math.max(0, 1 - s * couplingEffect * COUPLING_UNIT);
    const sp2 = verdict(sp, n, onset, shift, scale, x0Healthy, K, 0);
    onStep && onStep();
    return sp2;
  };
  if (at(S_MAX)) return null;
  let lo = 0, hi = S_MAX;
  for (let i = 0; i < iters; i++) {
    const mid = (lo + hi) / 2;
    if (at(mid)) lo = mid; else hi = mid;
  }
  return Math.round(((lo + hi) / 2) * 1000) / 1000;
}

/** The containment boundary at full coupling: the smallest x0 reduction that stops the spread. */
function boundaryX0(sp, n, onset, x0Healthy, K, iters, onStep) {
  const at = (shift) => {
    const sp2 = verdict(sp, n, onset, -shift, 1, x0Healthy, K, 0);
    onStep && onStep();
    return sp2;
  };
  const top = S_MAX * 2.541 * X0_UNIT;          // the strongest drug in the panel, at S_MAX
  if (at(top)) return null;
  let lo = 0, hi = top;
  for (let i = 0; i < iters; i++) {
    const mid = (lo + hi) / 2;
    if (at(mid)) lo = mid; else hi = mid;
  }
  return (lo + hi) / 2;
}
