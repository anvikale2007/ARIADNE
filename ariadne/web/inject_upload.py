"""Add the 'Load a connectome' feature to the ARIADNE dashboard.

The physics is a verified port of the Python pipeline (web/engine.mjs, checked against
results/stanford_hardi_onsetP39.json: calibration matches exactly, thresholds within one
bisection step). It runs in Web Workers so the page stays responsive, and the whole thing stays
in one self-contained HTML file so the dashboard still works offline from a laptop.
"""
import os, re

HERE = os.path.dirname(os.path.abspath(__file__))
DASH = "/tmp/claude-0/-home-claude/4dd36763-85b0-510b-8425-34461c4e8d03/scratchpad/story/dash.html"

worker_src = open(f"{HERE}/engine.worker.js", encoding="utf-8").read()

WORKER_TAIL = r"""
/* ---- worker protocol: one job per message ---- */
let SP = null, N = 0;
self.onmessage = (e) => {
  const m = e.data;
  if (m.cmd === "init") {
    const b = build(m.rows); SP = b.sp; N = b.n;
    postMessage({ ok: true, n: N, edges: countEdges(b.W, b.n), W: Array.from(b.W) });
    return;
  }
  if (m.cmd === "calibrate") {
    postMessage({ job: m.job, result: calibrate(SP, N, m.onset, m.K, () => postMessage({ tick: 1 })) });
    return;
  }
  if (m.cmd === "boundary") {
    postMessage({ job: m.job, result: boundaryX0(SP, N, m.onset, m.x0h, m.K, m.iters, () => postMessage({ tick: 1 })) });
    return;
  }
  if (m.cmd === "threshold") {
    postMessage({ job: m.job, drug: m.drug,
      result: thresholdFor(SP, N, m.onset, m.x0h, m.K, m.xe, m.ce, m.iters, () => postMessage({ tick: 1 })) });
    return;
  }
  /* One point of a search. The orchestrator fans a whole grid of these across the pool, so a search
     that was a sequential bisection becomes two or three parallel rounds. */
  if (m.cmd === "probe") {
    postMessage({ job: m.job, i: m.i,
      spreads: probeSpread(SP, N, m.onset, m.shift, m.scale, m.x0h, m.K) });
    postMessage({ tick: 1 });
    return;
  }
  if (m.cmd === "mask") {
    const shift = -m.s * m.xe * X0_UNIT;
    const scale = Math.max(0, 1 - m.s * m.ce * COUPLING_UNIT);
    const r = simulate(SP, N, m.onset, shift, scale, m.x0h, m.K, 0);
    postMessage({ job: m.job, drug: m.drug, mask: Array.from(r.mask), count: r.count });
    postMessage({ tick: 1 });
  }
};
function countEdges(W, n) { let c = 0; for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) if (W[i * n + j] > 0) c++; return c; }
"""

CSS = """
.up{display:flex;flex-direction:column;gap:9px}
.up input[type=file]{display:none}
.upbtn{background:var(--panel);border:1px dashed var(--line2);border-radius:6px;padding:13px;
  color:var(--dim);font-size:12.5px;text-align:center;cursor:pointer;transition:border-color .15s,color .15s}
.upbtn:hover{border-color:var(--ok);color:var(--ok)}
.bar{height:4px;background:var(--panel2);border-radius:2px;overflow:hidden}
.bar i{display:block;height:100%;width:0;background:var(--ok);transition:width .25s}
.go{background:var(--ok);color:#06121a;border:0;border-radius:5px;padding:9px;font-weight:600;
  font-size:12.5px;cursor:pointer}
.go[disabled]{background:var(--line2);color:var(--faint);cursor:default}
.warnbox{background:rgba(240,180,71,.1);border-radius:5px;padding:8px 10px;font-size:11.5px;color:var(--warn);line-height:1.45}
"""

MARKUP = """
    <div class="blk">
      <h2>Load a connectome <span id="upstat"></span></h2>
      <div class="up">
        <label class="upbtn" for="upfile" id="uplabel">Choose a CSV file</label>
        <input type="file" id="upfile" accept=".csv,.txt,text/csv">
        <p style="font-size:11.5px;color:var(--faint);line-height:1.5">An n&times;n matrix of streamline
          counts from any tractography pipeline. Runs the real model in your browser &mdash; about a minute.</p>
        <div id="upform" hidden style="display:flex;flex-direction:column;gap:9px">
          <label class="lbl" for="uponset">Seizure onset region</label>
          <select id="uponset"></select>
          <button class="go" id="uprun">Run ARIADNE</button>
        </div>
        <div id="upprog" hidden style="display:flex;flex-direction:column;gap:6px">
          <div class="bar"><i id="upbar"></i></div>
          <p style="font-family:var(--mono);font-size:11px;color:var(--dim)" id="upmsg"></p>
        </div>
        <div id="upnote" hidden class="warnbox"></div>
      </div>
    </div>"""

JS = r"""
/* ---------------- load a connectome and run the real model in the browser ---------------- */
const WORKER_SRC = __WORKER_SRC__;
/* Use the whole machine. Every job is pure CPU in a worker, the page itself stays idle while they
   run, and the searches below are built to keep all of them busy. */
const NWORK = Math.max(2, Math.min(8, navigator.hardwareConcurrency || 4));
let UPROWS = null, UPN = 0;

function parseMatrix(text){
  const rows = text.trim().split(/\r?\n/).filter(l=>l.trim())
    .map(l=>l.split(/[,;\t]+/).map(Number));
  const n = rows.length;
  if (n < 8) throw new Error("need at least 8 regions, found " + n);
  if (n > 200) throw new Error("at most 200 regions, found " + n);
  for (const r of rows){
    if (r.length !== n) throw new Error(`not square: ${n} rows but a row has ${r.length} columns`);
    for (const v of r) if (!isFinite(v) || v < 0) throw new Error("values must be non-negative numbers");
  }
  return rows;
}
function fibSphere(n){
  const out=[], ga=Math.PI*(3-Math.sqrt(5));
  for(let i=0;i<n;i++){
    const y=1-(i/(n-1))*2, r=Math.sqrt(Math.max(0,1-y*y)), th=ga*i;
    out.push([Math.cos(th)*r, y, Math.sin(th)*r]);
  }
  return out;
}
function topEdges(W,n,k){
  const e=[];
  for(let i=0;i<n;i++) for(let j=i+1;j<n;j++){ const w=W[i*n+j]; if(w>0) e.push([i,j,w]); }
  e.sort((a,b)=>b[2]-a[2]);
  const keep=e.slice(0,k), mx=keep.length?keep[0][2]:1;
  return keep.map(([i,j,w])=>[i,j,Math.round(w/mx*1000)/1000]);
}

$("upfile").onchange = async (ev) => {
  const f = ev.target.files[0]; if(!f) return;
  const note=$("upnote");
  try{
    UPROWS = parseMatrix(await f.text()); UPN = UPROWS.length;
    let best=0,bs=-1;
    for(let i=0;i<UPN;i++){ let s=0; for(let j=0;j<UPN;j++) s+=UPROWS[i][j]; if(s>bs){bs=s;best=i;} }
    $("uponset").innerHTML = UPROWS.map((_,i)=>
      `<option value="${i}" ${i===best?"selected":""}>Region ${i+1}${i===best?" (most connected)":""}</option>`).join("");
    $("uplabel").textContent = f.name;
    $("upstat").textContent = `${UPN} regions`;
    $("upform").hidden=false; note.hidden=true;
  }catch(err){
    note.hidden=false; note.textContent = "Could not read that file: " + err.message;
    $("upform").hidden=true; $("upstat").textContent="";
  }
};

$("uprun").onclick = () => {
  if(!UPROWS) return;
  const onset = +$("uponset").value, K = -5.0;
  const names = Object.keys(DRUGS);
  /* 15 calibration probes, then one top check plus three NWORK-wide rounds for the boundary and for
     each release drug, then one recruitment mask per drug plus the untreated baseline. More total
     simulations than a bisection, far fewer sequential waits. */
  const ROUND = 1 + 3*NWORK;
  const TOTAL = 15 + ROUND + names.filter(k=>DRUGS[k].ce>0).length*ROUND + names.length + 1;
  let done = 0;
  const bump = () => { done++; $("upbar").style.width = Math.min(100, done/TOTAL*100) + "%";
    $("upmsg").textContent = `${done} / ~${TOTAL} simulations`; };

  $("uprun").disabled = true; $("upprog").hidden = false; $("upnote").hidden = true;
  $("upmsg").textContent = "starting " + NWORK + " workers…";

  const url = URL.createObjectURL(new Blob([WORKER_SRC], {type:"text/javascript"}));
  const pool = Array.from({length:NWORK}, ()=> new Worker(url));
  let Winfo = null, ready = 0;

  const call = (w, msg) => new Promise(res=>{
    const h = (e)=>{ if(e.data.tick){ bump(); return; } w.removeEventListener("message",h); res(e.data); };
    w.addEventListener("message", h); w.postMessage(msg);
  });

  (async () => {
   try{
    const infos = await Promise.all(pool.map(w=>call(w,{cmd:"init", rows:UPROWS})));
    Winfo = infos[0];
    /* Run a queue of probe jobs across every worker and return their verdicts in order. */
    const fan = async (jobs) => {
      const out = new Array(jobs.length);
      const q = jobs.map((j,i)=>[j,i]);
      await Promise.all(pool.map(async w=>{
        while(q.length){ const [j,i]=q.shift(); out[i] = (await call(w, j)).spreads; }
      }));
      return out;
    };

    /* CALIBRATION, in parallel. Spreading is monotone in background excitability, and Python takes the
       lowest grid value that spreads, so evaluating the whole 15-point grid at once gives exactly the
       Python answer in two waves instead of walking a bisection five deep. */
    const grid = [];
    for(let v=-2.20; v<=-2.0585; v+=0.01) grid.push(Math.round(v*1000)/1000);
    const calRes = await fan(grid.map((v,i)=>({cmd:"probe", onset, shift:0, scale:1, x0h:v, K, i})));
    const ci = calRes.indexOf(true);
    if(ci < 0) throw new Error("this connectome never spreads a seizure across the tested excitability range, so there is nothing for a drug to stop");
    const cal = {ok:true, x0Healthy:grid[ci]};

    /* The smallest value in (lo,hi] at which the seizure stops spreading, found by narrowing a grid
       instead of bisecting. Each round divides the bracket by NWORK+1 using one parallel wave, so
       three rounds on 8 workers resolve to hi/729 in three waits where a bisection of equal precision
       would take ten. `mk(v)` turns a candidate value into a probe job. */
    const searchSmallest = async (lo, hi, rounds, mk) => {
      if((await call(pool[0], mk(hi))).spreads) return null;       // never contains, even at the top
      let a=lo, b=hi;
      for(let r=0;r<rounds;r++){
        const pts = Array.from({length:NWORK},(_,i)=> a + (b-a)*(i+1)/(NWORK+1));
        const res = await fan(pts.map(v=>mk(v)));
        let k = pts.length;
        for(let i=0;i<pts.length;i++) if(!res[i]){ k=i; break; }
        if(k < pts.length) b = pts[k];
        if(k > 0) a = pts[k-1];
      }
      return (a+b)/2;
    };

    const X0U = 0.004, CU = 0.06, SMAXU = 4;
    const bnd = await searchSmallest(0, SMAXU*2.541*X0U, 3,
      v => ({cmd:"probe", onset, shift:-v, scale:1, x0h:cal.x0Healthy, K}));

    const th = {};
    const mixed = names.filter(k=>DRUGS[k].ce>0);
    const pure  = names.filter(k=>DRUGS[k].ce===0);
    /* Drugs with no coupling effect need no search at all: the project's verified identity gives
       threshold = boundary / (x0_effect x X0_UNIT) exactly. Six of the eight land for free. */
    for(const k of pure) th[k] = bnd===null ? null : Math.round(bnd/(DRUGS[k].xe*X0U)*1000)/1000;

    /* The two release drugs move both dials, so they do need a search - but a parallel one. */
    for(const k of mixed){
      const xe=DRUGS[k].xe, ce=DRUGS[k].ce;
      const s = await searchSmallest(0, SMAXU, 3,
        v => ({cmd:"probe", onset, shift:-v*xe*X0U, scale:Math.max(0,1-v*ce*CU), x0h:cal.x0Healthy, K}));
      th[k] = s===null ? null : Math.round(s*1000)/1000;
    }
    for(const k of names) if(th[k]!==null && th[k]!==undefined && th[k] > 4) th[k]=null;

    const masks = {};
    const mq = names.slice();
    await Promise.all(pool.map(async w=>{
      while(mq.length){ const k=mq.shift();
        const r = await call(w,{cmd:"mask", onset, x0h:cal.x0Healthy, K, xe:DRUGS[k].xe, ce:DRUGS[k].ce, s:1.0, drug:k});
        masks[k]=r.mask; }
    }));
    masks.baseline = (await call(pool[0],{cmd:"mask", onset, x0h:cal.x0Healthy, K, xe:0, ce:0, s:0, drug:"baseline"})).mask;

    const W = Winfo.W, n = Winfo.n;
    const regions = Array.from({length:n},(_,i)=>({i, name:"Region "+(i+1), lobe:"other", h:""}));
    DATA.push({ id:"upload", name:"Your scan", onset, onsetName:"Region "+(onset+1),
                regions, masks, th, ranking: Object.keys(th).sort((a,b)=>(th[a]??99)-(th[b]??99)) });
    GEO.push({ pos:fibSphere(n), names:regions.map(r=>r.name), lobes:new Array(n).fill(null),
               edges:topEdges(W,n,Math.min(400,n*5)), hemi:null, layout:"upload", nEdges:Winfo.edges });

    pool.forEach(w=>w.terminate()); URL.revokeObjectURL(url);
    $("upprog").hidden = true; $("uprun").disabled = false;
    $("upnote").hidden = false;
    $("upnote").innerHTML = `Done. Calibrated at x<sub>0</sub> = ${cal.x0Healthy}; added as <b>Your scan</b> in the
      connectome selector. Same equations as the Python pipeline, but a different noise seed, so a threshold
      can sit one bisection step (~0.06) from a full run.`;
    S.scan = DATA.length-1; S.untreated = false;
    if(!drugsFor().includes(S.drug)) S.drug = sortedDrugs()[0];
    chrome(); render();
   }catch(err){
    pool.forEach(w=>w.terminate()); URL.revokeObjectURL(url);
    $("upprog").hidden = true; $("uprun").disabled = false;
    $("upnote").hidden = false; $("upnote").textContent = "Run failed: " + err.message;
   }
  })();
};
"""


def main():
    d = open(DASH, encoding="utf-8").read()

    # the worker needs each drug's two model numbers on the client side
    assert '"cls": ' in d or 'cls' in d
    d = d.replace("const S={scan:0,", "const UPREADY=true;\nconst S={scan:0,", 1)

    # css
    assert ".addon{color:#6fa8ff}" in d
    d = d.replace(".addon{color:#6fa8ff}", ".addon{color:#6fa8ff}\n" + CSS.strip(), 1)

    # replace the static explainer card with the real thing
    start = d.index('    <div class="blk">\n      <h2>Use your own scan</h2>')
    end = d.index("</div>", d.index("arrived exactly that way", start)) + len("</div>")
    d = d[:start] + MARKUP.strip("\n") + d[end:]

    js = JS.replace("__WORKER_SRC__", repr_js(worker_src + WORKER_TAIL))
    anchor = "function render(){"
    assert anchor in d
    d = d.replace(anchor, js + "\n" + anchor, 1)

    open(DASH, "w", encoding="utf-8").write(d)
    print("injected; dash.html now", len(d) // 1024, "KiB")


def repr_js(s):
    """Embed arbitrary JS source as a JS string literal."""
    return "`" + s.replace("\\", "\\\\").replace("`", "\\`").replace("${", "\\${") + "`"


if __name__ == "__main__":
    main()
