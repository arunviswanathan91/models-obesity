// Interactive equation explainers for the PDAC Obesity Atlas.
// Pure arithmetic drawn to SVG: no Pyodide, no network, no build step.
// Colours come from the site's own CSS custom properties, with fallbacks
// so the module also renders standalone.

const NS = 'http://www.w3.org/2000/svg';
const BLUE = 'var(--blue,#2C5AA0)';
const MAGENTA = 'var(--magenta,#D81B60)';
const TEAL = 'var(--teal,#0E7C7B)';
const GREY = 'var(--grey,#6F7378)';
const LINE = 'var(--line,#DCDDE0)';
const INK = 'var(--ink,#1F1F1F)';

/* ---------- small helpers ---------- */
const s = (tag, attrs = {}) => {
  const e = document.createElementNS(NS, tag);
  for (const [k, v] of Object.entries(attrs)) e.setAttribute(k, v);
  return e;
};
const h = (tag, text, attrs = {}) => {
  const e = document.createElement(tag);
  if (text !== undefined) e.textContent = text;
  for (const [k, v] of Object.entries(attrs)) e.setAttribute(k, v);
  return e;
};
const minus = x => (x < 0 ? '−' : '') + Math.abs(x).toFixed(2);
const minus3 = x => (x < 0 ? '−' : '') + Math.abs(x).toFixed(3);
const clamp = (x, a, b) => Math.max(a, Math.min(b, x));

// Inverse standard normal CDF (Acklam), used for the interval-widening factor.
function probit(p) {
  const a = [-39.69683028665376, 220.9460984245205, -275.9285104469687,
             138.3577518672690, -30.66479806614716, 2.506628277459239];
  const b = [-54.47609879822406, 161.5858368580409, -155.6989798598866,
             66.80131188771972, -13.28068155288572];
  const c = [-0.007784894002430293, -0.3223964580411365, -2.400758277161838,
             -2.549732539343734, 4.374664141464968, 2.938163982698783];
  const d = [0.007784695709041462, 0.3224671290700398, 2.445134137142996,
             3.754408661907416];
  const pl = 0.02425;
  let q, r;
  if (p < pl) {
    q = Math.sqrt(-2 * Math.log(p));
    return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
           ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
  }
  if (p > 1 - pl) return -probit(1 - p);
  q = p - 0.5; r = q * q;
  return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q /
         (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1);
}
const widenFactor = cov => 1.959964 / probit((1 + cov) / 2);

// Deterministic standard normals, so points never jitter while a slider moves.
function normals(n, seed) {
  const out = [];
  let x = seed;
  const u = () => { x = (x * 1664525 + 1013904223) % 4294967296; return (x + 0.5) / 4294967296; };
  for (let i = 0; i < n; i++) {
    const r = Math.sqrt(-2 * Math.log(u()));
    out.push(r * Math.cos(2 * Math.PI * u()));
  }
  return out;
}

/* ---------- illustrative data ----------
   Replace DEMO_FOREST / DEMO_GRID with rows from atlas_web_rows to show real
   posterior means. Shapes: {label, est, se} and {cell, sig, est, se}. */
const DEMO_FOREST = [
  { label: 'Macrophage lipid handling',      est:  0.34, se: 0.11 },
  { label: 'Macrophage iron handling',       est:  0.29, se: 0.12 },
  { label: 'Macrophage tissue repair',       est:  0.22, se: 0.14 },
  { label: 'CD4+ helper programme',          est:  0.26, se: 0.10 },
  { label: 'CD4+ regulatory programme',      est:  0.09, se: 0.13 },
  { label: 'CD8+ metabolic programme',       est:  0.24, se: 0.16 },
  { label: 'CD8+ cytotoxicity',              est:  0.07, se: 0.18 },
  { label: 'CD8+ exhaustion',                est: -0.05, se: 0.15 },
  { label: 'Fibroblast ECM assembly',        est: -0.18, se: 0.13 },
  { label: 'Fibroblast contractility',       est: -0.11, se: 0.17 },
  { label: 'Epithelial proliferation',       est:  0.04, se: 0.12 },
  { label: 'Endothelial activation',         est:  0.13, se: 0.19 }
];
const CELLS = ['Macrophage', 'CD4+ T', 'CD8+ T', 'B cell', 'Fibroblast', 'Epithelial'];
const SIGS  = ['Lipid', 'Iron', 'Repair', 'Metabolic', 'Cytotoxic', 'ECM'];
const DEMO_GRID = (() => {
  const z = normals(CELLS.length * SIGS.length, 77003);
  const out = [];
  let i = 0;
  for (const cell of CELLS) for (const sig of SIGS) {
    out.push({ cell, sig, est: Math.round(z[i] * 22) / 100, se: 0.10 + (i % 5) * 0.025 });
    i++;
  }
  return out;
})();

/* ---------- shared shrinkage maths ---------- */
// Precision-weighted posterior under a Normal(mean, tau^2) prior on the slope.
function shrink(est, se, mean, tau) {
  const pd = 1 / (se * se), pp = 1 / (tau * tau + 1e-9);
  const sd = Math.sqrt(1 / (pd + pp));
  return { est: (est * pd + mean * pp) / (pd + pp), sd };
}
const precisionMean = rows =>
  rows.reduce((a, r) => a + r.est / (r.se * r.se), 0) /
  rows.reduce((a, r) => a + 1 / (r.se * r.se), 0);

/* ---------- widget 1: parent-protein conditioning ---------- */
const conditioning = {
  id: 'conditioning',
  title: 'Parent-protein conditioning',
  eq: '<span data-term="bMP">β<sub>M|P</sub></span> = <span data-term="bM">β<sub>M</sub></span> − <span data-term="kappa">κ</span><span data-term="bP">β<sub>P</sub></span>, <span data-term="kappa">κ = ρ σ<sub>M</sub>/σ<sub>P</sub></span>',
  terms: {
    bMP: ['Conditional modification slope', 'The BMI association of a modification at a fixed amount of its parent protein. This is the conditional call in the PTM screen.'],
    bM: ['Marginal modification slope', 'The BMI slope of the modification on its own, ignoring what the parent protein does.'],
    bP: ['Parent-protein slope', 'The BMI slope of the protein that carries the modification.'],
    kappa: ['Coupling κ', 'How much of the protein slope is carried into the modification: the residual correlation ρ times the ratio of the two scales. With Student-t errors the conditional slope has ν + 1 degrees of freedom.']
  },
  ref: 'Equations 9–10',
  blurb: 'The BMI association of a modification at fixed parent protein is not the modification slope minus the protein slope. κ reweights the protein slope by the residual correlation and the ratio of the two Student-t scales. The presets are the paper’s own calibration scenarios.',
  note: 'Estimand: expected within-plex change in standardized modification per 5 kg/m² at a fixed standardized parent protein. Not occupancy, not a concentration ratio, not a causal direct effect. κ is recomputed within each joint posterior draw, preserving the dependence among the two slopes, their scales and their correlation. The scale matrix is not the covariance matrix, but the common ν/(ν−2) multiplier cancels in the ratio. The conditional distribution is Student-t with ν + 1 degrees of freedom: its mean is linear in p, so βₘ|ₚ is constant, while its scale is multiplied by √((ν + dₚ²)/(ν + 1)) — the flare in the shaded band. Of 367 eligible phosphofeatures and 1,281 glycofeatures, 14 and 67 became conditional candidates.',
  params: [
    { k: 'bm',    label: 'β<sub>M</sub> modification slope', min: -1,    max: 1,   step: 0.001, value: 0.341, colour: MAGENTA, fmt: x => minus3(x), hint: 'Marginal BMI association of the modified peptide.' },
    { k: 'bp',    label: 'β<sub>P</sub> parent protein slope', min: -1,  max: 1,   step: 0.001, value: 0.25, colour: BLUE, fmt: x => minus3(x), hint: 'BMI association of the unmodified parent protein.' },
    { k: 'rho',   label: 'ρ residual correlation', min: -0.95, max: 0.95, step: 0.01, value: 0.50, colour: TEAL, hint: 'How the two layers co-vary after covariates and plex.' },
    { k: 'ratio', label: 'σ<sub>M</sub>/σ<sub>P</sub> scale ratio', min: 0.2, max: 3, step: 0.005, value: 1.125, colour: TEAL, fmt: x => x.toFixed(3), hint: 'Student-t scale of the modification over the protein. The calibration templates fix these at 0.9 and 0.8.' },
    { k: 'nu',    label: 'ν degrees of freedom', min: 2.1, max: 30, step: 0.1, value: 8, colour: GREY, fmt: x => x.toFixed(1), hint: 'Sets how far the conditional spread flares away from the mean protein level.' }
  ],
  presets: [
    { label: 'Null', v: { bm: 0,     bp: 0,    rho: 0.5, ratio: 1.125, nu: 8 } },
    { label: 'Protein-only', v: { bm: 0.197, bp: 0.25, rho: 0.7, ratio: 1.125, nu: 8 } },
    { label: 'Conditional signal', v: { bm: 0.341, bp: 0.25, rho: 0.5, ratio: 1.125, nu: 8 } },
    { label: 'Heavy-tailed', v: { bm: 0.397, bp: 0.25, rho: 0.7, ratio: 1.125, nu: 3 } },
    { label: 'Sign reversal', v: { bm: 0.20, bp: -0.50, rho: 0.70, ratio: 1.60, nu: 8 } }
  ],
  usesData: true,
  real() {
    if (!DATA) return null;
    const groups = {};
    for (const r of DATA.ptm.rows) {
      const [gene, assay, tr, cand, bm, bp, cond, kappa, rho] = r;
      const ratio = rho ? kappa / rho : 1;
      (groups[tr] = groups[tr] || []).push({
        label: `${gene} ${assay}  \u03B2M ${minus(bm)} \u2192 ${minus(cond)}${cand ? '  \u2605' : ''}`,
        v: { bm, bp, rho, ratio: Math.abs(ratio) < 0.2 ? 0.2 : Math.min(3, Math.abs(ratio)), nu: 8 }
      });
    }
    return Object.entries(groups).map(([g, items]) => ({ group: g, items }));
  },
  init() { const z1 = normals(170, 20260109), z2 = normals(170, 91823); return { z1, z2 }; },
  draw(svg, v, mem) {
    const k = v.rho * v.ratio, cond = v.bm - k * v.bp, naive = v.bm - v.bp;
    svg.setAttribute('viewBox', '0 0 340 260');
    svg.replaceChildren();
    const X = p => 52 + ((clamp(p, -3, 3) + 3) / 6) * 268;
    const mMax = Math.min(6, Math.max(2, 3.2 * v.ratio));
    const Y = m => 126 - (clamp(m, -mMax, mMax) / mMax) * 100;

    for (const f of [0, 0.5, 1]) {
      svg.append(s('line', { x1: 52, x2: 320, y1: 26 + f * 200, y2: 26 + f * 200,
        stroke: LINE, 'stroke-dasharray': f === 0.5 ? '' : '3 3' }));
    }
    svg.append(s('line', { x1: X(0), x2: X(0), y1: 26, y2: 226, stroke: LINE, 'stroke-dasharray': '3 3' }));
    const sd = Math.sqrt(Math.max(0, 1 - v.rho * v.rho));
    // conditional scale flares away from the mean protein level: sqrt((nu + dP^2)/(nu + 1))
    const condSD = p => v.ratio * sd * Math.sqrt((v.nu + p * p) / (v.nu + 1));
    let up = '', dn = '';
    for (let i = 0; i <= 60; i++) {
      const p = -3 + (i / 60) * 6;
      up += (i ? 'L' : 'M') + X(p) + ' ' + Y(k * p + 1.96 * condSD(p));
      dn = 'L' + X(p) + ' ' + Y(k * p - 1.96 * condSD(p)) + dn;
    }
    svg.append(s('path', { d: up + dn.replace('L', 'L') + 'Z', fill: TEAL, opacity: 0.1 }));
    for (let i = 0; i < mem.z1.length; i++) {
      svg.append(s('circle', { cx: X(mem.z1[i]), cy: Y(v.ratio * (v.rho * mem.z1[i] + sd * mem.z2[i])),
        r: 2.6, fill: GREY, opacity: 0.45 }));
    }
    svg.append(s('line', { x1: X(-3), y1: Y(k * -3), x2: X(3), y2: Y(k * 3),
      stroke: TEAL, 'stroke-width': 2.4 }));
    svg.append(s('line', { x1: 52, x2: 320, y1: 226, y2: 226, stroke: INK, 'stroke-width': 1 }));
    svg.append(s('line', { x1: 52, x2: 52, y1: 26, y2: 226, stroke: INK, 'stroke-width': 1 }));
    const lab = (x, y, t, anchor = 'middle', fill = GREY) => {
      const n = s('text', { x, y, 'text-anchor': anchor, fill, 'font-size': 11 });
      n.textContent = t; return n;
    };
    svg.append(lab(186, 250, 'standardized parent protein'));
    svg.append(lab(52, 240, '−3'), lab(320, 240, '3'));
    svg.append(lab(47, 30, mMax.toFixed(1), 'end'), lab(47, 130, '0', 'end'), lab(47, 226, '−' + mMax.toFixed(1), 'end'));
    const yl = lab(16, 126, 'standardized modification');
    yl.setAttribute('transform', 'rotate(-90 16 126)');
    svg.append(yl);
    const slope = lab(X(2.2), Y(k * 2.2) - 8, 'κ = ' + minus(k), 'middle', TEAL);
    slope.setAttribute('font-weight', '700');
    svg.append(slope);

    let verdict, T = 0.15;
    const mHit = Math.abs(v.bm) >= T, cHit = Math.abs(cond) >= T;
    if (mHit && cHit && cond * v.bm < 0) verdict = ['Sign flip', 'Marginal and conditioned slopes point opposite ways; the apparent direction came from the parent protein.'];
    else if (mHit && !cHit) verdict = ['Dropped by conditioning', 'The modification tracked total protein abundance. At fixed parent protein the signal is gone.'];
    else if (!mHit && cHit) verdict = ['Revealed by conditioning', 'No marginal signal, but holding the parent protein fixed uncovers one. A marginal screen misses these.'];
    else if (mHit && cHit) verdict = ['Survives conditioning', 'Associated with BMI marginally and at fixed parent protein — the strongest class of candidate.'];
    else verdict = ['Not a candidate', 'No support either way at this illustrative cutoff.'];

    return {
      eq: `β<sub>M|P</sub> = ${minus3(v.bm)} − (${minus3(k)} × ${minus3(v.bp)}) = <b>${minus3(cond)}</b>`,
      cards: [
        ['κ', minus3(k), 'ρ × scale ratio'],
        ['β<sub>M|P</sub> conditioned', minus3(cond), 'what the paper reports', true],
        ['β<sub>M</sub> − β<sub>P</sub> naive', minus3(naive), 'what a reader assumes'],
        ['discrepancy', minus3(cond - naive), 'conditioned minus naive']
      ],
      verdict
    };
  }
};

/* ---------- widget 2: partial pooling, forest of 95% HDIs ---------- */
const pooling = {
  id: 'pooling',
  title: 'Partial pooling and the 95% HDI',
  eq: '<span data-term="beta">β<sub>j</sub></span> ~ Normal(<span data-term="m">μ<sub>c(j)</sub></span>, <span data-term="tau">τ<sup>2</sup></span>)',
  terms: {
    beta: ['Signature slope', 'The BMI slope of one signature.'],
    m: ['Cell-type mean', 'The average BMI slope of all signatures in that cell type. This is the cell-type estimate reported in Figure 3b.'],
    tau: ['Spread between signatures', 'How much signature slopes in one cell type are allowed to differ. Small τ pulls each slope towards the cell-type mean; prior HalfNormal(0.30).']
  },
  ref: 'Equation 4, decision rule',
  blurb: 'Signature slopes are shrunk towards their cell-type mean. Raise τ and the model trusts each signature on its own; lower it and weakly supported signatures are pulled in. The widening slider applies the paper’s correction for simulated undercoverage.',
  note: 'An HDI excluding zero is directional support under the model, with no multiplicity control and no minimum effect size. Simulated residual correlation among signatures gave about 84.6% coverage, so named associations were re-checked after widening each interval by 1.960/1.426 ≈ 1.37. Rows are 12 published slopes spread from the strongest to the weakest survivor; the noise multiplier is a what-if, because the published SEs are all small relative to τ and shrinkage is nearly invisible at ×1.',
  params: [
    { k: 'tau',  label: 'τ between-signature SD', min: 0.02, max: 1.2, step: 0.01, value: 0.10, colour: TEAL, hint: 'Prior spread of slopes. Small τ pulls every slope to the mean; large τ leaves them alone. The paper’s prior is HalfNormal(0.30).' },
    { k: 'mult', label: 'noisier data (× SE)', min: 1, max: 8, step: 0.25, value: 3, colour: BLUE, fmt: x => '×' + x.toFixed(2).replace(/0$/, ''), hint: 'What-if: inflate each published SE, as if fewer cells supported the signature. ×1 is the data as published.' },
    { k: 'cov',  label: 'achieved coverage', min: 0.70, max: 0.95, step: 0.005, value: 0.846, colour: MAGENTA, fmt: x => (x * 100).toFixed(1) + '%', hint: 'Simulated coverage of a nominal 95% interval. 84.6% is the paper’s value.' },
    { k: 'show', label: 'show unpooled estimates', min: 0, max: 1, step: 1, value: 1, colour: GREY, fmt: x => (x ? 'yes' : 'no'), hint: 'Overlay the raw per-signature estimates for comparison.' }
  ],
  presets: [
    { label: 'Published as is', v: { tau: 0.30, mult: 1, cov: 0.846, show: 1 } },
    { label: 'Sparse data, pooled', v: { tau: 0.10, mult: 3, cov: 0.846, show: 1 } },
    { label: 'Sparse, no pooling', v: { tau: 1.2, mult: 3, cov: 0.846, show: 1 } },
    { label: 'Sparse, heavy pooling', v: { tau: 0.03, mult: 4, cov: 0.846, show: 1 } },
    { label: 'Nominal coverage', v: { tau: 0.10, mult: 3, cov: 0.95, show: 1 } }
  ],
  usesData: true,
  init() { return {}; },
  draw(svg, v) {
    const W = DATA && DATA.widening && DATA.widening.rows;
    const pick = W && W.length ? Array.from({ length: 12 }, (_, i) => W[Math.round(i * (W.length - 1) / 11)]) : null;
    const realRows = pick ? pick.map(r => ({
          label: r[0] + ' · ' + r[1].replace(/_/g, ' '),
          est: r[2], se: (r[4] - r[3]) / 3.9199 * v.mult })) : null;
    const rows = realRows || DEMO_FOREST.map(r => ({ ...r, se: r.se * v.mult / 3 })), mean = precisionMean(rows), f = widenFactor(v.cov);
    const n = rows.length, rowH = 21, top = 30, H = top + n * rowH + 44;
    svg.setAttribute('viewBox', `0 0 520 ${H}`);
    svg.replaceChildren();
    const ext = Math.max(...rows.map(r => Math.abs(r.est) + 1.96 * r.se * Math.max(1, f)));
    const lim = [0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.2, 1.6, 2.4, 3.2].find(q => q >= ext * 1.02) || 3.2, L = 240, R = 506;
    const X = b => L + ((clamp(b, -lim, lim) + lim) / (2 * lim)) * (R - L);

    for (const t of [-lim, -lim / 2, 0, lim / 2, lim]) {
      svg.append(s('line', { x1: X(t), x2: X(t), y1: top - 10, y2: top + n * rowH,
        stroke: t === 0 ? INK : LINE, 'stroke-dasharray': t === 0 ? '' : '3 3' }));
      const tl = s('text', { x: X(t), y: top + n * rowH + 16, 'text-anchor': 'middle', fill: GREY, 'font-size': 11 });
      tl.textContent = minus(+t.toFixed(2)); svg.append(tl);
    }
    svg.append(s('line', { x1: X(mean), x2: X(mean), y1: top - 10, y2: top + n * rowH, stroke: TEAL, 'stroke-width': 1.5, opacity: 0.6 }));
    const ml = s('text', { x: X(mean), y: top - 15, 'text-anchor': 'middle', fill: TEAL, 'font-size': 11, 'font-weight': 700 });
    ml.textContent = 'cell-type mean'; svg.append(ml);

    let called = 0, calledWide = 0;
    rows.forEach((r, i) => {
      const y = top + i * rowH + rowH / 2;
      const p = shrink(r.est, r.se, mean, v.tau);
      const half = 1.959964 * p.sd, wide = half * f;
      const excl = Math.abs(p.est) > half, exclW = Math.abs(p.est) > wide;
      if (excl) called++;
      if (exclW) calledWide++;
      const col = p.est >= 0 ? MAGENTA : BLUE;

      const nm = s('text', { x: 232, y: y + 4, 'text-anchor': 'end', fill: excl ? INK : GREY, 'font-size': 11 });
      nm.textContent = r.label; svg.append(nm);

      if (v.show) {
        svg.append(s('line', { x1: X(r.est - 1.96 * r.se), x2: X(r.est + 1.96 * r.se), y1: y - 6, y2: y - 6, stroke: GREY, 'stroke-width': 1.6, opacity: 0.6 }));
        svg.append(s('circle', { cx: X(r.est), cy: y - 6, r: 3, fill: GREY, opacity: 0.85 }));
      }
      svg.append(s('line', { x1: X(p.est - wide), x2: X(p.est + wide), y1: y + 3, y2: y + 3, stroke: col, 'stroke-width': 1, opacity: 0.35 }));
      svg.append(s('line', { x1: X(p.est - half), x2: X(p.est + half), y1: y + 3, y2: y + 3, stroke: col, 'stroke-width': 2.2 }));
      svg.append(s('circle', { cx: X(p.est), cy: y + 3, r: 3.4, fill: col }));
    });
    const cap = s('text', { x: (L + R) / 2, y: H - 10, 'text-anchor': 'middle', fill: GREY, 'font-size': 11 });
    cap.textContent = 'BMI slope per 5 kg/m² (score SD)'; svg.append(cap);

    return {
      eq: `τ = ${v.tau.toFixed(2)} · SE ×${v.mult.toFixed(2).replace(/0$/, '')}${v.mult === 1 ? ' (as published)' : ' (what-if)'} · widening × <b>${f.toFixed(2)}</b> at ${(v.cov * 100).toFixed(1)}% coverage`,
      cards: [
        ['called, nominal 95%', String(called) + ' / ' + n, 'HDI excludes zero'],
        ['called, widened', String(calledWide) + ' / ' + n, 'after ×' + f.toFixed(2), true],
        ['lost to widening', String(called - calledWide), 'sensitive to coverage'],
        ['cell-type mean', minus(mean), 'precision-weighted']
      ],
      verdict: calledWide < called
        ? ['Widening changes the call', (called - calledWide) + ' of ' + called + ' associations no longer exclude zero once the interval is widened for the simulated undercoverage.']
        : ['Stable under widening', 'Every association that excluded zero still does after the correction.']
    };
  }
};

/* ---------- widget 3: Student-t and outlier influence ---------- */
const studentT = {
  id: 'student-t',
  title: 'Why a Student-t likelihood',
  eq: '<span data-term="w">w(r)</span> = (<span data-term="nu">ν</span> + 1) / (<span data-term="nu">ν</span> + <span data-term="r">r<sup>2</sup></span>)',
  terms: {
    w: ['Weight', 'How much a residual counts in the fit. A Gaussian likelihood gives every residual a weight of one.'],
    nu: ['Degrees of freedom ν', 'Tail heaviness. Low ν downweights large residuals strongly; as ν grows the curve flattens to the Gaussian. The model gives ν − 2 an Exponential(0.10) prior.'],
    r: ['Standardized residual', 'How far a patient sits from the fitted value, in units of the residual scale.']
  },
  ref: 'Likelihood, all three models',
  blurb: 'A Gaussian likelihood gives every residual the same weight, so one aberrant patient moves the slope. The Student-t downweights large residuals automatically. The weight a residual receives is the curve below; ν controls how fast it falls away.',
  note: 'All three models use a Student-t likelihood with ν − 2 given an Exponential(0.10) prior, estimated per cell type in the coarse immune model and shared in the fine immune and nonimmune models. As ν grows the weight curve flattens towards the Gaussian, where every residual counts equally.',
  params: [
    { k: 'nu', label: 'ν degrees of freedom', min: 2.1, max: 50, step: 0.1, value: 4, colour: TEAL, fmt: x => x.toFixed(1), hint: 'Low ν means heavy tails and strong downweighting.' },
    { k: 'r',  label: 'residual of one patient (SD)', min: 0, max: 6, step: 0.05, value: 3, colour: MAGENTA, fmt: x => x.toFixed(2), hint: 'How far that patient sits from the fitted value.' }
  ],
  presets: [
    { label: 'Generating \u03BD 8', v: { nu: 8, r: 3 } },
    { label: 'Heavy tails \u03BD 3', v: { nu: 3, r: 3 } },
    { label: 'Near-Gaussian', v: { nu: 30, r: 3 } },
    { label: 'Extreme outlier', v: { nu: 4, r: 5.5 } },
    { label: 'Ordinary residual', v: { nu: 4, r: 0.8 } }
  ],
  init() { return {}; },
  draw(svg, v) {
    svg.setAttribute('viewBox', '0 0 340 250');
    svg.replaceChildren();
    const L = 48, R = 322, T = 22, B = 206;
    const X = r => L + (r / 6) * (R - L);
    const Y = w => B - w * (B - T);
    const w = r => (v.nu + 1) / (v.nu + r * r);

    for (const g of [0, 0.25, 0.5, 0.75, 1]) {
      svg.append(s('line', { x1: L, x2: R, y1: Y(g), y2: Y(g), stroke: LINE, 'stroke-dasharray': g ? '3 3' : '' }));
      const t = s('text', { x: L - 6, y: Y(g) + 4, 'text-anchor': 'end', fill: GREY, 'font-size': 11 });
      t.textContent = g.toFixed(2); svg.append(t);
    }
    svg.append(s('line', { x1: L, x2: R, y1: Y(1), y2: Y(1), stroke: BLUE, 'stroke-width': 2, 'stroke-dasharray': '6 4' }));
    const gl = s('text', { x: R - 4, y: Y(1) - 7, 'text-anchor': 'end', fill: BLUE, 'font-size': 11, 'font-weight': 700 });
    gl.textContent = 'Gaussian: every residual weight 1'; svg.append(gl);

    let d = '';
    for (let i = 0; i <= 120; i++) { const r = (i / 120) * 6; d += (i ? 'L' : 'M') + X(r) + ' ' + Y(Math.min(1, w(r))); }
    svg.append(s('path', { d, fill: 'none', stroke: TEAL, 'stroke-width': 2.4 }));

    const wr = Math.min(1, w(v.r));
    svg.append(s('line', { x1: X(v.r), x2: X(v.r), y1: Y(0), y2: Y(wr), stroke: MAGENTA, 'stroke-width': 1.4, 'stroke-dasharray': '4 3' }));
    svg.append(s('circle', { cx: X(v.r), cy: Y(wr), r: 4.5, fill: MAGENTA }));
    const pl = s('text', { x: clamp(X(v.r), L + 30, R - 30), y: Y(wr) - 11, 'text-anchor': 'middle', fill: MAGENTA, 'font-size': 11.5, 'font-weight': 700 });
    pl.textContent = 'weight ' + wr.toFixed(2); svg.append(pl);

    svg.append(s('line', { x1: L, x2: R, y1: B, y2: B, stroke: INK }));
    svg.append(s('line', { x1: L, x2: L, y1: T, y2: B, stroke: INK }));
    for (const r of [0, 2, 4, 6]) {
      const t = s('text', { x: X(r), y: B + 16, 'text-anchor': 'middle', fill: GREY, 'font-size': 11 });
      t.textContent = String(r); svg.append(t);
    }
    const xl = s('text', { x: (L + R) / 2, y: 240, 'text-anchor': 'middle', fill: GREY, 'font-size': 11 });
    xl.textContent = 'residual, in residual SDs'; svg.append(xl);
    const yl = s('text', { x: 14, y: 114, 'text-anchor': 'middle', fill: GREY, 'font-size': 11, transform: 'rotate(-90 14 114)' });
    yl.textContent = 'weight in the fit'; svg.append(yl);

    return {
      eq: `w(${v.r.toFixed(2)}) = (${v.nu.toFixed(1)} + 1) / (${v.nu.toFixed(1)} + ${(v.r * v.r).toFixed(2)}) = <b>${wr.toFixed(3)}</b>`,
      cards: [
        ['weight under t', wr.toFixed(3), 'this patient', true],
        ['weight under Gaussian', '1.000', 'every patient'],
        ['influence removed', ((1 - wr) * 100).toFixed(0) + '%', 'relative to Gaussian'],
        ['ν', v.nu.toFixed(1), v.nu > 30 ? 'effectively Gaussian' : 'heavy-tailed']
      ],
      verdict: v.nu > 30
        ? ['Tails are effectively Gaussian', 'At this ν the weight curve is nearly flat, so the t buys almost nothing over a normal likelihood.']
        : wr < 0.5
          ? ['Outlier substantially downweighted', 'This patient contributes less than half what a Gaussian fit would give them.']
          : ['Residual treated as ordinary', 'Within the bulk of the distribution the t and the Gaussian behave much alike.']
    };
  }
};

/* ---------- widget 4: within-plex vs between-plex ---------- */
const plex = {
  id: 'plex',
  title: 'Within-plex and between-plex BMI',
  eq: '<span data-term="y">y<sub>ip</sub></span> = <span data-term="alpha">α<sub>p</sub></span> + <span data-term="within">β<sup>within</sup>(x<sub>ip</sub> − x̄<sub>p</sub>)</span> + <span data-term="between">β<sup>between</sup>x̄<sub>p</sub></span>',
  terms: {
    y: ['Score', 'The measured value for patient i in plex p.'],
    alpha: ['Plex intercept', 'The baseline level of one measurement plex (batch).'],
    within: ['Within-plex slope', 'BMI association using only differences between patients measured in the same plex, so plex-level batch effects cancel.'],
    between: ['Between-plex slope', 'Association with the plex-average BMI. It mixes the BMI signal with plex-to-plex batch differences, so it is read with more caution.']
  },
  ref: 'Equation 7',
  blurb: 'Patient BMI is split into the deviation from the plex mean and the plex mean itself. The reported coefficient is the within-plex slope, so a difference in average BMI between TMT plexes cannot leak into it. Pull the plexes apart and watch the pooled line tilt while the within-plex lines do not.',
  note: 'The RNA plex variable groups patients by their proteomics plex; it does not assert an RNA batch effect. Plex intercepts are non-centered with HalfNormal(0.3) scales. The same decomposition carries into the protein–PTM models.',
  params: [
    { k: 'within',  label: 'β within-plex slope', min: -0.8, max: 0.8, step: 0.01, value: 0.25, colour: TEAL, hint: 'The association the paper reports.' },
    { k: 'between', label: 'between-plex confounding', min: -1, max: 1, step: 0.01, value: 0.60, colour: MAGENTA, hint: 'How strongly plex intercepts track plex mean BMI.' },
    { k: 'spread',  label: 'spread of plex mean BMI', min: 0, max: 1.6, step: 0.01, value: 0.9, colour: BLUE, hint: 'How far apart the plexes sit on BMI.' }
  ],
  presets: [
    { label: 'Reference', v: { within: 0.25, between: 0.60, spread: 0.9 } },
    { label: 'Plexes aligned', v: { within: 0.25, between: 0.60, spread: 0 } },
    { label: 'Strong confounding', v: { within: 0.25, between: 1.0, spread: 1.4 } },
    { label: 'No plex effect', v: { within: 0.25, between: 0, spread: 0.9 } },
    { label: 'Null within-plex', v: { within: 0, between: 0.8, spread: 1.2 } }
  ],
  init() { return { zx: normals(96, 4242), zy: normals(96, 8181) }; },
  draw(svg, v, mem) {
    svg.setAttribute('viewBox', '0 0 340 250');
    svg.replaceChildren();
    const L = 46, R = 324, T = 20, B = 206;
    const X = x => L + ((clamp(x, -3, 3) + 3) / 6) * (R - L);
    const Y = y => clamp(114 - y * 42, T, B);
    const P = 4, per = 24, cols = [TEAL, BLUE, MAGENTA, '#8A6D1F'];
    const pts = [], centres = [];

    for (let p = 0; p < P; p++) {
      const centre = (p - (P - 1) / 2) * v.spread;
      centres.push(centre);
      for (let i = 0; i < per; i++) {
        const idx = p * per + i;
        const dev = mem.zx[idx] * 0.55;
        const x = centre + dev;
        const y = v.between * centre + v.within * dev + mem.zy[idx] * 0.22;
        pts.push({ x, y, p });
      }
    }
    const mx = pts.reduce((a, q) => a + q.x, 0) / pts.length;
    const my = pts.reduce((a, q) => a + q.y, 0) / pts.length;
    let sxy = 0, sxx = 0;
    for (const q of pts) { sxy += (q.x - mx) * (q.y - my); sxx += (q.x - mx) ** 2; }
    const pooled = sxx ? sxy / sxx : 0;

    for (const g of [-2, 0, 2]) svg.append(s('line', { x1: X(g), x2: X(g), y1: T, y2: B, stroke: LINE, 'stroke-dasharray': '3 3' }));
    svg.append(s('line', { x1: L, x2: R, y1: Y(0), y2: Y(0), stroke: LINE, 'stroke-dasharray': '3 3' }));
    for (const q of pts) svg.append(s('circle', { cx: X(q.x), cy: Y(q.y), r: 2.5, fill: cols[q.p], opacity: 0.5 }));
    for (let p = 0; p < P; p++) {
      const c = centres[p], a = v.between * c, lo = c - 1.1, hi = c + 1.1;
      svg.append(s('line', { x1: X(lo), y1: Y(a + v.within * (lo - c)), x2: X(hi), y2: Y(a + v.within * (hi - c)),
        stroke: cols[p], 'stroke-width': 2.2 }));
    }
    svg.append(s('line', { x1: X(-3), y1: Y(my + pooled * (-3 - mx)), x2: X(3), y2: Y(my + pooled * (3 - mx)),
      stroke: INK, 'stroke-width': 2, 'stroke-dasharray': '7 4' }));
    const pl = s('text', { x: R - 4, y: T + 12, 'text-anchor': 'end', fill: INK, 'font-size': 11, 'font-weight': 700 });
    pl.textContent = 'pooled (ignores plex)'; svg.append(pl);

    svg.append(s('line', { x1: L, x2: R, y1: B, y2: B, stroke: INK }));
    svg.append(s('line', { x1: L, x2: L, y1: T, y2: B, stroke: INK }));
    const xl = s('text', { x: (L + R) / 2, y: 232, 'text-anchor': 'middle', fill: GREY, 'font-size': 11 });
    xl.textContent = 'standardized BMI'; svg.append(xl);
    const yl = s('text', { x: 14, y: 114, 'text-anchor': 'middle', fill: GREY, 'font-size': 11, transform: 'rotate(-90 14 114)' });
    yl.textContent = 'standardized score'; svg.append(yl);
    const cl = s('text', { x: (L + R) / 2, y: 246, 'text-anchor': 'middle', fill: GREY, 'font-size': 11 });
    cl.textContent = 'four TMT plexes, 24 patients each'; svg.append(cl);

    return {
      eq: `β<sup>within</sup> = <b>${minus(v.within)}</b> · pooled slope = ${minus(pooled)}`,
      cards: [
        ['β within-plex', minus(v.within), 'reported', true],
        ['pooled slope', minus(pooled), 'plex ignored'],
        ['bias from pooling', minus(pooled - v.within), 'pooled minus within'],
        ['plex BMI spread', v.spread.toFixed(2), 'SD units']
      ],
      verdict: Math.abs(pooled - v.within) > 0.08
        ? ['Pooling is biased here', 'Ignoring plex gives a slope of ' + minus(pooled) + ' against a true within-plex ' + minus(v.within) + '. The decomposition is what keeps the two apart.']
        : ['Little to separate', 'With plexes this similar in mean BMI, pooled and within-plex slopes nearly agree. The decomposition costs nothing and protects you when they do not.']
    };
  }
};

/* ---------- widget 5: shrinkage across the catalog ---------- */
const heatmap = {
  id: 'heatmap',
  title: 'Shrinkage across the signature catalog',
  eq: '<span data-term="bhat">β̂<sub>j</sub></span> = (<span data-term="b">β<sub>j</sub></span>/<span data-term="s2">σ<sub>j</sub><sup>2</sup></span> + <span data-term="mu">μ<sub>c</sub></span>/<span data-term="tau">τ<sup>2</sup></span>) / (1/<span data-term="s2">σ<sub>j</sub><sup>2</sup></span> + 1/<span data-term="tau">τ<sup>2</sup></span>)',
  terms: {
    bhat: ['Shrunk slope', 'The posterior mean after pooling: a precision-weighted average of the signature’s own estimate and its cell-type mean.'],
    b: ['Own estimate', 'What the data alone say about this signature.'],
    s2: ['Own variance', 'The squared standard error of the signature’s own estimate. A noisy signature has a large σ² and is pulled in further.'],
    mu: ['Cell-type mean', 'The average slope of the cell type that the signature is pulled towards.'],
    tau: ['Spread between signatures', 'How far slopes in one cell type may differ. Small τ means strong pooling.']
  },
  ref: 'Equation 4',
  blurb: 'The same pooling, seen across every cell type at once. Lower τ and the grid cools towards each cell type’s own mean; raise it and each signature keeps its own estimate. Cells outlined in white are those whose 95% HDI excludes zero.',
  note: 'A cell-type mean summarises the catalog signatures for that cell type. It measures neither the abundance of those cells nor their total activation. Values shown are illustrative; wire the grid to atlas_web_rows to display real posterior means.',
  params: [
    { k: 'tau', label: 'τ between-signature SD', min: 0.02, max: 1.2, step: 0.01, value: 0.30, colour: TEAL, hint: 'The paper uses HalfNormal(0.30).' },
    { k: 'lim', label: 'colour scale limit', min: 0.1, max: 0.8, step: 0.01, value: 0.4, colour: GREY, hint: 'Slope at which the scale saturates.' }
  ],
  presets: [
    { label: 'Reference \u03C4 0.30', v: { tau: 0.30, lim: 0.4 } },
    { label: 'No pooling', v: { tau: 1.2, lim: 0.4 } },
    { label: 'Heavy pooling', v: { tau: 0.05, lim: 0.4 } },
    { label: 'Wide colour scale', v: { tau: 0.30, lim: 0.8 } }
  ],
  init() { return {}; },
  draw(svg, v) {
    const cw = 56, ch = 28, L = 96, T = 34;
    const W = L + SIGS.length * cw + 14, H = T + CELLS.length * ch + 52;
    svg.setAttribute('viewBox', `0 0 ${Math.max(W, 340)} ${H}`);
    svg.replaceChildren();
    const colour = b => {
      const t = clamp(Math.abs(b) / v.lim, 0, 1);
      return b >= 0 ? `rgba(216,27,96,${(0.12 + 0.82 * t).toFixed(3)})`
                    : `rgba(44,90,160,${(0.12 + 0.82 * t).toFixed(3)})`;
    };
    let called = 0, shrunkTotal = 0, rawTotal = 0;

    CELLS.forEach((cell, r) => {
      const rows = DEMO_GRID.filter(g => g.cell === cell);
      const mean = precisionMean(rows);
      rows.forEach((g, c) => {
        const p = shrink(g.est, g.se, mean, v.tau);
        const excl = Math.abs(p.est) > 1.959964 * p.sd;
        if (excl) called++;
        shrunkTotal += Math.abs(p.est); rawTotal += Math.abs(g.est);
        const x = L + c * cw, y = T + r * ch;
        svg.append(s('rect', { x, y, width: cw - 2, height: ch - 2, fill: colour(p.est),
          stroke: excl ? '#ffffff' : 'rgba(0,0,0,0.06)', 'stroke-width': excl ? 2 : 1, rx: 2 }));
        const t = s('text', { x: x + (cw - 2) / 2, y: y + ch / 2 + 3, 'text-anchor': 'middle',
          'font-size': 10, fill: Math.abs(p.est) / v.lim > 0.55 ? '#fff' : INK });
        t.textContent = p.est.toFixed(2).replace('-', '−'); svg.append(t);
      });
      const rl = s('text', { x: L - 8, y: T + r * ch + ch / 2 + 3, 'text-anchor': 'end', fill: INK, 'font-size': 11 });
      rl.textContent = cell; svg.append(rl);
    });
    SIGS.forEach((sig, c) => {
      const t = s('text', { x: L + c * cw + (cw - 2) / 2, y: T - 10, 'text-anchor': 'middle', fill: GREY, 'font-size': 11 });
      t.textContent = sig; svg.append(t);
    });
    // scale bar
    const bx = L, by = T + CELLS.length * ch + 16, bw = SIGS.length * cw - 2;
    for (let i = 0; i < 40; i++) {
      const b = -v.lim + (i / 39) * 2 * v.lim;
      svg.append(s('rect', { x: bx + (i / 40) * bw, y: by, width: bw / 40 + 0.6, height: 9, fill: colour(b) }));
    }
    const sl = (x, t, a) => { const n = s('text', { x, y: by + 24, 'text-anchor': a, fill: GREY, 'font-size': 11 }); n.textContent = t; return n; };
    svg.append(sl(bx, '−' + v.lim.toFixed(2), 'start'), sl(bx + bw / 2, 'BMI slope', 'middle'), sl(bx + bw, v.lim.toFixed(2), 'end'));

    const pull = rawTotal ? (1 - shrunkTotal / rawTotal) * 100 : 0;
    return {
      eq: `τ = ${v.tau.toFixed(2)} · mean magnitude pulled in by <b>${pull.toFixed(0)}%</b>`,
      cards: [
        ['τ', v.tau.toFixed(2), 'between-signature SD'],
        ['HDI excludes zero', String(called) + ' / ' + DEMO_GRID.length, 'white outline', true],
        ['shrinkage', pull.toFixed(0) + '%', 'reduction in magnitude'],
        ['cells', String(DEMO_GRID.length), CELLS.length + ' cell types']
      ],
      verdict: v.tau < 0.12
        ? ['Heavy pooling', 'At this τ the model treats signatures within a cell type as nearly interchangeable, and individual estimates almost vanish into the cell-type mean.']
        : v.tau > 0.8
          ? ['Almost no pooling', 'Each signature is estimated on its own evidence. This is close to fitting every feature separately, which is what the paired and PTM models do.']
          : ['Balanced pooling', 'Well-supported signatures keep their estimates while weakly supported ones are pulled towards the cell-type mean.']
    };
  }
};

/* ---------- shared: cohort and least squares, mirroring the simulation notebook ---------- */
// patient_frame(): three clipped Gaussian BMI groups, age confounded with BMI, logistic sex.
function cohort(counts, confounding, seed) {
  const rng = lcg(seed);
  const spec = [[22.2, 1.4, 18.5, 24.99], [27.4, 1.3, 25, 29.99], [34, 3, 30, 45]];
  const bmi = [], group = [];
  counts.forEach((n, g) => {
    const [m, sd, lo, hi] = spec[g];
    for (let i = 0; i < n; i++) { bmi.push(clamp(m + sd * rng.n(), lo, hi)); group.push(g); }
  });
  const mean = bmi.reduce((a, b) => a + b, 0) / bmi.length;
  const v = bmi.reduce((a, b) => a + (b - mean) ** 2, 0) / (bmi.length - 1);
  const sd = Math.sqrt(v);
  const z = bmi.map(b => (b - mean) / sd);
  const age = z.map(zi => confounding * zi + Math.sqrt(1 - confounding * confounding) * rng.n());
  const sex = z.map(zi => (rng.u() < 1 / (1 + Math.exp(0.15 - 0.3 * zi)) ? 1 : 0) - 0.5);
  return { bmi, z, group, age, sex, n: bmi.length, sd };
}
// seeded generator with normal, uniform and gamma draws
function lcg(seed) {
  let x = seed >>> 0;
  const u = () => { x = (Math.imul(x, 1664525) + 1013904223) >>> 0; return (x + 0.5) / 4294967296; };
  const n = () => Math.sqrt(-2 * Math.log(u())) * Math.cos(2 * Math.PI * u());
  // Marsaglia-Tsang gamma, for chi-square draws behind the Student-t scale mixture
  const gam = k => {
    if (k < 1) return gam(k + 1) * Math.pow(u(), 1 / k);
    const d = k - 1 / 3, c = 1 / Math.sqrt(9 * d);
    for (;;) {
      let xx, vv;
      do { xx = n(); vv = 1 + c * xx; } while (vv <= 0);
      vv = vv * vv * vv;
      const uu = u();
      if (uu < 1 - 0.0331 * xx ** 4) return d * vv;
      if (Math.log(uu) < 0.5 * xx * xx + d * (1 - vv + Math.log(vv))) return d * vv;
    }
  };
  // Student-t as a normal scale mixture: exact t marginals for any nu > 0
  const t = nu => n() / Math.sqrt(gam(nu / 2) * 2 / nu);
  return { u, n, gam, t };
}
// least squares by normal equations; X is an array of rows
function lstsq(X, y) {
  const p = X[0].length, A = [], b = [];
  for (let i = 0; i < p; i++) {
    A.push(new Array(p).fill(0)); b.push(0);
    for (let r = 0; r < X.length; r++) {
      b[i] += X[r][i] * y[r];
      for (let j = 0; j < p; j++) A[i][j] += X[r][i] * X[r][j];
    }
  }
  for (let i = 0; i < p; i++) {
    let piv = i;
    for (let r = i + 1; r < p; r++) if (Math.abs(A[r][i]) > Math.abs(A[piv][i])) piv = r;
    [A[i], A[piv]] = [A[piv], A[i]]; [b[i], b[piv]] = [b[piv], b[i]];
    const d = A[i][i] || 1e-12;
    for (let j = i; j < p; j++) A[i][j] /= d;
    b[i] /= d;
    for (let r = 0; r < p; r++) if (r !== i && A[r][i]) {
      const f = A[r][i];
      for (let j = i; j < p; j++) A[r][j] -= f * A[i][j];
      b[r] -= f * b[i];
    }
  }
  return b;
}
const axes = (svg, L, R, T, B, xlab, ylab) => {
  svg.append(s('line', { x1: L, x2: R, y1: B, y2: B, stroke: INK }));
  svg.append(s('line', { x1: L, x2: L, y1: T, y2: B, stroke: INK }));
  const a = s('text', { x: (L + R) / 2, y: B + 30, 'text-anchor': 'middle', fill: GREY, 'font-size': 11 });
  a.textContent = xlab; svg.append(a);
  const c = s('text', { x: 13, y: (T + B) / 2, 'text-anchor': 'middle', fill: GREY, 'font-size': 11,
    transform: `rotate(-90 13 ${(T + B) / 2})` });
  c.textContent = ylab; svg.append(c);
};

/* ---------- widget: the cellular model, Equations 2-5 ---------- */
const cellular = {
  id: 'cellular',
  title: 'The cellular RNA model (Student-t intercept)',
  eq: '<span data-term="Y">Y<sub>ij</sub></span> ~ <span data-term="t">t<sub>ν</sub></span>(<span data-term="mu">μ<sub>ij</sub></span>, <span data-term="sigma">σ<sub>j</sub></span>), <span data-term="mu">μ<sub>ij</sub></span> = <span data-term="alpha">α<sub>j</sub></span> + <span data-term="beta">β<sub>j</sub>BMI<sub>i</sub></span> + <span data-term="gamma">γ<sub>j</sub>age<sub>i</sub> + δ<sub>j</sub>sex<sub>i</sub></span> + <span data-term="v">v<sub>i,c(j)</sub></span>, <span data-term="v">v<sub>ic</sub> ~ t<sub>ν<sub>v</sub></sub>(0, ω<sub>c</sub>)</span>',
  terms: {
    Y: ['Signature score', 'Score of signature j in patient i, standardized across patients. Each signature belongs to one cell type c(j).'],
    t: ['Student-t likelihood', 'Scores scatter around their fitted value with heavier tails than a normal curve, so a few extreme patients pull the estimates less.'],
    mu: ['Fitted score', 'What the model expects patient i to score on signature j, given BMI, age, sex and the patient × cell-type intercept.'],
    sigma: ['Residual scale', 'How far scores typically fall from the fitted value. Each signature has its own scale, pooled within its cell type.'],
    alpha: ['Feature intercept', 'The average level of signature j. It does not change with BMI.'],
    beta: ['BMI slope', 'The quantity the paper reports: the change in the score per unit of standardized BMI, shown per 5 kg/m². Slopes in one cell type are pulled towards their cell-type mean (Equation 4).'],
    gamma: ['Age and sex', 'Adjustment terms. Age is related to BMI, so leaving it out would mix age effects into the BMI slope.'],
    v: ['Patient × cell-type intercept', 'One shift shared by every signature of cell type c in patient i. It enters once per patient and cell type, so a shift seen in many signatures is not counted as many separate pieces of evidence. Its scale ω<sub>c</sub> has a HalfNormal(0.50) prior and its tails a Student-t with ν<sub>v</sub> − 2 ~ Exponential(0.10).']
  },
  ref: 'Equations 2–5',
  blurb: 'One signature score, one patient at a time. Every term adds to the fitted mean, and only β is the BMI association. The term v is a patient × cell-type intercept: one shift shared by all signatures of that cell type in that patient. Move the highlighted patient along the cohort to see how the same slope gives very different scores.',
  note: 'BMI and age are standardized; sex is coded −0.5 female, +0.5 male. Reported slopes are posterior draws multiplied by 5/sᴮᴹᴵ for a 5 kg/m² increase, with the cohort SD 4.462 kg/m². BMI slopes are pooled within cell types around a cell-type mean m with spread τ, where m ~ N(0, 0.20²) and τ ~ HalfNormal(0.30). The intercept v is drawn here from a Student-t with 4 degrees of freedom. Deconvolution and score-construction uncertainty are not propagated.',
  params: [
    { k: 'beta',  label: 'β<sub>j</sub> BMI slope', min: -0.8, max: 0.8, step: 0.01, value: 0.30, colour: MAGENTA, hint: 'The quantity the paper reports.' },
    { k: 'gamma', label: 'γ<sub>j</sub> age coefficient', min: -0.5, max: 0.5, step: 0.01, value: 0.12, colour: BLUE, hint: 'Adjusted for, because age correlates with BMI.' },
    { k: 'delta', label: 'δ<sub>j</sub> sex coefficient', min: -0.5, max: 0.5, step: 0.01, value: -0.08, colour: BLUE, hint: 'Sex is centred, so this is the male minus female difference.' },
    { k: 'su',    label: 'ω<sub>c</sub> intercept scale', min: 0, max: 1, step: 0.01, value: 0.35, colour: TEAL, hint: 'Scale of the patient × cell-type intercept. Prior HalfNormal(0.5).' },
    { k: 'sigma', label: 'σ<sub>j</sub> residual scale', min: 0.1, max: 1.4, step: 0.01, value: 0.65, colour: GREY, hint: 'Student-t scale. The simulation centres this at 0.65.' },
    { k: 'who',   label: 'highlighted patient', min: 0, max: 126, step: 1, value: 96, colour: TEAL, fmt: x => '#' + (x + 1), hint: 'Pick a patient to decompose their fitted score.' }
  ],
  presets: [
    { label: 'Reference', v: { beta: 0.30, gamma: 0.12, delta: -0.08, su: 0.35, sigma: 0.65 } },
    { label: 'Null', v: { beta: 0, gamma: 0.12, delta: -0.08, su: 0.35, sigma: 0.65 } },
    { label: 'Noise dominates', v: { beta: 0.30, gamma: 0.12, delta: -0.08, su: 0.35, sigma: 1.3 } },
    { label: 'Strong intercept', v: { beta: 0.30, gamma: 0.12, delta: -0.08, su: 0.9, sigma: 0.45 } }
  ],
  init() {
    const c = cohort([52, 57, 18], 0.2, 20260926);
    const r = lcg(515151);
    // Student-t draws with 4 degrees of freedom for the patient x cell-type intercept.
    const t4 = () => r.n() / Math.sqrt((r.n() ** 2 + r.n() ** 2 + r.n() ** 2 + r.n() ** 2) / 4);
    return { c, u: Array.from({ length: c.n }, t4), e: Array.from({ length: c.n }, () => r.n()) };
  },
  draw(svg, v, mem) {
    const { c, u, e } = mem, who = Math.min(c.n - 1, Math.round(v.who));
    const alpha = 0;
    const mu = i => alpha + v.beta * c.z[i] + v.gamma * c.age[i] + v.delta * c.sex[i] + v.su * u[i];
    const y = i => mu(i) + v.sigma * e[i];
    svg.setAttribute('viewBox', '0 0 360 268');
    svg.replaceChildren();
    const L = 46, R = 348, T = 18, B = 214;
    const X = z => L + ((clamp(z, -2.6, 3.4) + 2.6) / 6) * (R - L);
    const lim = 2.6, Y = q => clamp(116 - (q / lim) * 96, T, B);
    for (const g of [-2, -1, 0, 1, 2]) svg.append(s('line', { x1: L, x2: R, y1: Y(g), y2: Y(g), stroke: LINE, 'stroke-dasharray': g ? '3 3' : '' }));
    const cols = [TEAL, BLUE, MAGENTA];
    for (let i = 0; i < c.n; i++) {
      if (i === who) continue;
      svg.append(s('circle', { cx: X(c.z[i]), cy: Y(y(i)), r: 2.5, fill: cols[c.group[i]], opacity: 0.4 }));
    }
    svg.append(s('line', { x1: X(-2.6), y1: Y(v.beta * -2.6), x2: X(3.4), y2: Y(v.beta * 3.4), stroke: INK, 'stroke-width': 2.2 }));
    const bl = s('text', { x: R - 4, y: Y(v.beta * 3.4) - 8, 'text-anchor': 'end', fill: INK, 'font-size': 11, 'font-weight': 700 });
    bl.textContent = 'α + β·BMI'; svg.append(bl);
    // highlighted patient: dashed drop from the line to the realised score
    svg.append(s('line', { x1: X(c.z[who]), x2: X(c.z[who]), y1: Y(v.beta * c.z[who]), y2: Y(y(who)),
      stroke: MAGENTA, 'stroke-width': 1.4, 'stroke-dasharray': '4 3' }));
    svg.append(s('circle', { cx: X(c.z[who]), cy: Y(y(who)), r: 5.5, fill: MAGENTA, stroke: '#fff', 'stroke-width': 1.6 }));
    axes(svg, L, R, T, B, 'standardized BMI', 'signature score');
    for (const g of [-2, 0, 2]) {
      const t = s('text', { x: X(g), y: B + 15, 'text-anchor': 'middle', fill: GREY, 'font-size': 11 });
      t.textContent = String(g); svg.append(t);
    }
    const key = s('text', { x: L, y: 262, fill: GREY, 'font-size': 10.5 });
    key.textContent = 'normal weight 52 · overweight 57 · obese 18'; svg.append(key);

    const parts = [
      ['β·BMI', v.beta * c.z[who], MAGENTA],
      ['γ·age', v.gamma * c.age[who], BLUE],
      ['δ·sex', v.delta * c.sex[who], BLUE],
      ['v (patient × cell type)', v.su * u[who], TEAL],
      ['residual', v.sigma * e[who], GREY]
    ];
    const biggest = parts.reduce((a, p) => Math.abs(p[1]) > Math.abs(a[1]) ? p : a);
    return {
      eq: `patient #${who + 1}: BMI ${c.bmi[who].toFixed(1)} · score <b>${minus(y(who))}</b> = ` +
          parts.map(p => minus(p[1])).join(' + ').replace(/\+ −/g, '− '),
      cards: [
        ['β·BMI', minus(parts[0][1]), 'the reported effect', true],
        ['γ·age + δ·sex', minus(parts[1][1] + parts[2][1]), 'adjusted away'],
        ['v<sub>i,c</sub> intercept', minus(parts[3][1]), 'shared within the cell type'],
        ['residual', minus(parts[4][1]), 'Student-t noise']
      ],
      verdict: Math.abs(parts[0][1]) >= Math.abs(biggest[1])
        ? ['BMI dominates this patient', 'The BMI term is the largest contribution to their fitted score. With 127 patients the slope is still estimated across the whole cohort, not from this one.']
        : ['BMI is not the largest term here', `For this patient the ${biggest[0]} term contributes more than BMI does. A real association can sit well below the per-patient noise and still be estimable across the cohort.`]
    };
  }
};


/* ---------- widget: the cellular model, orthogonal patient-effect version ---------- */
const orthCellular = {
  ...cellular,
  id: 'orth-cellular',
  title: 'The cellular RNA model (orthogonal patient effect)',
  eq: '<span data-term="Y">Y<sub>ij</sub></span> ~ <span data-term="t">t<sub>ν</sub></span>(<span data-term="mu">μ<sub>ij</sub></span>, <span data-term="sigma">σ<sub>j</sub></span>), <span data-term="mu">μ<sub>ij</sub></span> = <span data-term="alpha">α<sub>j</sub></span> + <span data-term="beta">β<sub>j</sub>BMI<sub>i</sub></span> + <span data-term="gamma">γ<sub>j</sub>age<sub>i</sub> + δ<sub>j</sub>sex<sub>i</sub></span> + <span data-term="u">u<sub>i</sub></span>, <span data-term="orth">Xᵀu = 0</span>',
  terms: {
    Y: cellular.terms.Y,
    t: ['Student-t likelihood', 'Scores scatter around their fitted value with heavier tails than a normal curve. This part is the same in both models.'],
    mu: ['Fitted score', 'What the model expects patient i to score on signature j, given BMI, age, sex and the single patient effect u<sub>i</sub>.'],
    sigma: cellular.terms.sigma,
    alpha: cellular.terms.alpha,
    beta: cellular.terms.beta,
    gamma: cellular.terms.gamma,
    u: ['One patient effect for all features', 'Each patient gets a single shift u<sub>i</sub>, shared by every signature in every cell type, with scale σ<sub>u</sub>. The Student-t intercept model gives each patient one shift per cell type instead.'],
    orth: ['Orthogonality constraint', 'The vector u is forced to be uncorrelated with BMI, age and sex (the design X), by building it on the SVD null space of X. It can only absorb variation the covariates do not explain.']
  },
  ref: 'Comparison model, Equations 2–3 of the earlier version',
  blurb: 'The same signature score, written for the orthogonal patient-effect model. There is one patient shift u for every feature, held orthogonal to the covariates, so it cannot soak up anything BMI explains. This is the model kept for comparison: its intervals for cell-type means were too narrow in simulation. Move the highlighted patient to see their score split into the same parts.',
  note: 'The continuous design has rank four, leaving 123 basis coefficients for 127 patients; u = σᵤBη with the scaling set so the average prior variance is σᵤ². Likelihood, covariates and BMI pooling are the same as in the Student-t intercept model. Orthogonality is exact in the full cohort and approximate among the patients kept for one cell type. The patients here are generated; patient-level data is controlled access.',
  params: [
    cellular.params[0], cellular.params[1], cellular.params[2],
    { k: 'su', label: 'σ<sub>u</sub> patient effect SD', min: 0, max: 1, step: 0.01, value: 0.5, colour: TEAL, hint: 'Scale of the single orthogonal patient effect. Prior HalfNormal(0.5).' },
    cellular.params[4], cellular.params[5]
  ],
  presets: [
    { label: 'Reference', v: { beta: 0.30, gamma: 0.12, delta: -0.08, su: 0.5, sigma: 0.65 } },
    { label: 'Null', v: { beta: 0, gamma: 0.12, delta: -0.08, su: 0.5, sigma: 0.65 } },
    { label: 'Noise dominates', v: { beta: 0.30, gamma: 0.12, delta: -0.08, su: 0.5, sigma: 1.3 } },
    { label: 'Large patient effect', v: { beta: 0.30, gamma: 0.12, delta: -0.08, su: 1.0, sigma: 0.45 } }
  ],
  init() {
    const c = cohort([52, 57, 18], 0.2, 20260926);
    const r = lcg(424242);
    let raw = Array.from({ length: c.n }, () => r.n());
    const X = c.z.map((z, i) => [1, z, c.age[i], c.sex[i]]);
    const b = lstsq(X, raw);
    raw = raw.map((q, i) => q - X[i].reduce((a, xv, j) => a + xv * b[j], 0));
    const sd = Math.sqrt(raw.reduce((a, q) => a + q * q, 0) / raw.length);
    return { c, u: raw.map(q => q / sd), e: Array.from({ length: c.n }, () => r.n()) };
  },
  draw(svg, v, mem) {
    const out = cellular.draw(svg, v, mem);
    out.cards[2] = ['u<sub>i</sub> patient effect', out.cards[2][1], 'one shift for every feature'];
    return out;
  }
};

/* ---------- widget: the patient x cell-type intercept and interval coverage ---------- */
// Normal CDF via the Abramowitz-Stegun erf approximation (error < 1.5e-7).
function Phi(x) {
  const t = 1 / (1 + 0.3275911 * Math.abs(x) / Math.SQRT2);
  const y = 1 - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t * Math.exp(-x * x / 2);
  return x >= 0 ? (1 + y) / 2 : (1 - y) / 2;
}
const intercept = {
  id: 'intercept',
  title: 'Why the patient × cell-type intercept matters',
  eq: '<span data-term="v">v<sub>ic</sub> ~ t<sub>ν<sub>v</sub></sub>(0, ω<sub>c</sub>)</span> shared by all <span data-term="K">K signatures</span> of cell type c in patient i, <span data-term="keff">K<sub>eff</sub> = (ω² + σ²) / (ω² + σ²/K)</span>',
  terms: {
    v: ['Patient × cell-type intercept', 'A shift common to all signatures of one cell type in one patient, for example because the deconvolved profile of that cell type is brighter overall in that tumour.'],
    K: ['Signatures per cell type', 'The number of signatures scored in the same deconvolved profile. Across the three analyses there are on average 28 to 38 signatures per cell type (914 features in 26 cell types).'],
    keff: ['Effective number of signatures', 'How many independent signatures the K correlated ones are worth. With a large shared shift ω, K signatures carry little more information than one.']
  },
  ref: 'Equations 3 and 5, simulation',
  blurb: 'Signatures scored in the same deconvolved profile move together within a patient. A model with no term for that shared shift treats them as independent, so its intervals for the cell-type mean are too narrow. Change the shift size and the number of signatures to see how far the coverage of a nominal 95% interval falls.',
  note: 'Exact arithmetic for the cell-type mean of K signatures, each with shared patient × cell-type variance ω² and independent residual variance σ². The correct variance per patient is ω² + σ²/K; a model that ignores the shared part uses (ω² + σ²)/K. The number of patients cancels. This isolates one mechanism and does not reproduce the full simulation, where the orthogonal patient-effect model covered 0.50 to 0.66.',
  params: [
    { k: 'om', label: 'ω<sub>c</sub> shared shift SD', min: 0, max: 1, step: 0.01, value: 0.35, colour: TEAL, hint: 'Size of the shift shared by all signatures of a cell type in a patient.' },
    { k: 'sg', label: 'σ residual SD', min: 0.2, max: 1.5, step: 0.01, value: 0.65, colour: GREY, hint: 'Independent noise of each signature.' },
    { k: 'K', label: 'K signatures in the cell type', min: 1, max: 60, step: 1, value: 35, colour: MAGENTA, fmt: x => String(Math.round(x)), hint: 'More signatures add little once the shared shift dominates.' }
  ],
  presets: [
    { label: 'Reference', v: { om: 0.35, sg: 0.65, K: 35 } },
    { label: 'No shared shift', v: { om: 0, sg: 0.65, K: 35 } },
    { label: 'Large shared shift', v: { om: 0.8, sg: 0.65, K: 35 } },
    { label: 'Few signatures', v: { om: 0.35, sg: 0.65, K: 4 } }
  ],
  init() { return {}; },
  draw(svg, v) {
    const cov = K => { const r = Math.sqrt(((v.om ** 2 + v.sg ** 2) / K) / (v.om ** 2 + v.sg ** 2 / K)); return 2 * Phi(1.959964 * r) - 1; };
    const K = Math.round(v.K), keff = (v.om ** 2 + v.sg ** 2) / (v.om ** 2 + v.sg ** 2 / K), c = cov(K);
    svg.setAttribute('viewBox', '0 0 360 268');
    svg.replaceChildren();
    const L = 46, R = 348, T = 18, B = 214;
    const X = k => L + (k - 1) / 59 * (R - L), Y = q => B - q * (B - T);
    for (const g of [0.25, 0.5, 0.75, 0.95]) {
      svg.append(s('line', { x1: L, x2: R, y1: Y(g), y2: Y(g), stroke: LINE, 'stroke-dasharray': g === 0.95 ? '' : '3 3' }));
      const t = s('text', { x: L - 6, y: Y(g) + 4, 'text-anchor': 'end', fill: GREY, 'font-size': 10 }); t.textContent = g.toFixed(2); svg.append(t);
    }
    let path = '';
    for (let k = 1; k <= 60; k++) path += (k === 1 ? 'M' : 'L') + X(k).toFixed(1) + ' ' + Y(cov(k)).toFixed(1);
    svg.append(s('path', { d: path, fill: 'none', stroke: MAGENTA, 'stroke-width': 2.4 }));
    svg.append(s('line', { x1: L, x2: R, y1: Y(0.95), y2: Y(0.95), stroke: TEAL, 'stroke-width': 2.4 }));
    const a = s('text', { x: R - 4, y: Y(0.95) - 7, 'text-anchor': 'end', fill: TEAL, 'font-size': 11, 'font-weight': 700 }); a.textContent = 'with the intercept: 0.95'; svg.append(a);
    const b = s('text', { x: R - 4, y: Y(cov(60)) + 15, 'text-anchor': 'end', fill: MAGENTA, 'font-size': 11, 'font-weight': 700 }); b.textContent = 'shared shift ignored'; svg.append(b);
    svg.append(s('circle', { cx: X(K), cy: Y(c), r: 5.5, fill: MAGENTA, stroke: '#fff', 'stroke-width': 1.6 }));
    axes(svg, L, R, T, B, 'signatures in the cell type (K)', 'coverage of a 95% interval');
    for (const k of [1, 20, 40, 60]) { const t = s('text', { x: X(k), y: B + 15, 'text-anchor': 'middle', fill: GREY, 'font-size': 11 }); t.textContent = String(k); svg.append(t); }
    const ratio = Math.sqrt(((v.om ** 2 + v.sg ** 2) / K) / (v.om ** 2 + v.sg ** 2 / K));
    return {
      eq: `K = ${K} signatures are worth <b>${keff.toFixed(1)}</b> independent ones · intervals that ignore the shared shift are <b>${(ratio * 100).toFixed(0)}%</b> of the correct width`,
      cards: [
        ['coverage, shift ignored', c.toFixed(2), 'nominal 0.95', true],
        ['coverage, with intercept', '0.95', 'by construction'],
        ['K<sub>eff</sub>', keff.toFixed(1), 'of ' + K + ' signatures'],
        ['interval width', (ratio * 100).toFixed(0) + '%', 'of the correct width']
      ],
      verdict: c < 0.9
        ? ['The intervals are too narrow', `Ignoring a shared shift of SD ${v.om.toFixed(2)} across ${K} signatures leaves the nominal 95% interval covering only ${(c * 100).toFixed(0)}% of the time. The patient × cell-type intercept puts that shared variation in one term, so the cell-type mean gets an interval of the right width.`]
        : ['Little is lost here', 'With a small shared shift or few signatures, treating signatures as independent costs little coverage. Raise ω or K to see the intervals shrink below their nominal level.']
    };
  }
};

/* ---------- widget: the orthogonal patient-effect model (comparison) ---------- */
const patient = {
  id: 'patient-effect',
  title: 'The orthogonal patient-effect model (comparison)',
  eq: '<span data-term="orth">Xᵀu = 0</span>, <span data-term="u">u = σ<sub>u</sub>·Bη</span>, <span data-term="B">B from the SVD null space of X</span>',
  terms: {
    orth: ['Orthogonality constraint', 'The patient effects u are forced to be uncorrelated with BMI, age and sex (the design X). They can only absorb variation the covariates do not explain.'],
    u: ['One patient effect for all features', 'Each patient gets a single shift u<sub>i</sub>, shared by every signature in every cell type, with scale σ<sub>u</sub>. The Student-t intercept model instead gives each patient one shift per cell type.'],
    B: ['Basis', 'B spans the directions left over after removing BMI, age and sex. 127 patients minus 4 design columns leaves 123 basis coefficients η.']
  },
  ref: 'Comparison model, Supplementary Methods',
  blurb: 'In the comparison model, Equation 3 carries one patient effect u shared by all features, built on a basis orthogonal to BMI, age and sex. It cannot absorb anything the covariates explain. The cost is shown here: a BMI-linked shift shared by every signature has nowhere to go but the slopes. In simulation this model covered the true effects only 50 to 66% of the time.',
  note: 'The continuous design has rank four, giving 123 basis coefficients for 127 patients. Each patient effect enters every feature with loading one, and the scaling makes the average prior variance equal σᵤ². The basis was fixed before low-read exclusions, so orthogonality holds exactly in the full cohort and only approximately among the patients retained for a given cell type.',
  params: [
    { k: 'su',    label: 'σ<sub>u</sub> patient effect SD', min: 0, max: 1.2, step: 0.01, value: 0.5, colour: TEAL, hint: 'Prior HalfNormal(0.5).' },
    { k: 'shift', label: 'BMI-linked shift shared by all signatures', min: 0, max: 0.6, step: 0.01, value: 0.2, colour: MAGENTA, hint: 'A technical or compositional effect that moves every signature together with BMI.' }
  ],
  presets: [
    { label: 'No shared shift', v: { su: 0.5, shift: 0 } },
    { label: 'Moderate shift', v: { su: 0.5, shift: 0.2 } },
    { label: 'Large shift', v: { su: 0.5, shift: 0.5 } },
    { label: 'Unrestricted patient effect', v: { su: 1.0, shift: 0.2 } }
  ],
  init() {
    const c = cohort([52, 57, 18], 0.2, 20260926);
    const r = lcg(31337);
    // orthogonalise raw patient draws against [1, bmi, age, sex] the way the SVD basis does
    let raw = Array.from({ length: c.n }, () => r.n());
    const X = c.z.map((z, i) => [1, z, c.age[i], c.sex[i]]);
    const b = lstsq(X, raw);
    raw = raw.map((q, i) => q - X[i].reduce((a, xv, j) => a + xv * b[j], 0));
    const sd = Math.sqrt(raw.reduce((a, q) => a + q * q, 0) / raw.length);
    return { c, u: raw.map(q => q / sd) };
  },
  draw(svg, v, mem) {
    const { c, u } = mem;
    const pts = c.z.map((z, i) => ({ z, u: v.su * u[i] }));
    const X = pts.map(p => [1, p.z]);
    const fit = lstsq(X, pts.map(p => p.u));
    svg.setAttribute('viewBox', '0 0 360 268');
    svg.replaceChildren();
    const L = 48, R = 348, T = 18, B = 214;
    const PX = z => L + ((clamp(z, -2.6, 3.4) + 2.6) / 6) * (R - L);
    const lim = 1.6, PY = q => clamp(116 - (q / lim) * 96, T, B);
    for (const g of [-1, 0, 1]) svg.append(s('line', { x1: L, x2: R, y1: PY(g), y2: PY(g), stroke: LINE, 'stroke-dasharray': g ? '3 3' : '' }));
    for (const p of pts) svg.append(s('circle', { cx: PX(p.z), cy: PY(p.u), r: 2.6, fill: TEAL, opacity: 0.45 }));
    svg.append(s('line', { x1: PX(-2.6), y1: PY(fit[0] + fit[1] * -2.6), x2: PX(3.4), y2: PY(fit[0] + fit[1] * 3.4),
      stroke: TEAL, 'stroke-width': 2.4 }));
    const fl = s('text', { x: L + 6, y: PY(0) - 8, fill: TEAL, 'font-size': 11, 'font-weight': 700 });
    fl.textContent = 'u against BMI: slope ' + fit[1].toExponential(1); svg.append(fl);
    // the shift u cannot represent
    svg.append(s('line', { x1: PX(-2.6), y1: PY(v.shift * -2.6), x2: PX(3.4), y2: PY(v.shift * 3.4),
      stroke: MAGENTA, 'stroke-width': 2.4, 'stroke-dasharray': '7 4' }));
    const sl = s('text', { x: PX(3.4) - 4, y: PY(v.shift * 3.4) - 8, 'text-anchor': 'end', fill: MAGENTA, 'font-size': 11, 'font-weight': 700 });
    sl.textContent = 'shared BMI-linked shift'; svg.append(sl);
    axes(svg, L, R, T, B, 'standardized BMI', 'patient effect u');
    for (const g of [-2, 0, 2]) {
      const t = s('text', { x: PX(g), y: B + 15, 'text-anchor': 'middle', fill: GREY, 'font-size': 11 });
      t.textContent = String(g); svg.append(t);
    }
    return {
      eq: `u against BMI: slope ${fit[1].toExponential(2)} (zero by construction) · shift carried by every feature slope instead: <b>+${minus(v.shift)}</b>`,
      cards: [
        ['slope of u on BMI', fit[1].toExponential(1), 'orthogonality holds', true],
        ['σ<sub>u</sub>', minus(v.su), 'spread of u'],
        ['shared shift', minus(v.shift), 'per BMI SD'],
        ['added to every β<sub>j</sub>', '+' + minus(v.shift), 'cannot go into u']
      ],
      verdict: v.shift > 0.02
        ? ['This is the limitation the paper names', 'The dashed line is a BMI-linked change common to every signature — compositional drift, a batch effect, anything global. Because u is orthogonal to BMI by construction, it cannot absorb that line. The whole of it is carried by the feature slopes and cell-type means instead, inflating every one of them by the same amount.']
        : ['Orthogonality doing its job', 'With no shared BMI-linked shift, u captures patient-to-patient variation while leaving the BMI slopes untouched. Raise the shift slider to see the cost of that design choice.']
    };
  }
};


/* ---------- widget: side-by-side model comparison ---------- */
const NICE = [0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.2, 1.6];
function cmpStats(v) {
  const K = Math.round(v.K), n = Math.round(v.n);
  const seT = Math.sqrt((v.om ** 2 + v.sg ** 2 / K) / n);
  const seO = Math.sqrt(((v.om ** 2 + v.sg ** 2) / K) / n);
  return { K, n, seT, seO, ratio: seO / seT, keff: (v.om ** 2 + v.sg ** 2) / (v.om ** 2 + v.sg ** 2 / K) };
}
function cmpPanel(svg, v, mem, model) {
  const st = cmpStats(v), isT = model === 'T';
  const se = isT ? st.seT : st.seO, hw = 1.959964 * se;
  const col = isT ? TEAL : MAGENTA;
  svg.replaceChildren();
  if (v.view < 0.5) {
    svg.setAttribute('viewBox', '0 0 360 286');
    const N = 40, L = 52, R = 348, T = 16, B = 244, row = (B - T) / N;
    const lim = NICE.find(q => q >= Math.abs(v.beta) + 3.3 * st.seT) || 1.6;
    const X = b => L + ((clamp(b, -lim, lim) + lim) / (2 * lim)) * (R - L);
    const z = mem.z[Math.round(v.seed) - 1];
    let cover = 0, called = 0, wsum = 0;
    for (const g of [-lim, -lim / 2, lim / 2, lim]) svg.append(s('line', { x1: X(g), x2: X(g), y1: T, y2: B, stroke: LINE, 'stroke-dasharray': '3 3' }));
    svg.append(s('line', { x1: X(0), x2: X(0), y1: T, y2: B, stroke: GREY, 'stroke-width': 1.2 }));
    for (let r = 0; r < N; r++) {
      const est = v.beta + st.seT * z[r], y = T + (r + 0.5) * row;
      const hit = Math.abs(est - v.beta) <= hw, excl = Math.abs(est) > hw;
      if (hit) cover++; if (excl) called++;
      wsum += hw;
      const c = hit ? col : INK;
      svg.append(s('line', { x1: X(est - hw), x2: X(est + hw), y1: y, y2: y, stroke: hit ? col : '#D81B60', 'stroke-width': 1.6, opacity: hit ? 0.55 : 0.95 }));
      svg.append(s('circle', { cx: X(est), cy: y, r: 1.9, fill: hit ? col : '#D81B60' }));
    }
    svg.append(s('line', { x1: X(v.beta), x2: X(v.beta), y1: T - 4, y2: B, stroke: INK, 'stroke-width': 1.8, 'stroke-dasharray': '6 3' }));
    const tl = s('text', { x: X(v.beta), y: T - 6, 'text-anchor': 'middle', fill: INK, 'font-size': 10.5, 'font-weight': 700 }); tl.textContent = 'true ' + minus(v.beta); svg.append(tl);
    axes(svg, L, R, T, B, 'estimated slope per BMI SD', '40 repeated cohorts');
    for (const g of [-lim, 0, lim]) { const q = s('text', { x: X(g), y: B + 14, 'text-anchor': 'middle', fill: GREY, 'font-size': 10.5 }); q.textContent = minus(+g.toFixed(2)); svg.append(q); }
    return { cover, called, hw, N };
  }
  svg.setAttribute('viewBox', '0 0 360 286');
  const L = 52, R = 348, T = 16, B = 244;
  const wAt = (om, which) => 1.959964 * Math.sqrt(which === 'T' ? (om ** 2 + v.sg ** 2 / st.K) / st.n : ((om ** 2 + v.sg ** 2) / st.K) / st.n);
  const ymax = wAt(1, 'T') * 1.08, X = om => L + om * (R - L), Y = w => B - (w / ymax) * (B - T);
  for (const g of [0.25, 0.5, 0.75, 1]) svg.append(s('line', { x1: L, x2: R, y1: Y(wAt(g, 'T')), y2: Y(wAt(g, 'T')), stroke: LINE, 'stroke-dasharray': '3 3', opacity: 0 }));
  const curve = which => { let d = ''; for (let i = 0; i <= 40; i++) { const om = i / 40; d += (i ? 'L' : 'M') + X(om).toFixed(1) + ' ' + Y(wAt(om, which)).toFixed(1); } return d; };
  if (!isT) svg.append(s('path', { d: curve('T'), fill: 'none', stroke: TEAL, 'stroke-width': 1.6, 'stroke-dasharray': '5 4', opacity: 0.7 }));
  svg.append(s('path', { d: curve(model), fill: 'none', stroke: col, 'stroke-width': 2.6 }));
  svg.append(s('circle', { cx: X(v.om), cy: Y(hw), r: 5.5, fill: col, stroke: '#fff', 'stroke-width': 1.6 }));
  axes(svg, L, R, T, B, 'ω shared shift SD', 'interval half-width');
  for (const g of [0, 0.5, 1]) { const q = s('text', { x: X(g), y: B + 14, 'text-anchor': 'middle', fill: GREY, 'font-size': 10.5 }); q.textContent = g.toFixed(1); svg.append(q); }
  for (const g of [0, 0.5, 1]) { const w = ymax * g * 0.9; const q = s('text', { x: L - 5, y: Y(w) + 4, 'text-anchor': 'end', fill: GREY, 'font-size': 10 }); q.textContent = w.toFixed(2); svg.append(q); }
  if (!isT) { const q = s('text', { x: R - 4, y: Y(wAt(1, 'T')) - 6, 'text-anchor': 'end', fill: TEAL, 'font-size': 10.5 }); q.textContent = 'Student-t intercept'; svg.append(q); }
  return { cover: null, called: null, hw, N: 0 };
}
const compare = {
  id: 'compare',
  compare: true,
  title: 'Model comparison: Student-t intercept vs orthogonal patient effect',
  ref: 'Comparison of the two cellular models',
  blurb: 'Both models see the same simulated cohorts. The only difference is how each treats the shift that all signatures of one cell type share within a patient. The Student-t intercept model has a term for it. The orthogonal model has one patient effect for all features, so it treats the signatures as more independent than they are and draws intervals that are too narrow. Change any slider and both graphs move together.',
  note: 'Exact arithmetic for the cell-type mean of K signatures: the correct variance per patient is ω² + σ²/K, and a model that ignores the shared shift uses (ω² + σ²)/K. Both panels reuse the same 40 draws of the sampling error, so the estimates are identical and only the intervals differ. This isolates one mechanism and does not reproduce the full simulation, where the orthogonal model covered the true effects 50 to 66% of the time and the Student-t intercept model stayed near nominal.',
  models: [
    { key: 'T', title: 'Student-t intercept model', tag: 'primary', terms: 'T',
      eq: '<span data-term="mu">μ<sub>ij</sub></span> = α<sub>j</sub> + β<sub>j</sub>BMI<sub>i</sub> + … + <span data-term="v">v<sub>i,c(j)</sub></span>' },
    { key: 'O', title: 'Orthogonal patient-effect model', tag: 'comparison', terms: 'O',
      eq: '<span data-term="mu">μ<sub>ij</sub></span> = α<sub>j</sub> + β<sub>j</sub>BMI<sub>i</sub> + … + <span data-term="u">u<sub>i</sub></span>, <span data-term="orth">Xᵀu = 0</span>' }
  ],
  terms: {
    T: {
      mu: ['Fitted score', 'The expected score of signature j in patient i, from BMI, age, sex and the patient shift.'],
      v: ['Patient × cell-type intercept', 'One shift per patient for each cell type, shared by every signature of that cell type. It takes the shared variation ω² out of the signature noise, so a cell-type mean gets an interval of the right width.']
    },
    O: {
      mu: ['Fitted score', 'The expected score of signature j in patient i, from BMI, age, sex and the single patient effect.'],
      u: ['One patient effect for all features', 'A single shift per patient for every signature in every cell type. It has no way to represent a shift that belongs to one cell type, so signatures of that cell type look more independent than they are.'],
      orth: ['Orthogonality constraint', 'u is held uncorrelated with BMI, age and sex, so it can only absorb what the covariates do not explain.']
    }
  },
  params: [
    { k: 'beta', label: 'true cell-type slope', min: -0.3, max: 0.6, step: 0.01, value: 0.10, colour: MAGENTA, hint: 'The effect the cohorts are simulated around, per BMI SD. Set to zero to count false calls.' },
    { k: 'om', label: 'ω<sub>c</sub> shared shift SD', min: 0, max: 1, step: 0.01, value: 0.35, colour: TEAL, hint: 'Size of the shift shared by all signatures of a cell type in a patient.' },
    { k: 'sg', label: 'σ residual SD', min: 0.2, max: 1.5, step: 0.01, value: 0.65, colour: GREY, hint: 'Independent noise of each signature.' },
    { k: 'K', label: 'K signatures in the cell type', min: 1, max: 60, step: 1, value: 35, colour: BLUE, fmt: x => String(Math.round(x)), hint: 'More signatures help little once the shared shift dominates.' },
    { k: 'n', label: 'patients', min: 40, max: 300, step: 1, value: 127, colour: BLUE, fmt: x => String(Math.round(x)), hint: 'The study has 127.' },
    { k: 'seed', label: 'set of repeated cohorts', min: 1, max: 20, step: 1, value: 1, colour: GREY, fmt: x => '#' + Math.round(x), hint: 'A different batch of 40 simulated cohorts.' },
    { k: 'view', label: 'graph', min: 0, max: 1, step: 1, value: 0, colour: GREY, fmt: x => (x ? 'width vs ω' : 'repeated cohorts'), hint: 'Switch between the repeated-cohort intervals and the interval width as the shared shift grows.' }
  ],
  presets: [
    { label: 'Reference', v: { beta: 0.10, om: 0.35, sg: 0.65, K: 35, n: 127, seed: 1 } },
    { label: 'No true effect', v: { beta: 0, om: 0.35, sg: 0.65, K: 35, n: 127, seed: 1 } },
    { label: 'No shared shift', v: { beta: 0.10, om: 0, sg: 0.65, K: 35, n: 127, seed: 1 } },
    { label: 'Large shared shift', v: { beta: 0.10, om: 0.8, sg: 0.65, K: 35, n: 127, seed: 1 } },
    { label: 'Few signatures', v: { beta: 0.10, om: 0.35, sg: 0.65, K: 4, n: 127, seed: 1 } },
    { label: 'Many patients', v: { beta: 0.10, om: 0.35, sg: 0.65, K: 35, n: 300, seed: 1 } }
  ],
  init() {
    const r = lcg(777001);
    return { z: Array.from({ length: 20 }, () => Array.from({ length: 40 }, () => r.n())) };
  },
  panel(key, svg, v, mem) {
    const st = cmpStats(v), out = cmpPanel(svg, v, mem, key), isT = key === 'T';
    const exp = isT ? 0.95 : 2 * Phi(1.959964 * st.ratio) - 1;
    const cards = [
      ['interval half-width', '±' + out.hw.toFixed(3), isT ? 'correct for the shared shift' : (st.ratio * 100).toFixed(0) + '% of the correct width', true],
      ['expected coverage', (exp * 100).toFixed(0) + '%', 'of a nominal 95% interval']
    ];
    if (out.cover !== null) cards.push([v.beta === 0 ? 'false calls in 40' : 'intervals missing the truth', String(v.beta === 0 ? out.called : out.N - out.cover), v.beta === 0 ? 'interval excludes zero' : 'of 40 repeated cohorts']);
    else cards.push(['K<sub>eff</sub>', st.keff.toFixed(1), 'of ' + st.K + ' signatures']);
    return { live: isT ? `K = ${st.K} signatures worth <b>${st.keff.toFixed(1)}</b> independent ones · half-width <b>±${out.hw.toFixed(3)}</b>` : `treats them as nearly independent · half-width <b>±${out.hw.toFixed(3)}</b>`, cards };
  },
  verdict(v) {
    const st = cmpStats(v), exp = 2 * Phi(1.959964 * st.ratio) - 1;
    return st.ratio < 0.92
      ? ['Same data, different intervals', `With a shared shift of ${v.om.toFixed(2)} across ${st.K} signatures the orthogonal model draws intervals ${(st.ratio * 100).toFixed(0)}% as wide as the Student-t intercept model, so a nominal 95% interval covers the truth about ${(exp * 100).toFixed(0)}% of the time. The estimates are identical; the difference is how much uncertainty each model admits to.`]
      : ['The two models agree here', 'With little shared shift, or very few signatures, treating signatures as independent costs almost nothing, so the two sets of intervals are nearly the same width. Raise ω or K to separate them.'];
  }
};

/* ---------- widget: protein minus RNA, Equation 8 ---------- */
const contrast = {
  id: 'contrast',
  title: 'Slope difference versus co-variation',
  eq: '<span data-term="delta">Δ = β<sub>protein</sub> − β<sub>RNA</sub></span> (within-plex, from joint draws) · <span data-term="rho">ρ = residual correlation</span>',
  terms: {
    delta: ['Slope difference', 'Protein slope minus RNA slope, taken draw by draw from the joint posterior so that the correlation between the two slopes is respected.'],
    rho: ['Residual correlation', 'How much the protein and RNA scores of one gene move together between patients beyond BMI. Strong co-variation can sit alongside no slope difference, and the reverse.']
  },
  ref: 'Equation 8, Figure 5',
  blurb: 'Two gene lists can track each other tightly across patients and still respond differently to BMI. Δ and ρ are separate quantities measuring separate things, and the ECM programme is the case where they disagree. Grey connectors join each patient’s RNA and protein score: short connectors mean high ρ.',
  note: 'The contrast compares layer-specific score SDs per 5 kg/m², not concentrations. It is formed within each retained joint posterior draw and summarised afterwards, so its interval carries the dependence between the two slopes. ρ is the residual correlation after covariates and plex, which is a different question from whether the two BMI slopes agree. 122 distinct measured gene lists were fitted, each on its own, with no hierarchy across programmes.',
  params: [
    { k: 'br',  label: 'β RNA BMI slope', min: -0.8, max: 0.8, step: 0.01, value: 0.34, colour: BLUE, hint: 'Within-plex slope for the RNA layer.' },
    { k: 'bp',  label: 'β protein BMI slope', min: -0.8, max: 0.8, step: 0.01, value: 0.05, colour: MAGENTA, hint: 'Within-plex slope for the protein layer.' },
    { k: 'rho', label: 'ρ residual correlation', min: -0.95, max: 0.95, step: 0.01, value: 0.78, colour: TEAL, hint: 'How the two layers co-vary across patients after adjustment.' }
  ],
  presets: [
    { label: 'ECM programme', v: { br: 0.34, bp: 0.05, rho: 0.78 } },
    { label: 'Concordant', v: { br: 0.30, bp: 0.28, rho: 0.75 } },
    { label: 'Independent layers', v: { br: 0.30, bp: 0.10, rho: 0.03 } },
    { label: 'Opposed', v: { br: 0.30, bp: -0.25, rho: 0.60 } }
  ],
  usesData: true,
  real() {
    if (!DATA) return null;
    return [{ group: 'published programmes, by |\u0394|', items: DATA.rnaprot.rows.map(r => ({
      label: `${r[0].replace('SIG_', '')}  \u0394 ${minus(r[3])}  \u03C1 ${minus(r[4])}`,
      v: { br: r[1], bp: r[2], rho: r[4] }
    })) }];
  },
  init() {
    const c = cohort([52, 57, 18], 0.2, 20260926);
    const r = lcg(606060);
    return { c, e1: Array.from({ length: c.n }, () => r.n()), e2: Array.from({ length: c.n }, () => r.n()) };
  },
  draw(svg, v, mem) {
    const { c, e1, e2 } = mem, sd = Math.sqrt(Math.max(0, 1 - v.rho * v.rho));
    svg.setAttribute('viewBox', '0 0 360 268');
    svg.replaceChildren();
    const L = 46, R = 348, T = 18, B = 214;
    const X = z => L + ((clamp(z, -2.6, 3.4) + 2.6) / 6) * (R - L);
    const lim = 2.2, Y = q => clamp(116 - (q / lim) * 96, T, B);
    for (const g of [-2, -1, 0, 1, 2]) svg.append(s('line', { x1: L, x2: R, y1: Y(g), y2: Y(g), stroke: LINE, 'stroke-dasharray': g ? '3 3' : '' }));
    const step = 2;
    for (let i = 0; i < c.n; i += step) {
      const rr = 0.52 * e1[i], pr = 0.52 * (v.rho * e1[i] + sd * e2[i]);
      const yr = v.br * c.z[i] + rr, yp = v.bp * c.z[i] + pr;
      svg.append(s('line', { x1: X(c.z[i]), x2: X(c.z[i]), y1: Y(yr), y2: Y(yp), stroke: GREY, 'stroke-width': 0.8, opacity: 0.35 }));
      svg.append(s('circle', { cx: X(c.z[i]), cy: Y(yr), r: 2.3, fill: BLUE, opacity: 0.55 }));
      svg.append(s('circle', { cx: X(c.z[i]), cy: Y(yp), r: 2.3, fill: MAGENTA, opacity: 0.55 }));
    }
    svg.append(s('line', { x1: X(-2.6), y1: Y(v.br * -2.6), x2: X(3.4), y2: Y(v.br * 3.4), stroke: BLUE, 'stroke-width': 2.4 }));
    svg.append(s('line', { x1: X(-2.6), y1: Y(v.bp * -2.6), x2: X(3.4), y2: Y(v.bp * 3.4), stroke: MAGENTA, 'stroke-width': 2.4 }));
    const l1 = s('text', { x: R - 4, y: Y(v.br * 3.4) - 7, 'text-anchor': 'end', fill: BLUE, 'font-size': 11, 'font-weight': 700 });
    l1.textContent = 'RNA'; svg.append(l1);
    const l2 = s('text', { x: R - 4, y: Y(v.bp * 3.4) + 14, 'text-anchor': 'end', fill: MAGENTA, 'font-size': 11, 'font-weight': 700 });
    l2.textContent = 'protein'; svg.append(l2);
    axes(svg, L, R, T, B, 'standardized BMI', 'standardized score');
    for (const g of [-2, 0, 2]) {
      const t = s('text', { x: X(g), y: B + 15, 'text-anchor': 'middle', fill: GREY, 'font-size': 11 });
      t.textContent = String(g); svg.append(t);
    }
    const key = s('text', { x: L, y: 262, fill: GREY, 'font-size': 10.5 });
    key.textContent = 'grey connectors join one patient’s two layers'; svg.append(key);

    const d = v.bp - v.br;
    let verdict;
    if (Math.abs(d) > 0.15 && v.rho > 0.5) verdict = ['Co-varying but discordant in BMI', 'This is the ECM pattern. The layers track each other closely across patients, yet the protein slope is well below the RNA slope. Strong ρ is no evidence that the two layers agree about BMI.'];
    else if (Math.abs(d) <= 0.15 && v.rho > 0.5) verdict = ['Concordant on both counts', 'The layers co-vary and their BMI slopes agree. Nothing here distinguishes Δ from ρ, which is exactly why the two are easy to conflate.'];
    else if (v.rho < 0.2) verdict = ['Independent layers', 'Low residual correlation: patients deviate separately in RNA and protein. Δ can still be anything — the two quantities are unrelated.'];
    else verdict = ['Opposed slopes', 'The layers move in opposite directions with BMI while still sharing residual variation.'];
    return {
      eq: `Δ = ${minus(v.bp)} − ${minus(v.br)} = <b>${minus(d)}</b> · ρ = ${minus(v.rho)}`,
      cards: [
        ['Δ protein − RNA', minus(d), 'slope difference', true],
        ['ρ residual correlation', minus(v.rho), 'co-variation'],
        ['β RNA', minus(v.br), 'within-plex'],
        ['β protein', minus(v.bp), 'within-plex']
      ],
      verdict
    };
  }
};

/* ---------- widget: the simulation, generating model and projection target ---------- */
const COUNTS = { '64': [26, 29, 9], '127': [52, 57, 18], '254': [104, 114, 36], 'obese 10': [52, 57, 10], 'balanced': [42, 42, 43] };
const COUNT_KEYS = Object.keys(COUNTS);
const SHAPES = ['null', 'linear', 'threshold', 'plateau'];
const ERRORS = ['student_t', 'normal', 'contaminated', 'heteroskedastic'];

const simulation = {
  id: 'simulation',
  title: 'The simulation and its target',
  eq: '<span data-term="B">B*</span> = (<span data-term="X">XᵀX</span>)⁻¹<span data-term="X">Xᵀ</span><span data-term="H">H</span> — the projection of the generating mean onto the fitted design',
  terms: {
    B: ['Target slope', 'What a fitted design is actually aiming at when the truth is not linear. Coverage is judged against this, not against the raw amplitude.'],
    X: ['Fitted design', 'Continuous BMI or the three WHO groups, plus age and sex.'],
    H: ['Generating mean', 'The true curve used to simulate the data: linear, threshold, plateau and so on.']
  },
  ref: 'Simulation design, Figure 2',
  blurb: 'Every scenario in the notebook, and the question they were built to answer: what is a fit even aiming at when the truth is not linear? The grey curve is the generating mean H. The two coloured fits are its least-squares projections onto the continuous and categorical designs. Those projections, not the raw amplitude, are what coverage is measured against.',
  note: 'Reproduces the notebook’s generating model: BMI from three clipped Gaussians, age confounded with BMI, logistic sex, a random quarter of features active with lognormal magnitudes and a shared sign per cell type, residual scales centred at 0.65, and Student-t errors with a shared within-cell component. Projecting a threshold or plateau onto a linear BMI term gives a nonzero but attenuated slope, so a continuous fit can hold nominal coverage for its own target while detecting step changes less often than the categorical contrasts. This draws one dataset; the paper averages 10 replicates across 20 condition-size combinations, 400 model assessments in all.',
  params: [
    { k: 'shape',  label: 'effect shape', min: 0, max: 3, step: 1, value: 1, colour: MAGENTA, fmt: x => SHAPES[x], hint: 'null, linear, threshold (obese only) or plateau (overweight or obese).' },
    { k: 'effect', label: 'A effect magnitude', min: 0, max: 0.6, step: 0.05, value: 0.20, colour: MAGENTA, hint: 'The notebook uses 0.10, 0.20 and 0.50.' },
    { k: 'counts', label: 'cohort', min: 0, max: 4, step: 1, value: 1, colour: TEAL, fmt: x => COUNT_KEYS[x], hint: 'Group sizes: normal, overweight, obese.' },
    { k: 'rho',    label: 'ρ within-cell-type residual correlation', min: 0, max: 0.9, step: 0.05, value: 0, colour: BLUE, hint: 'The correlated scenario uses 0.6. This is what drove the undercoverage.' },
    { k: 'df',     label: 'ν generating degrees of freedom', min: 2.2, max: 20, step: 0.2, value: 8, colour: GREY, fmt: x => x.toFixed(1), hint: 'Heavy-tails scenario uses 3.' },
    { k: 'err',    label: 'error type', min: 0, max: 3, step: 1, value: 0, colour: GREY, fmt: x => ERRORS[x].replace('_', '-'), hint: 'Contaminated multiplies 3% of residuals by five.' },
    { k: 'conf',   label: 'age–BMI confounding', min: 0, max: 0.8, step: 0.05, value: 0.2, colour: BLUE, hint: 'Reference 0.2; the confounding scenario uses 0.65.' },
    { k: 'rep',    label: 'replicate', min: 1, max: 10, step: 1, value: 1, colour: TEAL, fmt: x => '#' + x, hint: 'A new dataset from the same condition.' }
  ],
  presets: [
    { label: 'Reference', v: { shape: 1, effect: 0.20, counts: 1, rho: 0, df: 8, err: 0, conf: 0.2 } },
    { label: 'Null', v: { shape: 0, effect: 0, counts: 1, rho: 0, df: 8, err: 0, conf: 0.2 } },
    { label: 'Threshold', v: { shape: 2, effect: 0.20, counts: 1, rho: 0, df: 8, err: 0, conf: 0.2 } },
    { label: 'Plateau', v: { shape: 3, effect: 0.20, counts: 1, rho: 0, df: 8, err: 0, conf: 0.2 } },
    { label: 'Correlated 0.6', v: { shape: 1, effect: 0.20, counts: 1, rho: 0.6, df: 8, err: 0, conf: 0.2 } },
    { label: 'Heavy tails', v: { shape: 1, effect: 0.20, counts: 1, rho: 0, df: 3, err: 0, conf: 0.2 } },
    { label: 'Contaminated', v: { shape: 1, effect: 0.20, counts: 1, rho: 0, df: 8, err: 2, conf: 0.2 } },
    { label: 'Obese 10', v: { shape: 1, effect: 0.20, counts: 3, rho: 0, df: 8, err: 0, conf: 0.2 } },
    { label: 'Confounded', v: { shape: 1, effect: 0.20, counts: 1, rho: 0, df: 8, err: 0, conf: 0.65 } }
  ],
  usesData: true,
  scenarioKey(v) {
    const shape = SHAPES[Math.round(v.shape)];
    if (shape === 'null') return 'null';
    if (Math.round(v.err) === 2) return 'contaminated_0.20';
    if (v.df <= 3.2) return 'heavy_tails';
    if (v.rho >= 0.55) return 'correlated_0.20';
    const ci = Math.round(v.counts);
    if (ci === 0) return 'patients_64';
    if (ci === 2) return 'patients_254';
    if (ci === 3) return 'obese_10';
    if (shape === 'threshold') return 'threshold_0.20';
    if (shape === 'plateau') return 'plateau_0.20';
    return 'linear_' + (v.effect <= 0.125 ? '0.10' : v.effect >= 0.35 ? '0.50' : '0.20');
  },
  init() { return {}; },
  draw(svg, v) {
    const counts = COUNTS[COUNT_KEYS[Math.round(v.counts)]];
    const shape = SHAPES[Math.round(v.shape)], err = ERRORS[Math.round(v.err)];
    const rep = Math.round(v.rep);
    const c = cohort(counts, v.conf, 20260926 + rep * 7919);
    const r = lcg(424242 + rep * 104729);
    // one representative active feature: beta = A * lognormal magnitude * sign
    const beta = v.effect * Math.exp(-0.02 + 0.2 * r.n());
    const gA = 0.12 * r.n() + 0.08 * r.n(), gS = 0.12 * r.n() + 0.08 * r.n();
    const basis = i => shape === 'null' ? 0 : shape === 'linear' ? c.z[i]
      : shape === 'threshold' ? (c.group[i] === 2 ? 1 : 0) : (c.group[i] !== 0 ? 1 : 0);
    const H = c.z.map((_, i) => basis(i) * beta + c.age[i] * gA + c.sex[i] * gS);
    const scale = Math.exp(Math.log(0.65) + 0.15 * r.n() + 0.25 * r.n());
    const shared = Array.from({ length: c.n }, () => r.n());
    const y = H.map((hi, i) => {
      const indep = err === 'normal' ? r.n() : r.t(v.df);
      let e = Math.sqrt(1 - v.rho) * indep + Math.sqrt(v.rho) * shared[i];
      let sc = scale;
      if (err === 'contaminated' && r.u() < 0.03) e *= 5;
      if (err === 'heteroskedastic') sc *= Math.exp(0.35 * c.z[i]);
      return hi + e * sc;
    });
    // projection targets, exactly as targets() does
    const Xc = c.z.map((z, i) => [1, z, c.age[i], c.sex[i]]);
    const Xk = c.z.map((_, i) => [1, c.group[i] === 1 ? 1 : 0, c.group[i] === 2 ? 1 : 0, c.age[i], c.sex[i]]);
    const bc = lstsq(Xc, H), bk = lstsq(Xk, H);
    const contCoef = bc[1], ow = bk[1], ob = bk[2];

    svg.setAttribute('viewBox', '0 0 360 292');
    svg.replaceChildren();
    const L = 46, R = 348, T = 16, B = 214;
    const bmin = 18, bmax = 46;
    const X = b => L + ((clamp(b, bmin, bmax) - bmin) / (bmax - bmin)) * (R - L);
    const lim = Math.max(1.2, 2.4 * scale);
    const Y = q => clamp(115 - (q / lim) * 95, T, B);
    for (const g of [-1, 0, 1]) svg.append(s('line', { x1: L, x2: R, y1: Y(g), y2: Y(g), stroke: LINE, 'stroke-dasharray': g ? '3 3' : '' }));
    for (const b of [25, 30]) svg.append(s('line', { x1: X(b), x2: X(b), y1: T, y2: B, stroke: LINE, 'stroke-dasharray': '2 4' }));
    const cols = [TEAL, BLUE, MAGENTA];
    for (let i = 0; i < c.n; i++) svg.append(s('circle', { cx: X(c.bmi[i]), cy: Y(y[i]), r: 2.3, fill: cols[c.group[i]], opacity: 0.4 }));
    // generating mean H, averaged within narrow BMI bins
    const bins = 14, acc = Array.from({ length: bins }, () => [0, 0]);
    c.bmi.forEach((b, i) => { const k = clamp(Math.floor((b - bmin) / (bmax - bmin) * bins), 0, bins - 1); acc[k][0] += H[i]; acc[k][1]++; });
    let d = '';
    acc.forEach((a, k) => { if (!a[1]) return; const bx = bmin + (k + 0.5) / bins * (bmax - bmin); d += (d ? 'L' : 'M') + X(bx) + ' ' + Y(a[0] / a[1]); });
    svg.append(s('path', { d, fill: 'none', stroke: INK, 'stroke-width': 2.6 }));
    // continuous projection, drawn at mean age and sex
    const zAt = b => (b - (c.bmi.reduce((a, q) => a + q, 0) / c.n)) / c.sd;
    svg.append(s('line', { x1: X(bmin), y1: Y(bc[0] + contCoef * zAt(bmin)), x2: X(bmax), y2: Y(bc[0] + contCoef * zAt(bmax)),
      stroke: MAGENTA, 'stroke-width': 2.2, 'stroke-dasharray': '7 4' }));
    // categorical projection, three steps
    const seg = [[bmin, 25, bk[0]], [25, 30, bk[0] + ow], [30, bmax, bk[0] + ob]];
    for (const [a, b2, val] of seg) svg.append(s('line', { x1: X(a), y1: Y(val), x2: X(b2), y2: Y(val), stroke: TEAL, 'stroke-width': 2.6 }));
    axes(svg, L, R, T, B, 'BMI (kg/m²)', 'score');
    for (const b of [20, 25, 30, 35, 40, 45]) {
      const t = s('text', { x: X(b), y: B + 15, 'text-anchor': 'middle', fill: GREY, 'font-size': 10.5 });
      t.textContent = String(b); svg.append(t);
    }
    const sw = (x, y, col, dash, wd) => s('line', { x1: x, x2: x + 16, y1: y - 3.5, y2: y - 3.5, stroke: col, 'stroke-width': wd, 'stroke-dasharray': dash || '' });
    const lg = [[L, INK, '', 2.6, 'generating mean H'], [L + 112, MAGENTA, '5 3', 2.2, 'continuous target'], [L + 214, TEAL, '', 2.6, 'categorical target']];
    for (const [x, col, dash, wd, txt] of lg) {
      svg.append(sw(x, 266, col, dash, wd));
      const q = s('text', { x: x + 21, y: 266, fill: INK, 'font-size': 10 }); q.textContent = txt; svg.append(q);
    }
    const k4 = s('text', { x: L, y: 284, fill: GREY, 'font-size': 10 });
    k4.textContent = `${c.n} patients · σⱼ ${scale.toFixed(2)} · βⱼ ${minus(beta)}`;
    svg.append(k4);

    const atten = beta ? contCoef / beta : 0;
    // published operating characteristics for the resolved scenario
    const skey = simulation.scenarioKey(v);
    let real = null;
    if (DATA && DATA.simulation) {
      const B = DATA.simulation.rows.find(r => r[0] === skey && r[2] === 'bmi' && r[3] === 'B');
      const O = DATA.simulation.rows.find(r => r[0] === skey && r[2] === 'bmi' && r[3] === 'O');
      if (B && O) real = { bRmse: B[4], bCov: B[5], oRmse: O[4], oCov: O[5] };
    }
    if (real) {
      svg.setAttribute('viewBox', '0 0 360 356');
      const y0 = 328, SC = 600, XT = 176, bw = r => Math.min(120, Math.max(2, r * SC));
      svg.append(s('line', { x1: 46, x2: 348, y1: 296, y2: 296, stroke: LINE }));
      const k = s('text', { x: 46, y: 312, fill: INK, 'font-size': 10, 'font-weight': 700 });
      k.textContent = 'Published simulation: ' + skey; svg.append(k);
      [['Bayes', real.bRmse, real.bCov, TEAL, y0], ['OLS', real.oRmse, real.oCov, BLUE, y0 + 14]].forEach(([nm, r, cv, col, y]) => {
        svg.append(s('rect', { x: 46, y: y - 8, width: bw(r), height: 8, rx: 2, fill: col }));
        const q = s('text', { x: XT, y, fill: GREY, 'font-size': 10 });
        q.textContent = `${nm} RMSE ${r.toFixed(3)} · coverage ${(cv * 100).toFixed(1)}%`; svg.append(q);
      });
    } else svg.setAttribute('viewBox', '0 0 360 292');

    let verdict;
    if (shape === 'null') verdict = ['Nothing to find', 'Both targets are zero. This is the condition that measures the false-discovery proportion, not power.'];
    else if (shape === 'linear') verdict = ['Targets agree', 'With a genuinely linear effect the continuous projection recovers β directly, and the categorical contrasts describe the same line sampled at three points.'];
    else verdict = ['The targets differ, and that is the point',
      `A ${shape} effect projected onto a linear BMI term gives ${minus(contCoef)} against a generating amplitude of ${minus(beta)} — attenuated to ${(atten * 100).toFixed(0)}%. Coverage for the continuous fit is judged against that projection, not against the amplitude, which is why a continuous fit can stay nominally calibrated while detecting the step less often than the categorical contrasts do.`];
    if (real) {
      verdict = [real.bCov < 0.90 ? 'Better estimates, worse intervals' : 'Published result for this scenario',
        `In the published simulation this scenario gave Bayes RMSE ${real.bRmse.toFixed(3)} against OLS ${real.oRmse.toFixed(3)}, with coverage ${(real.bCov * 100).toFixed(1)}% against ${(real.oCov * 100).toFixed(1)}%. ` +
        (real.bCov < 0.90
          ? 'The pooled model estimates more accurately and still undercovers. The two criteria come apart, which is the whole reason both are reported.'
          : 'Accuracy and calibration agree in this scenario.')];
    }
    return {
      eq: `continuous target <b>${minus(contCoef)}</b> · overweight ${minus(ow)} · obese ${minus(ob)} · obese−overweight ${minus(ob - ow)}`,
      cards: [
        ['continuous target', minus(contCoef), 'projection of H', true],
        ['obese − normal', minus(ob), 'categorical'],
        ['obese − overweight', minus(ob - ow), 'derived in draws'],
        real ? ['Bayes vs OLS RMSE', real.bRmse.toFixed(3) + ' / ' + real.oRmse.toFixed(3), 'published']
             : ['generating βⱼ', minus(beta), shape === 'linear' ? 'recovered directly' : 'attenuated to ' + (atten * 100).toFixed(0) + '%']
      ],
      verdict
    };
  }
};

/* ---------- published data ---------- */
let DATA = null, dataPromise = null;
export function loadExplainerData() {
  if (!dataPromise) {
    dataPromise = fetch(new URL('./explainer-data.json', import.meta.url))
      .then(r => (r.ok ? r.json() : null))
      .then(d => { DATA = d; return d; })
      .catch(() => null);
  }
  return dataPromise;
}
const DATASRC = 'Published posterior summaries from the Atlas database.';

/* ---------- widget: what the screen intended and what it returned ---------- */
const funnel = {
  id: 'screen',
  title: 'What the screen returned',
  eq: '<span data-term="marg">marginal call</span> → <span data-term="cond">conditional call</span>',
  terms: {
    marg: ['Marginal call', 'Negative, uncertain or positive from the modification alone, ignoring its parent protein, by whether the 95% HDI excludes zero.'],
    cond: ['Conditional call', 'The same call after conditioning on the parent protein. Features can change class between the two.']
  },
  ref: 'PTM screen, Figure 6',
  blurb: 'The screen set out to find modified forms associated with BMI. Conditioning on parent protein then rebuilt the candidate list almost from scratch. These are the published counts: how many features entered, how many survived each gate, and how the marginal and conditional calls disagree.',
  note: 'A feature qualified with finite protein and PTM values for at least 102 patients across at least 16 plexes, nonzero variance in both responses and an available parent protein. Candidates require a numerically adequate primary fit and a conditional 95% HDI excluding zero. No hierarchy pools features within genes, signatures or pathways, and there is no multiplicity control across the screen. ' + DATASRC,
  usesData: true,
  params: [
    { k: 'assay', label: 'assay', min: 0, max: 1, step: 1, value: 1, colour: TEAL, fmt: x => (x ? 'glyco' : 'phospho'), hint: 'Phosphosites or glycopeptides.' },
    { k: 'logx',  label: 'funnel scale', min: 0, max: 1, step: 1, value: 1, colour: GREY, fmt: x => (x ? 'log' : 'linear'), hint: 'The first stage is three orders of magnitude above the last.' }
  ],
  init() { return {}; },
  draw(svg, v) {
    const assay = v.assay ? 'glyco' : 'phospho';
    svg.setAttribute('viewBox', '0 0 420 340');
    svg.replaceChildren();
    if (!DATA) {
      const t = s('text', { x: 210, y: 140, 'text-anchor': 'middle', fill: GREY, 'font-size': 12 });
      t.textContent = 'Loading published counts…'; svg.append(t);
      return { eq: 'Loading…', cards: [], verdict: ['Loading', 'Fetching the published screen counts.'] };
    }
    const col = assay === 'glyco' ? 2 : 1;
    const stages = DATA.funnel.rows.map(r => [r[0], r[col]]);
    const genes = DATA.funnel.genes[assay];
    const tr = DATA.transitions.rows.filter(r => r[0] === assay);
    const get = (m, c) => (tr.find(r => r[1] === m && r[2] === c) || [0, 0, 0, 0])[3];
    const dirs = ['negative', 'uncertain', 'positive'];
    const marginalCalled = get('negative','negative') + get('negative','uncertain') + get('positive','positive') + get('positive','uncertain');
    const conditionalCalled = get('negative','negative') + get('positive','positive') + get('uncertain','negative') + get('uncertain','positive');
    const kept = get('negative','negative') + get('positive','positive');
    const lost = get('negative','uncertain') + get('positive','uncertain');
    const gained = get('uncertain','negative') + get('uncertain','positive');

    // funnel bars
    const L = 188, R = 366, T = 14, rowH = 24;
    const max = stages[0][1];
    const w = n => {
      if (!v.logx) return (n / max) * (R - L);
      return (Math.log10(Math.max(1, n)) / Math.log10(max)) * (R - L);
    };
    stages.forEach(([name, n], i) => {
      const y = T + i * rowH;
      svg.append(s('rect', { x: L, y, width: Math.max(2, w(n)), height: rowH - 7, rx: 2,
        fill: i === stages.length - 1 ? TEAL : BLUE, opacity: i === stages.length - 1 ? 0.95 : 0.3 + 0.12 * i }));
      const lab = s('text', { x: L - 7, y: y + rowH / 2, 'text-anchor': 'end', fill: INK, 'font-size': 10.5 });
      lab.textContent = name; svg.append(lab);
      const num = s('text', { x: L + Math.max(2, w(n)) + 6, y: y + rowH / 2, fill: GREY, 'font-size': 10.5 });
      num.textContent = n.toLocaleString('en'); svg.append(num);
    });
    const gl = s('text', { x: L, y: T + stages.length * rowH + 10, fill: GREY, 'font-size': 10.5 });
    gl.textContent = `${genes.candidate} candidate genes of ${genes.eligible} eligible`; svg.append(gl);

    // transition matrix
    const MT = T + stages.length * rowH + 40, cw = 62, ch = 40, MX = 196;
    const head = s('text', { x: MX + cw * 1.5, y: MT - 8, 'text-anchor': 'middle', fill: GREY, 'font-size': 10.5 });
    head.textContent = 'conditional call →'; svg.append(head);
    const side = s('text', { x: 20, y: MT + 18 + ch * 1.5, 'text-anchor': 'middle', fill: GREY, 'font-size': 10.5,
      transform: `rotate(-90 20 ${MT + 18 + ch * 1.5})` });
    side.textContent = '← marginal call'; svg.append(side);
    dirs.forEach((c, j) => {
      const t = s('text', { x: MX + j * cw + (cw - 3) / 2, y: MT + 9, 'text-anchor': 'middle', fill: GREY, 'font-size': 10 });
      t.textContent = c; svg.append(t);
    });
    dirs.forEach((m, i) => {
      const t = s('text', { x: MX - 8, y: MT + 18 + i * ch + ch / 2, 'text-anchor': 'end', fill: GREY, 'font-size': 10 });
      t.textContent = m; svg.append(t);
      dirs.forEach((c, j) => {
        const n = get(m, c), x = MX + j * cw, y = MT + 18 + i * ch;
        const diag = m === c, call = c !== 'uncertain';
        const fill = n === 0 ? 'transparent'
          : diag && call ? 'rgba(14,124,123,0.75)'
          : call ? 'rgba(216,27,96,0.55)'
          : m !== 'uncertain' ? 'rgba(44,90,160,0.45)' : 'rgba(0,0,0,0.05)';
        svg.append(s('rect', { x, y, width: cw - 4, height: ch - 4, rx: 3, fill,
          stroke: LINE, 'stroke-width': 1 }));
        const t2 = s('text', { x: x + (cw - 4) / 2, y: y + (ch - 4) / 2 + 4, 'text-anchor': 'middle',
          'font-size': 11, 'font-weight': call || diag ? 700 : 400,
          fill: n && (call || (diag && call)) ? '#fff' : INK });
        t2.textContent = n ? String(n) : '·'; svg.append(t2);
      });
    });
    const key = s('text', { x: MX, y: MT + 18 + 3 * ch + 12, fill: GREY, 'font-size': 10 });
    key.textContent = 'teal kept · blue lost · pink gained'; svg.append(key);

    return {
      eq: `marginal screen would call <b>${marginalCalled}</b> · conditional screen calls <b>${conditionalCalled}</b> · only ${kept} are the same features`,
      cards: [
        ['kept', String(kept), 'called both ways', true],
        ['lost to conditioning', String(lost), 'tracked parent protein'],
        ['gained by conditioning', String(gained), 'hidden marginally'],
        ['turnover', Math.round(100 * (1 - kept / Math.max(1, conditionalCalled))) + '%', 'of the final list is new']
      ],
      verdict: ['The final list is not the list we set out with',
        `Of ${marginalCalled} ${assay} features a marginal screen would have called, ${kept} survive conditioning and ${lost} drop out. A further ${gained} appear that showed nothing marginally. ` +
        (assay === 'phospho'
          ? 'For phosphosites that is near-total turnover: the published 14 candidates share only 6 features with what the marginal screen would have reported.'
          : 'For glycopeptides 39 of the 67 published candidates were already marginal calls, so the overlap is larger but still partial.')]
    };
  }
};

/* ---------- widget: interval widening, and which signals survive ---------- */
const widening = {
  id: 'widening',
  title: 'Which associations survive widening',
  eq: '<span data-term="factor">factor = z<sub>0.975</sub> / z<sub>(1+c)/2</sub></span>, <span data-term="surv">survives if |β| > factor × (distance to the zero-side bound)</span>',
  terms: {
    factor: ['Widening factor', 'How far each interval is stretched so that its coverage matches 95%: 1.960 divided by the normal quantile for the achieved coverage c.'],
    surv: ['Survival rule', 'An association survives if its interval, once widened, still excludes zero. Intervals are asymmetric, so the distance is measured from the posterior mean to the bound on the side of zero.']
  },
  ref: 'Decision rule, simulation coverage',
  blurb: 'Simulated residual correlation among signatures pushed coverage of the nominal 95% intervals down to about 84.6%, so named associations were re-checked after widening each interval by 1.37. This applies that correction to every published cellular slope and shows exactly which ones hold.',
  note: 'Signatures within a cell type share genes, so their residuals are correlated and the intervals are narrower than their nominal level. The widening is an approximate, illustrative correction, not a re-fit. Applied across all three compartments. ' + DATASRC,
  usesData: true,
  params: [
    { k: 'cov',  label: 'achieved coverage', min: 0.70, max: 0.95, step: 0.002, value: 0.846, colour: MAGENTA, fmt: x => (x * 100).toFixed(1) + '%', hint: 'Simulation gave 84.6% for the correlated scenario. 95% means no correction.' },
    { k: 'show', label: 'rows shown', min: 10, max: 40, step: 1, value: 22, colour: GREY, fmt: x => String(Math.round(x)), hint: 'Highest-ranked associations by margin.' }
  ],
  presets: [
    { label: 'Paper: 84.6%', v: { cov: 0.846, show: 22 } },
    { label: 'No correction', v: { cov: 0.95, show: 22 } },
    { label: 'Nominal 93%', v: { cov: 0.93, show: 22 } },
    { label: 'Severe 75%', v: { cov: 0.75, show: 22 } }
  ],
  init() { return {}; },
  draw(svg, v) {
    if (!DATA) {
      svg.setAttribute('viewBox', '0 0 460 120'); svg.replaceChildren();
      const t = s('text', { x: 230, y: 60, 'text-anchor': 'middle', fill: GREY, 'font-size': 12 });
      t.textContent = 'Loading published estimates…'; svg.append(t);
      return { eq: 'Loading…', cards: [], verdict: ['Loading', 'Fetching published cellular slopes.'] };
    }
    const W = DATA.widening, f = widenFactor(v.cov);
    const n = Math.round(v.show);
    const rows = W.rows.slice(0, n);
    const survive = W.rows.filter(r => r[5] > f).length;
    const extra = W.excludes_zero - W.rows.length; // below the stored cut, never survives f>1.20
    const totalSurvive = f >= 1.2 ? survive : survive + extra;

    const rowH = 19, top = 26, H = top + rows.length * rowH + 38;
    svg.setAttribute('viewBox', `0 0 470 ${H}`);
    svg.replaceChildren();
    const lim = 0.28, L = 232, R = 458;
    const X = b => L + (clamp(b, -lim, lim) / lim) * ((R - L) / 2) + (R - L) / 2;
    svg.append(s('line', { x1: X(0), x2: X(0), y1: top - 8, y2: top + rows.length * rowH, stroke: INK }));
    for (const t of [-0.2, 0.2]) {
      svg.append(s('line', { x1: X(t), x2: X(t), y1: top - 8, y2: top + rows.length * rowH, stroke: LINE, 'stroke-dasharray': '3 3' }));
      const tl = s('text', { x: X(t), y: top + rows.length * rowH + 15, 'text-anchor': 'middle', fill: GREY, 'font-size': 10.5 });
      tl.textContent = minus(t); svg.append(tl);
    }
    rows.forEach((r, i) => {
      const [ct, sig, est, lo, hi, fmax] = r;
      const y = top + i * rowH + rowH / 2;
      const ok = fmax > f;
      const col = est >= 0 ? MAGENTA : BLUE;
      const nm = s('text', { x: 226, y: y + 4, 'text-anchor': 'end', fill: ok ? INK : GREY, 'font-size': 10.5 });
      nm.textContent = ct + ' · ' + sig.replace(/_/g, ' ');
      if (!ok) nm.setAttribute('opacity', '0.55');
      svg.append(nm);
      const wlo = est - f * (est - lo), whi = est + f * (hi - est);
      svg.append(s('line', { x1: X(wlo), x2: X(whi), y1: y, y2: y, stroke: ok ? col : GREY, 'stroke-width': 1, opacity: 0.4 }));
      svg.append(s('line', { x1: X(lo), x2: X(hi), y1: y, y2: y, stroke: col, 'stroke-width': 2.2, opacity: ok ? 1 : 0.45 }));
      svg.append(s('circle', { cx: X(est), cy: y, r: 3.2, fill: col, opacity: ok ? 1 : 0.45 }));
      if (!ok) svg.append(s('circle', { cx: X(0), cy: y, r: 2.6, fill: 'none', stroke: GREY, 'stroke-width': 1.4 }));
    });
    const cap = s('text', { x: (L + R) / 2, y: H - 8, 'text-anchor': 'middle', fill: GREY, 'font-size': 10.5 });
    cap.textContent = 'BMI slope per 5 kg/m², pale bar = widened interval'; svg.append(cap);

    return {
      eq: `coverage ${(v.cov * 100).toFixed(1)}% → widening ×<b>${f.toFixed(2)}</b> · <b>${totalSurvive}</b> of ${W.excludes_zero} associations still exclude zero`,
      cards: [
        ['exclude zero at 95%', String(W.excludes_zero), 'of ' + W.total + ' signatures'],
        ['survive widening', String(totalSurvive), '×' + f.toFixed(2), true],
        ['lost', String(W.excludes_zero - totalSurvive), Math.round(100 * (1 - totalSurvive / W.excludes_zero)) + '% of calls'],
        ['widening factor', '×' + f.toFixed(2), 'z ratio']
      ],
      verdict: f <= 1.01
        ? ['No correction applied', `At nominal coverage all ${W.excludes_zero} associations stand. The simulation says this is optimistic, because correlated residuals among overlapping signatures make these intervals narrower than 95%.`]
        : [`${totalSurvive} of ${W.excludes_zero} hold`,
           `At the paper's 84.6% coverage the correction is ×1.37 and 27 associations survive — the TAM iron, lipid, lysosomal and repair programmes, CD4+ helper and costimulatory function, and CD8+ costimulation. CD8+ metabolic programming sits at 1.369, just under the threshold, which is why the manuscript calls it sensitive to uncertainty assumptions. Cytotoxicity never excluded zero at any level.`]
    };
  }
};

/* ---------- widget: measurement support ---------- */
const support = {
  id: 'support',
  title: 'Measurement support and what it costs',
  eq: 'profiles with fewer than <span data-term="cut">1,000</span> <span data-term="reads">allocated reads</span> were dropped',
  terms: {
    reads: ['Allocated reads', 'Sequencing reads assigned to one cell type by deconvolution in one patient.'],
    cut: ['Read cutoff', 'The minimum for a patient × cell-type profile to enter the fit in the earlier fits. The sensitivity analysis shows what changes without it.']
  },
  ref: 'Read-depth sensitivity, component eligibility',
  blurb: 'Deconvolution assigns each bulk sample’s reads across cell types. For some types almost nothing is assigned in some patients, and a score built on 19 reads is noise. Each point is one component: the horizontal axis is its median allocated reads, the bar runs down to the worst patient it had.',
  note: 'The exclusion rule does not use BMI, but it changes which patients represent each cell type, so estimates and interval widths were compared before and after. Those before-and-after estimates are not published in the Atlas; what is published, and shown here, are the support metrics that drove the gate. strict_qc_clear passes the strict gate; borderline_supported fails it but is retained by the tiered gate used for the primary analysis; exploratory_zero_heavy fails both. ' + DATASRC,
  usesData: true,
  params: [
    { k: 'cutoff', label: 'read cutoff', min: 1, max: 6, step: 0.1, value: 3, colour: MAGENTA, fmt: x => Math.round(Math.pow(10, x)).toLocaleString('en'), hint: 'The analysis used 1,000 reads per patient-by-cell-type profile.' },
    { k: 'comp',   label: 'compartment', min: 0, max: 2, step: 1, value: 0, colour: TEAL, fmt: x => ['coarse immune', 'fine immune', 'non-immune'][x], hint: 'Three separate deconvolutions of the same bulk data.' }
  ],
  presets: [
    { label: 'Coarse immune', v: { cutoff: 3, comp: 0 } },
    { label: 'Fine immune', v: { cutoff: 3, comp: 1 } },
    { label: 'Non-immune', v: { cutoff: 3, comp: 2 } },
    { label: 'Stricter: 10,000', v: { cutoff: 4, comp: 0 } }
  ],
  init() { return {}; },
  draw(svg, v) {
    if (!DATA) {
      svg.setAttribute('viewBox', '0 0 420 120'); svg.replaceChildren();
      const t = s('text', { x: 210, y: 60, 'text-anchor': 'middle', fill: GREY, 'font-size': 12 });
      t.textContent = 'Loading support metrics…'; svg.append(t);
      return { eq: 'Loading…', cards: [], verdict: ['Loading', 'Fetching component support metrics.'] };
    }
    const comp = ['immune_coarse', 'immune_fine', 'non_immune'][Math.round(v.comp)];
    const rows = DATA.support.rows.filter(r => r[1] === comp)
      .slice().sort((a, b) => a[3] - b[3]);
    const cut = Math.pow(10, v.cutoff);
    const rowH = 17, top = 24, H = top + rows.length * rowH + 44;
    svg.setAttribute('viewBox', `0 0 470 ${H}`);
    svg.replaceChildren();
    const L = 212, R = 452;
    const X = n => L + (Math.log10(Math.max(1, n)) / 7.2) * (R - L);
    for (const d of [1, 3, 5, 7]) {
      svg.append(s('line', { x1: X(Math.pow(10, d)), x2: X(Math.pow(10, d)), y1: top - 6, y2: top + rows.length * rowH, stroke: LINE, 'stroke-dasharray': '3 3' }));
      const t = s('text', { x: X(Math.pow(10, d)), y: top + rows.length * rowH + 15, 'text-anchor': 'middle', fill: GREY, 'font-size': 10 });
      t.textContent = '10' + ['', '', '', '³', '', '⁵', '', '⁷'][d]; svg.append(t);
    }
    svg.append(s('line', { x1: X(cut), x2: X(cut), y1: top - 10, y2: top + rows.length * rowH, stroke: MAGENTA, 'stroke-width': 1.8 }));
    const cl = s('text', { x: X(cut), y: top - 14, 'text-anchor': 'middle', fill: MAGENTA, 'font-size': 10.5, 'font-weight': 700 });
    cl.textContent = 'cutoff'; svg.append(cl);

    const tierCol = { strict_qc_clear: TEAL, borderline_supported: '#B07A12', exploratory_zero_heavy: MAGENTA };
    let below = 0, borderline = 0;
    rows.forEach((r, i) => {
      const [name, , tier, med, min] = r;
      const y = top + i * rowH + rowH / 2;
      if (min < cut) below++;
      if (tier === 'borderline_supported') borderline++;
      const nm = s('text', { x: 206, y: y + 3.5, 'text-anchor': 'end', fill: INK, 'font-size': 10 });
      nm.textContent = name.replace(/_/g, ' ').toLowerCase(); svg.append(nm);
      svg.append(s('line', { x1: X(min), x2: X(med), y1: y, y2: y, stroke: tierCol[tier], 'stroke-width': 1.6, opacity: 0.45 }));
      svg.append(s('circle', { cx: X(min), cy: y, r: 2.4, fill: 'none', stroke: tierCol[tier], 'stroke-width': 1.4 }));
      svg.append(s('circle', { cx: X(med), cy: y, r: 3.4, fill: tierCol[tier] }));
    });
    const xl = s('text', { x: (L + R) / 2, y: H - 10, 'text-anchor': 'middle', fill: GREY, 'font-size': 10.5 });
    xl.textContent = 'allocated reads: open circle = worst patient, filled = median'; svg.append(xl);

    const kept = rows.filter(r => r[9]).length;
    return {
      eq: `${comp.replace('_', ' ')}: <b>${kept}</b> of ${rows.length} components retained · ${below} have at least one patient below the cutoff`,
      cards: [
        ['components', String(rows.length), 'in this compartment'],
        ['retained', String(kept), 'tiered primary gate', true],
        ['borderline', String(borderline), 'fail strict, kept tiered'],
        ['with a patient below cutoff', String(below), 'profiles excluded']
      ],
      verdict: comp === 'immune_coarse'
        ? ['The two cell types carrying the headline findings are the borderline ones',
           'CD8+ T cells and CD4+ T cells both fail the strict gate and are retained only by the tiered one. Their median allocation looks healthy — 189,354 and 548,097 reads — but their worst patients contribute 19 and 57. That is the asymmetry the cutoff exists to handle, and it is also why the CD8+ metabolic association is the one that moves when the uncertainty assumptions change.']
        : comp === 'immune_fine'
          ? ['Fine subsets mostly have nothing to measure',
             'Seven fine T-cell subsets have median allocations in single digits and a median deconvolution fraction of zero. No modelling choice rescues a score built on four reads, so these are reported as exploratory and excluded from the primary analysis.']
          : ['Most non-immune components fail both gates',
             'Tumour epithelial, normal ductal, stellate and adipocyte components are well supported. Fibroblast subtypes, acinar, Schwann and islet components have median allocations under 20 reads, which is why the non-immune associations in the manuscript stay imprecise.']
    };
  }
};


/* ---------- data provenance, shown beside every figure ---------- */
const PROVENANCE = {
  'screen':        ['published', 'Published counts', 'Funnel, transition matrix and gene counts are the study\u2019s own figures.'],
  'widening':      ['published', 'Orthogonal patient-effect estimates', 'All 914 cellular slopes and their HDIs from the orthogonal patient-effect model, with the correction applied live. Student-t intercept estimates are in the RNA collections.'],
  'support':       ['published', 'Published metrics', 'Allocated-read totals and gate tiers per deconvolution component.'],
  'pooling':       ['published', 'Orthogonal patient-effect estimates', 'Twelve slopes from the orthogonal patient-effect model, strongest to weakest survivor. The noise multiplier is a what-if.'],
  'conditioning':  ['mixed', 'Published estimates on request', 'Opens on the paper\u2019s calibration scenarios. Load a real feature from the selector above.'],
  'contrast':      ['mixed', 'Published estimates on request', 'Opens on illustrative values. Load a real programme from the selector above.'],
  'simulation':    ['mixed', 'Published metrics, simulated draw', 'The dataset is generated live from the notebook\u2019s model; the RMSE and coverage bars are published.'],
  'student-t':     ['exact', 'Exact formula', 'The weight curve is the likelihood itself, not a data summary.'],
  'intercept':     ['exact', 'Exact formula', 'Coverage is computed from the variance of a cell-type mean; no patient data are used.'],
  'orth-cellular': ['illustrative', 'Simulated cohort', 'Comparison model. Patient-level data is controlled access, so the 127 patients here are generated; orthogonality of u is exact arithmetic on them.'],
  'compare':       ['exact', 'Exact formula, simulated draws', 'Interval widths are closed-form. The 40 repeated cohorts are seeded draws of the sampling error, identical in both panels.'],
  'cellular':      ['illustrative', 'Simulated cohort', 'Patient-level data is controlled access, so the 127 patients here are generated from the study\u2019s cohort structure.'],
  'patient-effect':['illustrative', 'Simulated cohort', 'Comparison model. Patient-level data is controlled access; the orthogonality shown is exact arithmetic on generated patients.'],
  'plex':          ['illustrative', 'Simulated cohort', 'Patient-level data is controlled access, so plex membership and scores here are generated.'],
  'heatmap':       ['illustrative', 'Illustrative grid', 'The published cell-type by signature grid is too sparse to fill; shrinkage shown is the exact conjugate result.']
};

export const EXPLAINERS = [cellular, orthCellular, compare, intercept, pooling, heatmap, support, patient, plex, contrast, conditioning, funnel, widening, studentT, simulation];

/* ---------- mount ---------- */
function mountCompare(root, spec) {
  const v = {};
  for (const p of spec.params) v[p.k] = p.value;
  const mem = spec.init();
  root.replaceChildren();
  const head = h('div', undefined, { class: 'xp-head' });
  head.append(h('p', spec.ref, { class: 'xp-ref' }), h('p', spec.blurb, { class: 'xp-blurb' }));
  root.append(head);

  const bar = h('div', undefined, { class: 'xp-cmp-controls' });
  const inputs = {};
  for (const p of spec.params) {
    const wrap = h('div', undefined, { class: 'xp-ctrl' });
    const lab = h('label', undefined, { for: 'xp-cmp-' + p.k }); lab.innerHTML = p.label;
    const val = h('span', undefined, { class: 'xp-val' }); val.style.color = p.colour;
    const inp = h('input', undefined, { type: 'range', id: 'xp-cmp-' + p.k, min: p.min, max: p.max, step: p.step, value: p.value });
    inp.style.accentColor = p.colour;
    const row = h('div', undefined, { class: 'xp-ctrl-row' }); row.append(lab, val);
    wrap.append(row, inp, h('span', p.hint, { class: 'xp-hint' }));
    bar.append(wrap);
    inputs[p.k] = { inp, val, p };
    inp.addEventListener('input', () => { v[p.k] = parseFloat(inp.value); update(); });
  }
  const apply = vals => {
    for (const p of spec.params) { const nv = vals[p.k] !== undefined ? vals[p.k] : p.value; if (p.k === 'view' && vals.view === undefined) continue; v[p.k] = nv; inputs[p.k].inp.value = nv; }
    update();
  };
  const pr = h('div', undefined, { class: 'xp-presets' });
  pr.append(h('span', 'Scenarios', { class: 'xp-presets-label' }));
  for (const preset of spec.presets) { const b = h('button', preset.label, { type: 'button', class: 'xp-pre' }); b.addEventListener('click', () => apply(preset.v)); pr.append(b); }
  const reset = h('button', 'Reset values', { type: 'button', class: 'text-button xp-reset' });
  reset.addEventListener('click', () => { v.view = 0; inputs.view.inp.value = 0; apply({}); });
  const ctl = h('div', undefined, { class: 'xp-cmp-bar' });
  ctl.append(bar, pr, reset);
  root.append(ctl);

  const grid = h('div', undefined, { class: 'xp-cmp-grid' });
  const panels = spec.models.map(m => {
    const box = h('section', undefined, { class: 'xp-cmp-panel xp-cmp-' + m.key });
    const top = h('div', undefined, { class: 'xp-cmp-top' });
    const t = h('h3', undefined); t.append(m.title + ' ', h('small', m.tag, { class: 'xp-cmp-tag' }));
    const eq = h('p', undefined, { class: 'xp-cmp-eq' }); eq.innerHTML = m.eq;
    top.append(t, eq);
    bindTerms(top, spec.terms[m.terms], [eq]);
    const live = h('p', undefined, { class: 'xp-live' });
    const svg = s('svg', { role: 'img', 'aria-label': m.title + ' diagram' });
    const cards = h('div', undefined, { class: 'xp-cards' });
    box.append(top, live, svg, cards); grid.append(box);
    return { key: m.key, live, svg, cards };
  });
  const pv = PROVENANCE[spec.id];
  if (pv) { const tag = h('div', undefined, { class: 'xp-prov xp-prov-' + pv[0] }); tag.append(h('strong', pv[1]), h('span', pv[2])); root.append(tag); }
  const verdict = h('div', undefined, { class: 'xp-verdict' });
  root.append(grid, verdict, h('p', spec.note, { class: 'xp-note' }));

  function update() {
    for (const p of spec.params) inputs[p.k].val.textContent = p.fmt ? p.fmt(v[p.k]) : minus(v[p.k]);
    for (const pn of panels) {
      const out = spec.panel(pn.key, pn.svg, v, mem);
      pn.live.innerHTML = out.live;
      pn.cards.replaceChildren();
      for (const [label, big, sub, lead] of out.cards) {
        const c = h('div', undefined, { class: 'xp-card' + (lead ? ' lead' : '') });
        const hh = h('h4'); hh.innerHTML = label;
        c.append(hh, h('strong', big), h('span', sub)); pn.cards.append(c);
      }
    }
    const vd = spec.verdict(v);
    verdict.replaceChildren(h('strong', vd[0] + '. '), h('span', vd[1]));
  }
  update();
  return update;
}

/* Each marked part of an equation opens a speech bubble explaining it. */
function bindTerms(head, terms, eqEls) {
  const bubble = h('div', undefined, { class: 'xp-bubble', role: 'dialog', 'aria-live': 'polite' });
  bubble.hidden = true;
  let open = null;
  const close = () => { bubble.hidden = true; if (open) open.setAttribute('aria-expanded', 'false'); open = null; };
  const show = term => {
    const [title, text] = terms[term.dataset.term] || [];
    if (!title) return;
    if (open === term) { close(); return; }
    close();
    bubble.replaceChildren(h('strong', title));
    const body = h('span'); body.innerHTML = text; bubble.append(body);
    bubble.hidden = false; open = term; term.setAttribute('aria-expanded', 'true');
    const hr = head.getBoundingClientRect(), tr = term.getBoundingClientRect();
    const left = Math.max(0, Math.min(tr.left - hr.left + tr.width / 2 - 150, hr.width - 300));
    bubble.style.left = left + 'px';
    bubble.style.top = (tr.bottom - hr.top + 10) + 'px';
    bubble.style.setProperty('--tail', Math.max(14, Math.min(286, tr.left - hr.left + tr.width / 2 - left)) + 'px');
  };
  for (const el of eqEls) for (const term of el.querySelectorAll('[data-term]')) {
    term.classList.add('xp-term');
    term.setAttribute('tabindex', '0');
    term.setAttribute('role', 'button');
    term.setAttribute('aria-expanded', 'false');
    term.title = 'Click for an explanation';
    term.addEventListener('click', e => { e.stopPropagation(); show(term); });
    term.addEventListener('keydown', e => { if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); show(term); } });
  }
  document.addEventListener('click', e => { if (open && !bubble.contains(e.target)) close(); });
  document.addEventListener('keydown', e => { if (e.key === 'Escape') close(); });
  head.append(bubble, h('p', 'Click any underlined part of the equation for an explanation.', { class: 'xp-termhint' }));
}

export function mountExplainer(root, id) {
  const spec = EXPLAINERS.find(e => e.id === id) || EXPLAINERS[0];
  if (spec.compare) return mountCompare(root, spec);
  const v = {};
  for (const p of spec.params) v[p.k] = p.value;
  const mem = spec.init();

  root.replaceChildren();
  const head = h('div', undefined, { class: 'xp-head' });
  head.append(h('p', spec.ref, { class: 'xp-ref' }));
  const eqMain = h('p', undefined, { class: 'xp-eqmain' });
  eqMain.innerHTML = spec.eq;
  head.append(eqMain);
  if (spec.terms) bindTerms(head, spec.terms, [eqMain]);
  head.append(h('p', spec.blurb, { class: 'xp-blurb' }));
  root.append(head);

  const realWrap = h('div', undefined, { class: 'xp-real' });
  realWrap.hidden = true;
  root.append(realWrap);
  const stage = h('div', undefined, { class: 'xp-stage' });
  const figure = h('div', undefined, { class: 'xp-figure' });
  const pv = PROVENANCE[spec.id];
  if (pv) {
    const tag = h('div', undefined, { class: 'xp-prov xp-prov-' + pv[0] });
    tag.append(h('strong', pv[1]), h('span', pv[2]));
    figure.append(tag);
  }
  const live = h('p', undefined, { class: 'xp-live' });
  const svg = s('svg', { role: 'img', 'aria-label': spec.title + ' diagram' });
  figure.append(live, svg);

  const controls = h('div', undefined, { class: 'xp-controls' });
  const inputs = {};
  for (const p of spec.params) {
    const wrap = h('div', undefined, { class: 'xp-ctrl' });
    const lab = h('label', undefined, { for: 'xp-' + spec.id + '-' + p.k });
    lab.innerHTML = p.label;
    const val = h('span', undefined, { class: 'xp-val' });
    val.style.color = p.colour;
    const inp = h('input', undefined, {
      type: 'range', id: 'xp-' + spec.id + '-' + p.k,
      min: p.min, max: p.max, step: p.step, value: p.value
    });
    inp.style.accentColor = p.colour;
    const row = h('div', undefined, { class: 'xp-ctrl-row' });
    row.append(lab, val);
    wrap.append(row, inp, h('span', p.hint, { class: 'xp-hint' }));
    controls.append(wrap);
    inputs[p.k] = { inp, val, p };
    inp.addEventListener('input', () => { v[p.k] = parseFloat(inp.value); update(); });
  }
  const apply = vals => {
    for (const p of spec.params) {
      const nv = vals[p.k] !== undefined ? vals[p.k] : p.value;
      v[p.k] = nv; inputs[p.k].inp.value = nv;
    }
    update();
  };
  if (spec.presets) {
    const pr = h('div', undefined, { class: 'xp-presets' });
    pr.append(h('span', 'Scenarios', { class: 'xp-presets-label' }));
    for (const preset of spec.presets) {
      const b = h('button', preset.label, { type: 'button', class: 'xp-pre' });
      b.addEventListener('click', () => apply(preset.v));
      pr.append(b);
    }
    controls.append(pr);
  }
  const reset = h('button', 'Reset values', { type: 'button', class: 'text-button xp-reset' });
  reset.addEventListener('click', () => apply({}));
  controls.append(reset);

  stage.append(figure, controls);
  const cards = h('div', undefined, { class: 'xp-cards' });
  const verdict = h('div', undefined, { class: 'xp-verdict' });
  root.append(stage, cards, verdict, h('p', spec.note, { class: 'xp-note' }));

  function update() {
    for (const p of spec.params) {
      const { val } = inputs[p.k];
      val.textContent = p.fmt ? p.fmt(v[p.k]) : minus(v[p.k]);
    }
    const out = spec.draw(svg, v, mem);
    live.innerHTML = out.eq;
    cards.replaceChildren();
    for (const [label, big, sub, lead] of out.cards) {
      const c = h('div', undefined, { class: 'xp-card' + (lead ? ' lead' : '') });
      const t = h('h4');
      t.innerHTML = label;
      c.append(t, h('strong', big), h('span', sub));
      cards.append(c);
    }
    verdict.replaceChildren(h('strong', out.verdict[0] + '. '), h('span', out.verdict[1]));
  }
  function buildReal() {
    if (!spec.real) return;
    const groups = spec.real();
    if (!groups || !groups.length) return;
    realWrap.replaceChildren();
    realWrap.append(h('label', 'Load a published estimate', { for: 'xp-real-' + spec.id }));
    const sel = h('select', undefined, { id: 'xp-real-' + spec.id });
    sel.append(h('option', 'illustrative values', { value: '' }));
    const flat = [];
    for (const g of groups) {
      const og = document.createElement('optgroup');
      og.label = g.group;
      for (const it of g.items) {
        og.append(h('option', it.label, { value: String(flat.length) }));
        flat.push(it);
      }
      sel.append(og);
    }
    const badge = h('span', '', { class: 'xp-badge' });
    sel.addEventListener('change', () => {
      if (sel.value === '') { badge.textContent = ''; apply({}); return; }
      badge.textContent = 'published';
      apply(flat[+sel.value].v);
    });
    realWrap.append(sel, badge);
    realWrap.hidden = false;
  }
  update();
  if (spec.usesData) loadExplainerData().then(() => { buildReal(); update(); });
  return update;
}
