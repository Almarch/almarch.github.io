/* Breeder's equation interactive widget.
   Hover sets the truncation point, click freezes the cycle and draws the
   offspring distribution. Standardized scale: mu = 0, sigma_P = 1. */
(function () {
  "use strict";

  const CFG = {
    width: 660, height: 340,
    margin: { top: 20, right: 20, bottom: 42, left: 20 },
    xMin: -3, xMax: 3,          // display and clamping domain
    yMax: 0.60,                 // fixed so the plot never rescales on click
    nPoints: 241,
    colorBase: "#E8734A",       // coral orange
    colorSel: "#B98FC4",        // light plum
    colorMuted: "#B9BCC0"
  };

  const SVG_NS = "http://www.w3.org/2000/svg";

  // ---- statistics -------------------------------------------------------

  function pdf(x, mu, sigma) {
    const z = (x - mu) / sigma;
    return Math.exp(-0.5 * z * z) / (sigma * Math.sqrt(2 * Math.PI));
  }

  /* Standard normal CDF, Abramowitz & Stegun 26.2.17 (|err| < 7.5e-8). */
  function cdf(x) {
    const b = [0.319381530, -0.356563782, 1.781477937, -1.821255978, 1.330274429];
    const sign = x < 0 ? -1 : 1;
    const ax = Math.abs(x);
    const t = 1 / (1 + 0.2316419 * ax);
    let poly = 0;
    for (let k = 4; k >= 0; k--) poly = (poly + b[k]) * t;
    const tail = pdf(ax, 0, 1) * poly;
    return 0.5 * (1 + sign * (1 - 2 * tail));
  }

  /* Selection outcome for a truncation point x on a standardized trait. */
  function selection(x, h2) {
    const p = 1 - cdf(x);                       // selected proportion
    const i = p > 1e-12 ? pdf(x, 0, 1) / p : 0; // selection intensity
    const S = i;                                // selection differential
    const R = h2 * i;                           // response to selection
    const k = i * (i - x);                      // variance reduction factor
    const varOff = 1 - k * h2 * h2 / 2;         // phenotypic, standardized
    const varAoff = h2 * (1 - k * h2 / 2);      // additive, standardized
    return {
      p: p, i: i, S: S, R: R,
      sigmaOff: Math.sqrt(varOff),
      h2Off: varAoff / varOff
    };
  }

  // ---- svg helpers ------------------------------------------------------

  function el(tag, attrs) {
    const node = document.createElementNS(SVG_NS, tag);
    for (const key in attrs) node.setAttribute(key, attrs[key]);
    return node;
  }

  function makeScales() {
    const m = CFG.margin;
    const w = CFG.width - m.left - m.right;
    const h = CFG.height - m.top - m.bottom;
    return {
      x: v => m.left + (v - CFG.xMin) / (CFG.xMax - CFG.xMin) * w,
      y: v => m.top + h - (v / CFG.yMax) * h,
      xInv: px => CFG.xMin + (px - m.left) / w * (CFG.xMax - CFG.xMin),
      baseline: m.top + h
    };
  }

  /* Filled area under a normal density between lo and hi. */
  function areaPath(sc, lo, hi, mu, sigma) {
    if (hi <= lo) return "";
    const n = Math.max(2, Math.round(CFG.nPoints * (hi - lo) / (CFG.xMax - CFG.xMin)));
    const pts = [];
    for (let k = 0; k <= n; k++) {
      const x = lo + (hi - lo) * k / n;
      pts.push(sc.x(x).toFixed(2) + "," + sc.y(pdf(x, mu, sigma)).toFixed(2));
    }
    return "M" + sc.x(lo).toFixed(2) + "," + sc.baseline.toFixed(2) +
           " L" + pts.join(" L") +
           " L" + sc.x(hi).toFixed(2) + "," + sc.baseline.toFixed(2) + " Z";
  }

  /* Open curve (no baseline closure) for a normal density. */
  function curvePath(sc, mu, sigma) {
    const pts = [];
    for (let k = 0; k <= CFG.nPoints; k++) {
      const x = CFG.xMin + (CFG.xMax - CFG.xMin) * k / CFG.nPoints;
      pts.push(sc.x(x).toFixed(2) + "," + sc.y(pdf(x, mu, sigma)).toFixed(2));
    }
    return "M" + pts.join(" L");
  }

  // ---- widget -----------------------------------------------------------

  function build(root) {
    const sc = makeScales();
    const state = { x: 0.8, h2: 0.4, frozen: false };

    root.classList.add("bw-root");
    root.innerHTML =
      '<div class="bw-layout">' +
        '<div class="bw-plot"></div>' +
        '<div class="bw-panel">' +
          '<label class="bw-label">Heritability <span class="bw-h2">0.40</span></label>' +
          '<input class="bw-slider" type="range" min="0" max="1" step="0.01" value="0.4">' +
          '<dl class="bw-readout">' +
            '<dt>Selected</dt><dd class="bw-p"></dd>' +
            '<dt>Intensity <em>i</em></dt><dd class="bw-i"></dd>' +
            '<dt>Response <em>R</em></dt><dd class="bw-R"></dd>' +
            '<dt>Offspring &sigma;<sub>P</sub></dt><dd class="bw-sd"></dd>' +
            '<dt>Offspring <em>h&sup2;</em></dt><dd class="bw-h2off"></dd>' +
          '</dl>' +
        '</div>' +
      '</div>';

    const svg = el("svg", {
      viewBox: "0 0 " + CFG.width + " " + CFG.height,
      class: "bw-svg", preserveAspectRatio: "xMidYMid meet"
    });
    root.querySelector(".bw-plot").appendChild(svg);

    const clip = el("clipPath", { id: "bw-clip" });
    clip.appendChild(el("rect", {
      x: CFG.margin.left, y: CFG.margin.top,
      width: CFG.width - CFG.margin.left - CFG.margin.right,
      height: CFG.height - CFG.margin.top - CFG.margin.bottom
    }));
    const defs = el("defs", {});
    defs.appendChild(clip);
    svg.appendChild(defs);

    const nodes = {
      areaOut:  el("path",  { class: "bw-area-out" }),
      areaSel:  el("path",  { class: "bw-area-sel" }),
      curve:    el("path",  { class: "bw-curve" }),
      offArea:  el("path",  { class: "bw-off-area", "clip-path": "url(#bw-clip)" }),
      offCurve: el("path",  { class: "bw-off-curve", "clip-path": "url(#bw-clip)" }),
      axis:     el("line",  { class: "bw-axis", x1: CFG.margin.left, x2: CFG.width - CFG.margin.right,
                              y1: sc.baseline, y2: sc.baseline }),
      vBase:    el("line",  { class: "bw-vbase", x1: sc.x(0), x2: sc.x(0),
                              y1: CFG.margin.top, y2: sc.baseline }),
      vResp:    el("line",  { class: "bw-vresp" }),
      vTrunc:   el("line",  { class: "bw-vtrunc" }),
      labelR:   el("text",  { class: "bw-label-r", "text-anchor": "middle" })
    };

    [nodes.areaOut, nodes.areaSel, nodes.offArea, nodes.curve, nodes.offCurve,
     nodes.axis, nodes.vBase, nodes.vResp, nodes.vTrunc, nodes.labelR]
      .forEach(n => svg.appendChild(n));

    for (let t = CFG.xMin; t <= CFG.xMax + 1e-9; t += 1) {
      const tick = el("text", {
        class: "bw-tick", x: sc.x(t), y: sc.baseline + 20, "text-anchor": "middle"
      });
      tick.textContent = t.toFixed(0);
      svg.appendChild(tick);
    }
    const xTitle = el("text", {
      class: "bw-axis-title", x: (CFG.margin.left + CFG.width - CFG.margin.right) / 2,
      y: CFG.height - 6, "text-anchor": "middle"
    });
    xTitle.textContent = "Phenotype (standardized)";
    svg.appendChild(xTitle);

    const ui = {
      h2:     root.querySelector(".bw-h2"),
      slider: root.querySelector(".bw-slider"),
      p:      root.querySelector(".bw-p"),
      i:      root.querySelector(".bw-i"),
      R:      root.querySelector(".bw-R"),
      sd:     root.querySelector(".bw-sd"),
      h2off: root.querySelector(".bw-h2off"),
    };

    function render() {
      const res = selection(state.x, state.h2);
      const px = sc.x(state.x);

      nodes.areaOut.setAttribute("d", areaPath(sc, CFG.xMin, state.x, 0, 1));
      nodes.areaSel.setAttribute("d", areaPath(sc, state.x, CFG.xMax, 0, 1));
      nodes.curve.setAttribute("d", curvePath(sc, 0, 1));

      nodes.vTrunc.setAttribute("x1", px);
      nodes.vTrunc.setAttribute("x2", px);
      nodes.vTrunc.setAttribute("y1", CFG.margin.top);
      nodes.vTrunc.setAttribute("y2", sc.baseline);

      if (state.frozen) {
        nodes.offArea.setAttribute("d", areaPath(sc, CFG.xMin, CFG.xMax, res.R, res.sigmaOff));
        nodes.offCurve.setAttribute("d", curvePath(sc, res.R, res.sigmaOff));
        const pr = sc.x(Math.min(res.R, CFG.xMax));
        nodes.vResp.setAttribute("x1", pr);
        nodes.vResp.setAttribute("x2", pr);
        nodes.vResp.setAttribute("y1", CFG.margin.top);
        nodes.vResp.setAttribute("y2", sc.baseline);
        nodes.labelR.setAttribute("x", (sc.x(0) + pr) / 2);
        nodes.labelR.setAttribute("y", CFG.margin.top + 12);
        nodes.labelR.textContent = "R = " + res.R.toFixed(2);
      } else {
        nodes.offArea.setAttribute("d", "");
        nodes.offCurve.setAttribute("d", "");
        nodes.vResp.setAttribute("x1", -10);
        nodes.vResp.setAttribute("x2", -10);
        nodes.labelR.textContent = "";
      }
      root.classList.toggle("bw-frozen", state.frozen);

      ui.h2.textContent = state.h2.toFixed(2);
      ui.p.textContent = (100 * res.p).toFixed(1) + " %";
      ui.i.textContent = res.i.toFixed(2);
      ui.R.textContent = res.R.toFixed(2);
      ui.sd.textContent = res.sigmaOff.toFixed(3);
      ui.h2off.textContent = res.h2Off.toFixed(3);
    }

    /* Pointer events cover mouse, pen and touch alike. */
    function setTruncation(ev) {
      if (state.frozen) return;
      const box = svg.getBoundingClientRect();
      const px = (ev.clientX - box.left) / box.width * CFG.width;
      const raw = sc.xInv(px);
      state.x = Math.min(CFG.xMax, Math.max(CFG.xMin, raw));
      render();
    }

    svg.addEventListener("pointermove", setTruncation);
    svg.addEventListener("pointerdown", setTruncation);

    svg.addEventListener("click", function () {
      state.frozen = !state.frozen;
      render();
    });

    ui.slider.addEventListener("input", function () {
      state.h2 = parseFloat(this.value);
      render();
    });
    render();
  }

  function init() {
    const root = document.getElementById("breeder-widget");
    if (root) build(root);
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", init);
  } else {
    init();
  }
})();

/* Selection index widget (Smith-Hazel with desired gains).
   Three standardized traits, n breeding values on the BLUP scale.
   The correlation quadrant defines G; the orthonormal frame stays fixed
   and the cloud deforms. Drag the cloud to rotate, drag arrow heads to
   change the genetic correlations. */
(function () {
  "use strict";

  const CFG = {
    n: 300,
    cloud: { width: 620, height: 400 },
    quad: { size: 150, pad: 22 },
    axisLen: 3.2,                       // 3D axis length, trait units
    scale: 62,                          // px per trait unit at zero depth
    persp: 0.09,                        // depth foreshortening strength
    traitColors: ["#E8734A", "#B98FC4", "#2E6B45"],  // coral, plum, forest
    traitNames: ["Trait 1", "Trait 2", "Trait 3"],
    ptBase: "#1A1A1A",
    ptSel: "#D62828",
    seed: 20240517
  };

  const SVG_NS = "http://www.w3.org/2000/svg";

  // ---- random numbers ---------------------------------------------------

  /* Mulberry32: deterministic PRNG so every reader sees the same cloud. */
  function makeRng(seed) {
    let a = seed >>> 0;
    return function () {
      a = (a + 0x6D2B79F5) >>> 0;
      let t = Math.imul(a ^ (a >>> 15), 1 | a);
      t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }

  /* n x 3 matrix of independent standard normal deviates (Box-Muller). */
  function sampleZ(n, rng) {
    const out = [];
    for (let r = 0; r < n; r++) {
      const row = [];
      for (let c = 0; c < 3; c++) {
        const u = Math.max(rng(), 1e-12);
        row.push(Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * rng()));
      }
      out.push(row);
    }
    return out;
  }

  // ---- linear algebra ---------------------------------------------------

  /* Unit vectors from quadrant coordinates; the third component is the
     residual needed to reach norm one, taken positive by convention. */
  function arrowsToV(arrows) {
    return arrows.map(function (a) {
      const z2 = Math.max(0, 1 - a[0] * a[0] - a[1] * a[1]);
      return [a[0], a[1], Math.sqrt(z2)];
    });
  }

  /* G = V'V. Unit vectors give a unit diagonal, so G is a correlation
     matrix and is positive semi-definite by construction. */
  function gramMatrix(V) {
    const G = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
    for (let j = 0; j < 3; j++) {
      for (let k = 0; k < 3; k++) {
        let s = 0;
        for (let m = 0; m < 3; m++) s += V[j][m] * V[k][m];
        G[j][k] = s;
      }
    }
    return G;
  }

  function det3(M) {
    return M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1])
         - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
         + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]);
  }

  /* Index weights b = G^-1 d, by cofactor expansion. Returns null when the
     traits are collinear and the index is undefined. */
  function solveIndex(G, d) {
    const D = det3(G);
    if (Math.abs(D) < 1e-6) return null;
    const C = [
      [G[1][1] * G[2][2] - G[1][2] * G[2][1],
       G[0][2] * G[2][1] - G[0][1] * G[2][2],
       G[0][1] * G[1][2] - G[0][2] * G[1][1]],
      [G[1][2] * G[2][0] - G[1][0] * G[2][2],
       G[0][0] * G[2][2] - G[0][2] * G[2][0],
       G[0][2] * G[1][0] - G[0][0] * G[1][2]],
      [G[1][0] * G[2][1] - G[1][1] * G[2][0],
       G[0][1] * G[2][0] - G[0][0] * G[2][1],
       G[0][0] * G[1][1] - G[0][1] * G[1][0]]
    ];
    const b = [0, 0, 0];
    for (let r = 0; r < 3; r++) {
      let s = 0;
      for (let c = 0; c < 3; c++) s += C[r][c] * d[c];
      b[r] = s / D;
    }
    return b;
  }

  /* Rotate a point by yaw (around the vertical) then pitch. */
  function rotate(p, yaw, pitch) {
    const cy = Math.cos(yaw), sy = Math.sin(yaw);
    const cp = Math.cos(pitch), sp = Math.sin(pitch);
    const x1 =  p[0] * cy + p[2] * sy;
    const z1 = -p[0] * sy + p[2] * cy;
    const y2 =  p[1] * cp - z1 * sp;
    const z2 =  p[1] * sp + z1 * cp;
    return [x1, y2, z2];
  }

  // ---- svg helper -------------------------------------------------------

  function el(tag, attrs) {
    const node = document.createElementNS(SVG_NS, tag);
    for (const key in attrs) node.setAttribute(key, attrs[key]);
    return node;
  }

  // ---- widget -----------------------------------------------------------

  function build(root) {
    const rng = makeRng(CFG.seed);
    const state = {
      Z: sampleZ(CFG.n, rng),
      arrows: [[0, 0.99], [0.99, 0], [0, 0]],   // orthogonal traits at start
      d: [1, 1, 1],
      nSel: Math.round(CFG.n * 0.2),
      yaw: 0.65, pitch: -0.30,
      dragCloud: false, dragArrow: -1,
      lastX: 0, lastY: 0
    };

    root.classList.add("sx-root");
    root.innerHTML =
      '<div class="sx-layout">' +
        '<div class="sx-plot"></div>' +
        '<div class="sx-panel">' +
          '<label class="sx-label sx-label-first">Genetic correlations</label>' +
          '<div class="sx-quad"></div>' +
          '<div class="sx-gains"></div>' +
          '<label class="sx-label">Selected ' +
            '<span class="sx-nsel"></span> / ' + CFG.n + '</label>' +
          '<input class="sx-slider" type="range" min="1" max="' + CFG.n +
            '" step="1" value="' + state.nSel + '">' +
          '<table class="sx-table"><thead><tr>' +
            '<th></th><th><em>b</em></th><th><em>d</em></th><th>&Delta;</th>' +
          '</tr></thead><tbody></tbody></table>' +
          '<p class="sx-warn"></p>' +
        '</div>' +
      '</div>';

    // --- desired gain steppers
    const gainsBox = root.querySelector(".sx-gains");
    for (let j = 0; j < 3; j++) {
      const row = document.createElement("div");
      row.className = "sx-stepper";
      row.innerHTML =
        '<span class="sx-swatch" style="background:' + CFG.traitColors[j] + '"></span>' +
        '<span class="sx-sname">' + CFG.traitNames[j] + '</span>' +
        '<button type="button" data-j="' + j + '" data-dir="-1">&minus;</button>' +
        '<span class="sx-dval" data-j="' + j + '"></span>' +
        '<button type="button" data-j="' + j + '" data-dir="1">+</button>';
      gainsBox.appendChild(row);
    }

    // --- main cloud
    const cw = CFG.cloud.width, ch = CFG.cloud.height;
    const svgC = el("svg", {
      viewBox: "0 0 " + cw + " " + ch,
      class: "sx-svg-cloud", preserveAspectRatio: "xMidYMid meet"
    });
    root.querySelector(".sx-plot").appendChild(svgC);

    const gAxesBack = el("g", {});
    const gPts = el("g", {});
    const gAxesFront = el("g", {});
    svgC.appendChild(gAxesBack);
    svgC.appendChild(gPts);
    svgC.appendChild(gAxesFront);

    const axNodes = [];
    for (let j = 0; j < 3; j++) {
      const line = el("line", { class: "sx-ax3d", stroke: CFG.traitColors[j] });
      const lab = el("text", {
        class: "sx-axlab", fill: CFG.traitColors[j],
        "text-anchor": "middle", "dominant-baseline": "middle"
      });
      lab.textContent = "T" + (j + 1);
      axNodes.push({ line: line, lab: lab });
    }

    const ptNodes = [];
    for (let r = 0; r < CFG.n; r++) {
      const c = el("circle", { class: "sx-pt" });
      gPts.appendChild(c);
      ptNodes.push(c);
    }

    // --- correlation quadrant
    const qs = CFG.quad.size, qp = CFG.quad.pad;
    const qSpan = qs - qp;
    const svgQ = el("svg", {
      viewBox: "0 0 " + (qs + 10) + " " + (qs + 10), class: "sx-svg-quad"
    });
    root.querySelector(".sx-quad").appendChild(svgQ);

    const qx = v => qp + v * qSpan;
    const qy = v => qs - v * qSpan;

    /* sweep-flag 0: with y pointing down, this draws the convex arc that
       bulges away from the origin. */
    svgQ.appendChild(el("path", {
      class: "sx-arc",
      d: "M " + qx(1) + " " + qy(0) +
         " A " + qSpan + " " + qSpan + " 0 0 0 " + qx(0) + " " + qy(1)
    }));
    svgQ.appendChild(el("line", {
      class: "sx-qax", x1: qx(0), y1: qy(0), x2: qx(1.12), y2: qy(0)
    }));
    svgQ.appendChild(el("line", {
      class: "sx-qax", x1: qx(0), y1: qy(0), x2: qx(0), y2: qy(1.12)
    }));

    const arrowNodes = [];
    for (let j = 0; j < 3; j++) {
      const line = el("line", { class: "sx-arrow", stroke: CFG.traitColors[j] });
      const head = el("circle", {
        class: "sx-head", r: 5.5, fill: CFG.traitColors[j], "data-j": j
      });
      svgQ.appendChild(line);
      svgQ.appendChild(head);
      arrowNodes.push({ line: line, head: head });
    }

    const ui = {
      slider: root.querySelector(".sx-slider"),
      nsel: root.querySelector(".sx-nsel"),
      tbody: root.querySelector(".sx-table tbody"),
      warn: root.querySelector(".sx-warn"),
      dvals: Array.prototype.slice.call(root.querySelectorAll(".sx-dval"))
    };

    for (let j = 0; j < 3; j++) {
      const tr = document.createElement("tr");
      tr.innerHTML =
        '<td><span class="sx-swatch" style="background:' + CFG.traitColors[j] +
        '"></span>T' + (j + 1) + '</td>' +
        '<td class="sx-b"></td><td class="sx-d"></td><td class="sx-delta"></td>';
      ui.tbody.appendChild(tr);
    }
    const cellsB = root.querySelectorAll(".sx-b");
    const cellsD = root.querySelectorAll(".sx-d");
    const cellsDelta = root.querySelectorAll(".sx-delta");

    // ---- render ---------------------------------------------------------

    function project(p) {
      const r = rotate(p, state.yaw, state.pitch);
      const f = 1 / (1 + CFG.persp * (1 - r[2] / CFG.axisLen));
      return {
        px: cw / 2 + r[0] * CFG.scale * f,
        py: ch / 2 - r[1] * CFG.scale * f,
        z: r[2], f: f
      };
    }

    function render() {
      const V = arrowsToV(state.arrows);
      const G = gramMatrix(V);
      const b = solveIndex(G, state.d);

      /* Trait values: project each latent point onto the trait unit vector.
         Var(T_j) = 1 and Cov(T_j, T_k) = G_jk, so the frame stays orthonormal
         while the cloud takes the shape imposed by G. */
      const T = state.Z.map(function (z) {
        const t = [0, 0, 0];
        for (let j = 0; j < 3; j++) {
          t[j] = z[0] * V[j][0] + z[1] * V[j][1] + z[2] * V[j][2];
        }
        return t;
      });

      const sel = new Array(CFG.n).fill(false);
      let delta = [0, 0, 0];

      if (b) {
        const s = T.map(t => t[0] * b[0] + t[1] * b[1] + t[2] * b[2]);
        const order = s.map((v, k) => k).sort((p, q) => s[q] - s[p]);
        for (let k = 0; k < state.nSel; k++) {
          sel[order[k]] = true;
          for (let j = 0; j < 3; j++) delta[j] += T[order[k]][j];
        }
        delta = delta.map(v => v / state.nSel);
        ui.warn.textContent = "";
      } else {
        ui.warn.textContent = "Traits are collinear: G is singular.";
      }

      // --- fixed orthonormal frame, drawn behind or in front of the cloud
      const axes = [];
      for (let j = 0; j < 3; j++) {
        const dir = [0, 0, 0];
        dir[j] = CFG.axisLen;
        const tip = project(dir);
        const org = project([0, 0, 0]);
        axNodes[j].line.setAttribute("x1", org.px.toFixed(1));
        axNodes[j].line.setAttribute("y1", org.py.toFixed(1));
        axNodes[j].line.setAttribute("x2", tip.px.toFixed(1));
        axNodes[j].line.setAttribute("y2", tip.py.toFixed(1));
        axNodes[j].lab.setAttribute("x", (org.px + (tip.px - org.px) * 1.09).toFixed(1));
        axNodes[j].lab.setAttribute("y", (org.py + (tip.py - org.py) * 1.09).toFixed(1));
        axes.push({ j: j, z: tip.z });
      }
      axes.forEach(function (a) {
        const target = a.z < 0 ? gAxesBack : gAxesFront;
        target.appendChild(axNodes[a.j].line);
        target.appendChild(axNodes[a.j].lab);
      });

      // --- cloud, painted back to front
      const proj = T.map(function (t, k) {
        const p = project(t);
        p.k = k;
        return p;
      });
      proj.sort((a, z) => a.z - z.z);
      proj.forEach(function (p, idx) {
        const node = ptNodes[idx];
        const on = sel[p.k];
        node.setAttribute("cx", p.px.toFixed(1));
        node.setAttribute("cy", p.py.toFixed(1));
        node.setAttribute("r", ((on ? 3.0 : 2.4) * p.f).toFixed(2));
        node.setAttribute("fill", on ? CFG.ptSel : CFG.ptBase);
        node.setAttribute("opacity", on ? 1 : 0.5);
      });

      // --- quadrant arrows
      for (let j = 0; j < 3; j++) {
        const a = state.arrows[j];
        arrowNodes[j].line.setAttribute("x1", qx(0));
        arrowNodes[j].line.setAttribute("y1", qy(0));
        arrowNodes[j].line.setAttribute("x2", qx(a[0]));
        arrowNodes[j].line.setAttribute("y2", qy(a[1]));
        arrowNodes[j].head.setAttribute("cx", qx(a[0]));
        arrowNodes[j].head.setAttribute("cy", qy(a[1]));
      }

      // --- readouts
      ui.nsel.textContent = state.nSel;
      const bn = b ? Math.sqrt(b[0] * b[0] + b[1] * b[1] + b[2] * b[2]) : 1;
      for (let j = 0; j < 3; j++) {
        ui.dvals[j].textContent = state.d[j].toFixed(1);
        cellsB[j].textContent = b ? (b[j] / bn).toFixed(2) : "\u2014";
        cellsD[j].textContent = state.d[j].toFixed(1);
        cellsDelta[j].textContent = b ? delta[j].toFixed(2) : "\u2014";
      }
    }

    // ---- interaction ----------------------------------------------------

    function svgPoint(svg, ev) {
      const box = svg.getBoundingClientRect();
      const vb = svg.viewBox.baseVal;
      return {
        x: (ev.clientX - box.left) / box.width * vb.width,
        y: (ev.clientY - box.top) / box.height * vb.height
      };
    }

    /* Pointer events cover mouse, pen and touch alike. */
    svgC.addEventListener("pointerdown", function (ev) {
      state.dragCloud = true;
      state.lastX = ev.clientX;
      state.lastY = ev.clientY;
      ev.preventDefault();
    });

    svgQ.addEventListener("pointerdown", function (ev) {
      const attr = ev.target.getAttribute && ev.target.getAttribute("data-j");
      if (attr !== null && attr !== undefined) {
        state.dragArrow = parseInt(attr, 10);
        ev.preventDefault();
      }
    });

    window.addEventListener("pointermove", function (ev) {
      if (state.dragCloud) {
        state.yaw += (ev.clientX - state.lastX) * 0.009;
        state.pitch += (ev.clientY - state.lastY) * 0.009;
        state.pitch = Math.max(-1.4, Math.min(1.4, state.pitch));
        state.lastX = ev.clientX;
        state.lastY = ev.clientY;
        render();
      } else if (state.dragArrow >= 0) {
        const p = svgPoint(svgQ, ev);
        let ux = Math.max(0, (p.x - qp) / qSpan);
        let uy = Math.max(0, (qs - p.y) / qSpan);
        const nrm = Math.sqrt(ux * ux + uy * uy);
        if (nrm > 1) { ux /= nrm; uy /= nrm; }
        state.arrows[state.dragArrow] = [ux, uy];
        render();
      }
    });

    function endDrag() {
      state.dragCloud = false;
      state.dragArrow = -1;
    }
    window.addEventListener("pointerup", endDrag);
    window.addEventListener("pointercancel", endDrag);

    ui.slider.addEventListener("input", function () {
      state.nSel = parseInt(this.value, 10);
      render();
    });

    gainsBox.addEventListener("click", function (ev) {
      const btn = ev.target.closest("button");
      if (!btn) return;
      const j = parseInt(btn.getAttribute("data-j"), 10);
      const dir = parseInt(btn.getAttribute("data-dir"), 10);
      const next = Math.round((state.d[j] + dir * 0.1) * 10) / 10;
      state.d[j] = Math.max(-3, Math.min(3, next));
      render();
    });

    render();
  }

  function init() {
    const root = document.getElementById("index-widget");
    if (root) build(root);
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", init);
  } else {
    init();
  }
})();

/* Legendre polynomial widget.
   Draws L_0 up to L_order on [-1, 1] using Bonnet's recurrence, with an
   order stepper on the side panel. Orthonormal normalisation:
   L_n = sqrt((2n + 1) / 2) * P_n, as used in random regression models. */
(function () {
  "use strict";

  const CFG = {
    width: 620, height: 400,
    margin: { top: 18, right: 18, bottom: 46, left: 46 },
    maxOrder: 10,
    nPoints: 401,
    yMax: 3.4,                    // covers L_10(1) = sqrt(21/2) = 3.24
    /* Ordered rainbow: hue runs from red to violet so the colour encodes
       the degree. Hand-tuned lightness keeps the yellows readable. */
    colors: [
      "#D32F2F", "#E65100", "#EF9A00", "#C2A500", "#7CB342",
      "#2E7D32", "#00897B", "#0288D1", "#3949AB", "#6A3AB2", "#9C27B0"
    ]
  };

  const SVG_NS = "http://www.w3.org/2000/svg";

  // ---- mathematics ------------------------------------------------------

  /* Normalized Legendre polynomials, L_j = sqrt((2j+1)/2) * P_j, obtained
     from Bonnet's recurrence on the standard P_j:
     (n+1) P_{n+1}(x) = (2n+1) x P_n(x) - n P_{n-1}(x).
     This is the orthonormal basis used in random regression models
     (Kirkpatrick et al. 1990): the integral of L_j squared over [-1, 1]
     equals one. Returns the vector [L_0(x), ..., L_maxOrder(x)]. */
  function legendreAt(x, maxOrder) {
    const P = [1];
    if (maxOrder >= 1) P.push(x);
    for (let n = 1; n < maxOrder; n++) {
      P.push(((2 * n + 1) * x * P[n] - n * P[n - 1]) / (n + 1));
    }
    return P.map((v, n) => v * Math.sqrt((2 * n + 1) / 2));
  }

  // ---- svg helper -------------------------------------------------------

  function el(tag, attrs) {
    const node = document.createElementNS(SVG_NS, tag);
    for (const key in attrs) node.setAttribute(key, attrs[key]);
    return node;
  }

  // ---- widget -----------------------------------------------------------

  function build(root) {
    const state = { order: 1 };

    const m = CFG.margin;
    const plotW = CFG.width - m.left - m.right;
    const plotH = CFG.height - m.top - m.bottom;
    const sx = v => m.left + (v + 1) / 2 * plotW;
    const sy = v => m.top + plotH / 2 - v / CFG.yMax * (plotH / 2) * 0.96;

    /* Precompute every curve once; the stepper only toggles visibility. */
    const paths = [];
    for (let n = 0; n <= CFG.maxOrder; n++) paths.push([]);
    for (let k = 0; k < CFG.nPoints; k++) {
      const x = -1 + 2 * k / (CFG.nPoints - 1);
      const vals = legendreAt(x, CFG.maxOrder);
      for (let n = 0; n <= CFG.maxOrder; n++) {
        paths[n].push(sx(x).toFixed(2) + "," + sy(vals[n]).toFixed(2));
      }
    }

    root.classList.add("lp-root");
    root.innerHTML =
      '<div class="lp-layout">' +
        '<div class="lp-plot"></div>' +
        '<div class="lp-panel">' +
          '<label class="lp-label lp-label-first">Order</label>' +
          '<div class="lp-stepper">' +
            '<button type="button" data-dir="-1">&minus;</button>' +
            '<span class="lp-oval"></span>' +
            '<button type="button" data-dir="1">+</button>' +
          '</div>' +
          '<ul class="lp-legend"></ul>' +
        '</div>' +
      '</div>';

    const svg = el("svg", {
      viewBox: "0 0 " + CFG.width + " " + CFG.height,
      class: "lp-svg", preserveAspectRatio: "xMidYMid meet"
    });
    root.querySelector(".lp-plot").appendChild(svg);

    // --- axes
    svg.appendChild(el("line", {
      class: "lp-axis", x1: m.left, y1: sy(0), x2: m.left + plotW, y2: sy(0)
    }));
    svg.appendChild(el("line", {
      class: "lp-axis", x1: sx(0), y1: m.top, x2: sx(0), y2: m.top + plotH
    }));

    [-1, -0.5, 0.5, 1].forEach(function (t) {
      const lab = el("text", {
        class: "lp-tick", x: sx(t), y: sy(0) + 17, "text-anchor": "middle"
      });
      lab.textContent = t.toString();
      svg.appendChild(lab);
    });
    [-3, -2, -1, 1, 2, 3].forEach(function (t) {
      const lab = el("text", {
        class: "lp-tick", x: sx(0) - 8, y: sy(t) + 4, "text-anchor": "end"
      });
      lab.textContent = t.toString();
      svg.appendChild(lab);
    });

    const xTitle = el("text", {
      class: "lp-axis-title", x: m.left + plotW / 2,
      y: CFG.height - 8, "text-anchor": "middle"
    });
    xTitle.textContent = "x";
    svg.appendChild(xTitle);

    const yTitle = el("text", {
      class: "lp-axis-title", x: 14, y: m.top + plotH / 2,
      "text-anchor": "middle",
      transform: "rotate(-90 14 " + (m.top + plotH / 2) + ")"
    });
    const yt1 = el("tspan", {});
    yt1.textContent = "L";
    const yt2 = el("tspan", { "baseline-shift": "sub", "font-size": "9" });
    yt2.textContent = "m";
    const yt3 = el("tspan", { "baseline-shift": "0" });
    yt3.textContent = "(x)";
    yTitle.appendChild(yt1);
    yTitle.appendChild(yt2);
    yTitle.appendChild(yt3);
    svg.appendChild(yTitle);

    // --- curves
    const curveNodes = [];
    for (let n = 0; n <= CFG.maxOrder; n++) {
      const p = el("path", {
        class: "lp-curve", stroke: CFG.colors[n],
        d: "M" + paths[n].join(" L")
      });
      svg.appendChild(p);
      curveNodes.push(p);
    }

    // --- legend
    const legend = root.querySelector(".lp-legend");
    const legendNodes = [];
    for (let n = 0; n <= CFG.maxOrder; n++) {
      const li = document.createElement("li");
      li.innerHTML =
        '<span class="lp-swatch" style="background:' + CFG.colors[n] + '"></span>' +
        '<em>L</em><sub>' + n + '</sub>(<em>x</em>)';
      legend.appendChild(li);
      legendNodes.push(li);
    }

    const ui = {
      oval: root.querySelector(".lp-oval"),
      stepper: root.querySelector(".lp-stepper")
    };

    function render() {
      ui.oval.textContent = state.order;
      for (let n = 0; n <= CFG.maxOrder; n++) {
        const on = n <= state.order;
        curveNodes[n].setAttribute("visibility", on ? "visible" : "hidden");
        legendNodes[n].style.display = on ? "" : "none";
      }
      const btns = ui.stepper.querySelectorAll("button");
      btns[0].disabled = state.order <= 0;
      btns[1].disabled = state.order >= CFG.maxOrder;
    }

    ui.stepper.addEventListener("click", function (ev) {
      const btn = ev.target.closest("button");
      if (!btn) return;
      const dir = parseInt(btn.getAttribute("data-dir"), 10);
      state.order = Math.max(0, Math.min(CFG.maxOrder, state.order + dir));
      render();
    });

    render();
  }

  function init() {
    const root = document.getElementById("legendre-widget");
    if (root) build(root);
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", init);
  } else {
    init();
  }
})();