---
marp: true
size: 16:9
footer: |
  <img src="assets/otar-logo.png" alt="" /><span>Open Targets</span>
style: |
  :root {
    --ot-navy: #163A5F;
    --ot-blue: #2C77B5;
    --ot-blue-mid: #4A93CB;
    --ot-triangle: #6FA1CE;
    --ot-gray: #58595B;
    --ot-text: #141414;
  }

  section {
    font-family: Arial, Helvetica, sans-serif;
    color: var(--ot-text);
    background: #ffffff;
    padding: 56px 64px;
    position: relative;
    overflow: hidden;
  }

  section h1 {
    font-size: 36px;
    font-weight: 700;
    color: var(--ot-text);
    margin: 0 0 0.3em;
  }

  section h2 {
    font-size: 22px;
    font-weight: 400;
    color: var(--ot-blue);
    margin: 0;
  }

  section::before {
    content: "";
    position: absolute;
    left: 0;
    bottom: 0;
    width: 30%;
    height: 44%;
    background: var(--ot-triangle);
    clip-path: polygon(0 100%, 0 28%, 68% 100%);
    z-index: 0;
  }

  section footer {
    position: absolute;
    right: 36px;
    bottom: 28px;
    left: auto;
    display: flex;
    align-items: center;
    gap: 8px;
    font-size: 13px;
    color: var(--ot-gray);
    z-index: 1;
  }

  section footer img {
    height: 26px;
    width: auto;
    display: block;
  }

  section > * {
    position: relative;
    z-index: 2;
  }

  /* ---- Status slide ---- */
  section.status {
    padding: 48px 56px 64px;
    align-content: start !important;
  }

  section.status h1 {
    font-size: 30px;
    margin-bottom: 2px;
  }

  section.status h2 {
    font-size: 17px;
    margin-bottom: 16px;
  }

  section.status table {
    width: 100%;
    border-collapse: collapse;
    font-size: 17px;
    line-height: 1.3;
    background: #ffffff;
  }

  section.status th {
    text-align: left;
    font-size: 12px;
    text-transform: uppercase;
    letter-spacing: 0.04em;
    color: var(--ot-gray);
    background: #ffffff;
    border: none;
    border-bottom: 2px solid var(--ot-navy);
    padding: 4px 10px;
  }

  section.status td {
    vertical-align: top;
    background: #ffffff;
    border: none;
    border-bottom: 1px solid #dde3ea;
    padding: 11px 10px;
  }

  section.status td:first-child {
    white-space: nowrap;
  }

  section.status td strong {
    color: var(--ot-navy);
  }

  section.status .src {
    display: block;
    font-size: 12px;
    color: var(--ot-gray);
  }

  section.status .pill {
    display: inline-block;
    padding: 2px 9px;
    border-radius: 10px;
    font-size: 12px;
    font-weight: 700;
    white-space: nowrap;
  }

  .pill.done    { background: #DDF1E4; color: #1D6B3A; }
  .pill.wait    { background: #E3EEF8; color: var(--ot-navy); }
  .pill.blocked { background: #FCEBD5; color: #8A4B00; }
  .pill.scope   { background: #ECEAF4; color: #4A3F7A; }

  section.status td.num {
    white-space: nowrap;
  }

  section.status .na {
    color: #9AA3AD;
  }

  /* ---- Coverage slide ---- */
  section.coverage {
    padding: 40px 56px 60px;
    align-content: start !important;
  }

  section.coverage h1 {
    font-size: 30px;
    margin-bottom: 2px;
  }

  section.coverage h2 {
    font-size: 17px;
    margin-bottom: 14px;
  }

  .cov-grid {
    display: grid;
    grid-template-columns: 1.15fr 1fr;
    gap: 28px;
  }

  .cov-panel h3 {
    font-size: 15px;
    font-weight: 700;
    color: var(--ot-navy);
    margin: 0 0 2px;
  }

  .cov-panel p {
    font-size: 12px;
    color: var(--ot-gray);
    margin: 0 0 6px;
  }

  .cov-panel img {
    width: 100%;
    margin-bottom: 6px;
    display: block;
    background: #ffffff;
  }

  .tiles {
    display: flex;
    gap: 10px;
    margin: 4px 0 14px;
  }

  .tile {
    flex: 1;
    border-top: 3px solid var(--ot-blue);
    background: #F5F8FB;
    padding: 6px 10px;
  }

  .tile b {
    display: block;
    font-size: 24px;
    line-height: 1.1;
    color: var(--ot-text);
  }

  .tile span {
    display: block;
    font-size: 12px;
    color: var(--ot-gray);
  }
---

<!-- _class: status -->

# Brain & Mind: new datasets for the Platform

## Two sources ready for the PPP 26.12 run, three in the pipeline

<table>
<thead>
<tr><th>Dataset</th><th>Status</th><th>Biosamples</th><th>Genes</th><th>Diseases</th><th>Next</th></tr>
</thead>
<tbody>
<tr>
<td><strong>BigBrain</strong> eQTL + sQTL<span class="src">#4462 · 104k credible sets</span></td>
<td><span class="pill done">Ready</span></td>
<td>Bulk brain<span class="src">n = 10.7k, EUR</span></td>
<td class="num"><strong>19.2k</strong> genes<span class="src">4.0k new to brain QTLs</span></td>
<td class="na">n/a (molQTL)</td>
<td>In the PPP 26.12 run</td>
</tr>
<tr>
<td><strong>PGC</strong> psychiatric GWAS<span class="src">#4580 · 3,854 credible sets</span></td>
<td><span class="pill done">Ready</span></td>
<td class="na">n/a (GWAS)</td>
<td class="num"><strong>2.0k</strong> L2G genes<span class="src">5.0k gene–disease pairs</span></td>
<td><strong>27</strong> diseases<span class="src">with credible sets · 47 mapped</span></td>
<td>In the PPP 26.12 run. 41 GWAS Catalog studies to replace</td>
</tr>
<tr>
<td><strong>PsychAD</strong> single-cell eQTL<span class="src">via eQTL Catalogue · 34 datasets</span></td>
<td><span class="pill wait">Pending upstream</span></td>
<td><strong>26</strong> cell subtypes<span class="src">8 classes · 373–1,655 donors</span></td>
<td class="na">After ingestion</td>
<td class="na">n/a (molQTL)</td>
<td>Comes in with the next eQTL Catalogue release</td>
</tr>
<tr>
<td><strong>PsychENCODE</strong> molQTL<span class="src">#4564 · not in eQTL Catalogue</span></td>
<td><span class="pill blocked">Blocked</span></td>
<td>Brain, incl. fetal</td>
<td class="na">Scope TBD</td>
<td class="na">n/a (molQTL)</td>
<td>Licence answer from the fetal dataset authors, then our own ingestion</td>
</tr>
<tr>
<td><strong>CSF proteomics</strong> pQTL<span class="src">#4461 · GWAS Catalog pQTLs</span></td>
<td><span class="pill scope">To decide</span></td>
<td>CSF</td>
<td class="na">Scope TBD</td>
<td class="na">n/a (pQTL)</td>
<td>Decide if CSF belongs in the mental-disorders scope</td>
</tr>
</tbody>
</table>

---

<!-- _class: coverage -->

# What the new data covers

## PPP 26.12 run (ppp-2612-1), after validation

<div class="cov-grid">
<div class="cov-panel">
<h3>PGC: L2G gene–disease pairs, top 12 psychiatric diseases</h3>
<p>5,034 pairs, 1,998 genes · <b>1,933 pairs found by no other GWAS</b> · 28% of psychiatric GWAS credible sets</p>
<img src="figures/pgc_disease_pairs.svg" alt="" />
</div>
<div class="cov-panel">
<h3>BigBrain: genes with a brain credible set</h3>
<p>Bulk brain (UBERON_0000955) · 41k eQTL + 62k sQTL credible sets</p>
<div class="tiles">
<div class="tile"><b>16,325</b><span>eQTL genes</span></div>
<div class="tile"><b>8,187</b><span>sQTL genes</span></div>
<div class="tile"><b>+3,986</b><span>new to brain QTLs</span></div>
</div>
<h3>BigBrain is the 2nd largest brain QTL source</h3>
<p>Brain QTL genes in the run: 22,019 → <b>26,005</b> with BigBrain</p>
<img src="figures/brain_qtl_projects.svg" alt="" />
</div>
</div>
