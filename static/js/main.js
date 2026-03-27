/* TargetSNAP – main client-side logic */

"use strict";

// Cache SNP info between steps so we don't refetch it.
let cachedSnpInfo = null;
let cachedGenes = [];

// ── Helpers ──────────────────────────────────────────────────────────────

function show(id) {
  document.getElementById(id).classList.remove("hidden");
}

function hide(id) {
  document.getElementById(id).classList.add("hidden");
}

function setError(id, msg) {
  const el = document.getElementById(id);
  el.textContent = msg;
  show(id);
}

function clearError(id) {
  const el = document.getElementById(id);
  el.textContent = "";
  hide(id);
}

async function apiFetch(url, options) {
  const resp = await fetch(url, options);
  const data = await resp.json();
  if (!resp.ok) {
    throw new Error(data.error || `Request failed (${resp.status})`);
  }
  return data;
}

// ── Step 1: Search SNP ──────────────────────────────────────────────────

document.getElementById("snp-form").addEventListener("submit", async (e) => {
  e.preventDefault();
  const rsId = document.getElementById("rs-input").value.trim();
  if (!rsId) return;

  clearError("search-error");
  hide("gene-section");
  hide("results-section");

  try {
    // Fetch SNP info
    cachedSnpInfo = await apiFetch(`/api/snp/${encodeURIComponent(rsId)}`);

    const alleles = cachedSnpInfo.alleles || [];
    document.getElementById("snp-info").innerHTML =
      `<strong>${cachedSnpInfo.id}</strong> &mdash; ` +
      `chr${cachedSnpInfo.chr}:${cachedSnpInfo.position} ` +
      `(${alleles.join(" / ")})`;
    show("snp-info");

    // Fetch overlapping genes
    cachedGenes = await apiFetch(`/api/genes/${encodeURIComponent(rsId)}`);
    renderGeneList(cachedGenes, cachedSnpInfo);
  } catch (err) {
    setError("search-error", err.message);
  }
});

// ── Step 2: Gene / transcript list ──────────────────────────────────────

function renderGeneList(genes, snpInfo) {
  const ul = document.getElementById("gene-list");
  ul.innerHTML = "";

  genes.forEach((g) => {
    const li = document.createElement("li");
    const strandClass = g.strand === "+" ? "strand-plus" : "strand-minus";
    li.innerHTML =
      `<span class="gene-name">${esc(g.gene_symbol)}</span>` +
      `<span class="tx-id">${esc(g.transcript_id)}</span>` +
      `<span class="strand-tag ${strandClass}">${esc(g.strand)} strand</span>`;

    li.addEventListener("click", () => runAnalysis(g, snpInfo));
    ul.appendChild(li);
  });

  show("gene-section");
  clearError("gene-error");
}

// ── Step 3: Run analysis ────────────────────────────────────────────────

async function runAnalysis(gene, snpInfo) {
  clearError("results-error");
  show("loading");

  const alleles = snpInfo.alleles || [];
  const refAllele = alleles[0] || "";
  const altAllele = alleles.length > 1 ? alleles[1] : "";

  const body = {
    rs_id: snpInfo.id,
    transcript_id: gene.transcript_id,
    gene_symbol: gene.gene_symbol,
    strand: gene.strand,
    utr_start: gene.utr_start,
    utr_end: gene.utr_end,
    snp_position: snpInfo.position,
    ref_allele: refAllele,
    alt_allele: altAllele,
    chr: snpInfo.chr,
  };

  try {
    const res = await apiFetch("/api/analyze", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(body),
    });
    renderResults(res);
  } catch (err) {
    setError("results-error", err.message);
  } finally {
    hide("loading");
  }
}

// ── Render results ──────────────────────────────────────────────────────

function renderResults(res) {
  // Meta info
  const effRef = res.eff_ref_allele || res.ref_allele;
  const effAlt = res.eff_alt_allele || res.alt_allele;
  document.getElementById("results-meta").innerHTML = `
    <table class="results-meta-table">
      <tr><td>Gene</td><td>${esc(res.gene_symbol)} (${esc(res.transcript_id)})</td></tr>
      <tr><td>Strand</td><td>${esc(res.strand)}</td></tr>
      <tr><td>Allele change</td><td>${esc(res.ref_allele)} &rarr; ${esc(res.alt_allele)}
          ${res.strand === "-" ? `(effective: ${esc(effRef)} &rarr; ${esc(effAlt)})` : ""}</td></tr>
      <tr><td>UTR length</td><td>${res.utr_length} nt</td></tr>
      <tr><td>SNP offset in UTR</td><td>${res.snp_offset}</td></tr>
    </table>
    ${res.note ? `<p class="hint">${esc(res.note)}</p>` : ""}`;

  // Counts
  document.getElementById("lof-count").textContent = res.lof.length;
  document.getElementById("gof-count").textContent = res.gof.length;
  document.getElementById("unchanged-count").textContent = res.unchanged.length;

  // Lists
  fillMirnaList("lof-list", res.lof);
  fillMirnaList("gof-list", res.gof);
  fillMirnaList("unchanged-list", res.unchanged);

  show("results-section");
}

function fillMirnaList(containerId, items) {
  const container = document.getElementById(containerId);
  container.innerHTML = "";

  if (items.length === 0) {
    container.innerHTML = '<p class="hint">None</p>';
    return;
  }

  items.forEach((item) => {
    const div = document.createElement("div");
    div.className = "mirna-item";

    // Summary row (always visible)
    const summary = document.createElement("div");
    summary.className = "mirna-summary";
    summary.innerHTML =
      `<span class="family">${esc(item.miRNA_family)}</span>` +
      `<span class="scores">ref: ${fmt(item.ref_score)} &rarr; alt: ${fmt(item.alt_score)} (&Delta; ${fmt(item.delta)})</span>`;

    // Details panel (hidden until click)
    const details = document.createElement("div");
    details.className = "mirna-details";
    details.innerHTML = buildDetailsHTML(item);

    summary.addEventListener("click", () => {
      details.classList.toggle("open");
    });

    div.appendChild(summary);
    div.appendChild(details);
    container.appendChild(div);
  });
}

function buildDetailsHTML(item) {
  let html = "";

  if (item.ref_details && item.ref_details.length) {
    html += "<h4>Reference allele sites</h4>";
    html += detailsTable(item.ref_details);
  }

  if (item.alt_details && item.alt_details.length) {
    html += "<h4 style='margin-top:0.5rem'>Alternate allele sites</h4>";
    html += detailsTable(item.alt_details);
  }

  if (!html) {
    html = "<p>No site-level details available.</p>";
  }
  return html;
}

function detailsTable(rows) {
  if (!rows.length) return "";
  const keys = Object.keys(rows[0]);
  let html = "<table><thead><tr>";
  keys.forEach((k) => { html += `<th>${esc(k)}</th>`; });
  html += "</tr></thead><tbody>";
  rows.forEach((r) => {
    html += "<tr>";
    keys.forEach((k) => { html += `<td>${esc(String(r[k] ?? ""))}</td>`; });
    html += "</tr>";
  });
  html += "</tbody></table>";
  return html;
}

// ── Utilities ───────────────────────────────────────────────────────────

function esc(s) {
  const el = document.createElement("span");
  el.textContent = s;
  return el.innerHTML;
}

function fmt(n) {
  if (n === 0 || n === null || n === undefined) return "0";
  return Number(n).toFixed(4);
}
