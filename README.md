# TargetSNAP

**SNP Impact Analysis on miRNA Binding Sites**

TargetSNAP is a locally deployable web tool that predicts how single-nucleotide polymorphisms (SNPs) in 3′ UTRs affect microRNA binding. It integrates **TargetScan 7** context++ scoring with **RNAhybrid** thermodynamic validation and multi-source clinical annotations.

---

## Features

- **Dual-tool prediction** — TargetScan 7 (seed + context++) and RNAhybrid (MFE + p-value)
- **LOF / GOF / Neutral classification** with composite priority scoring
- **ClinVar** clinical significance (NCBI eUtils)
- **GWAS Catalog** trait associations (EBI REST)
- **Conservation** score — phastCons 46-way vertebrate (UCSC REST)
- **Tissue expression filter** — 15 tissues (TCGA / miRmine data)
- **eQTL / DAE evidence** from GTEx v10
- **Interactive duplex visualisation** with SNP highlighting
- **Export** to CSV, TSV, or JSON
- **Persistent caching** — first run ~30–60 s, cached replay ~100 ms

---

## Prerequisites

| Tool | Version | Notes |
|------|---------|-------|
| **Python** | 3.10+ | Tested on 3.11–3.13 |
| **Perl** | 5.x | Required by TargetScan scripts |
| **RNAhybrid** | ≥ 2.1 | Via WSL on Windows, or native on Linux/macOS |

### Windows — RNAhybrid via WSL

```bash
# Install WSL (if not already)
wsl --install

# Inside WSL, install RNAhybrid
sudo apt update && sudo apt install -y rnahybrid
```

Verify: `wsl RNAhybrid -h` should print usage info.

### Linux / macOS

Install RNAhybrid via your package manager or from source:
```bash
sudo apt install rnahybrid        # Debian/Ubuntu
# or
brew install rnahybrid             # macOS (if available)
```

---

## Installation

```bash
# Clone the repository
git clone https://github.com/shafin07/TargetSNAP.git
cd TargetSNAP

# Install Python dependencies
pip install -r requirements_web.txt
```

---

## Running

```bash
python app.py
```

Open **http://127.0.0.1:5000** in your browser.

### Quick Start

1. Enter an **rs ID** (e.g., `rs10318`) in the search box
2. Select a **gene** from the results
3. Select a **transcript** — the analysis runs automatically
4. View **LOF / GOF / Neutral** miRNA predictions with scores and MFE values
5. Use the **tissue filter** to highlight tissue-relevant miRNAs
6. **Download** results as CSV, TSV, or JSON

---

## Project Structure

```
TargetSNAP/
├── app.py                          # Flask backend (API endpoints)
├── targetsnap_web_utils_clean.py   # Core genomic data + TargetScan pipeline logic
├── requirements_web.txt            # Python dependencies
├── run_server.bat                  # Windows batch launcher
├── templates/
│   └── web_index.html              # Frontend (single-page app)
├── genomic_data/
│   ├── snp_gene_mapping.json       # Local SNP→gene cache
│   └── gene_transcripts.json       # Local transcript cache
├── TargetScan/
│   ├── targetscan_70/              # TargetScan 7 Perl scripts + data
│   ├── TargetScan7_context_scores/ # Context++ scoring scripts + data
│   └── TSHuman_7_hg19_3UTRs.gff   # 3′ UTR annotations (GRCh37)
├── targetsnap_cache/               # Persistent prediction cache (auto-created)
├── targetsnap_logs/                # Execution logs (auto-created)
├── test_full_flow.py               # Integration test
├── test_rs10318.py                 # rs10318 validation test
└── test_ui.py                      # UI smoke test
```

---

## API Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/api/search-rs` | POST | Map rs ID → genes + alleles |
| `/api/gene-transcripts` | POST | Get transcripts for a gene |
| `/api/compare-alleles` | POST | **Main analysis** — TargetScan + RNAhybrid + classification |
| `/api/snp-annotations` | POST | ClinVar, GWAS Catalog, conservation |
| `/api/eqtl-dae` | POST | GTEx eQTL + DAE evidence |
| `/api/export-results` | POST | Download as CSV / TSV / JSON |
| `/api/rnahybrid` | POST | Standalone RNAhybrid for one miRNA pair |
| `/api/debug-targetscan` | POST | Pipeline diagnostics |
| `/api/health` | GET | Health check |

---

## How It Works

1. **SNP Lookup** — Maps the rs ID to genes via Ensembl VEP (GRCh37)
2. **Sequence Retrieval** — Gets 3′ UTR via local TargetScan index → BioMart → REST API (multi-stage fallback)
3. **Allele Substitution** — Applies REF and ALT alleles with strand-aware handling
4. **TargetScan Pipeline** — Runs seed matching + context++ scoring on both alleles (batched, parallel)
5. **RNAhybrid** — Computes MFE for each miRNA–UTR duplex (8 parallel threads)
6. **Classification** — Compares REF vs ALT to classify LOF / GOF / Neutral
7. **Priority Scoring** — Ranks by Δcontext++, site type, binding change, ΔMFE, and p-value
8. **Annotations** — ClinVar, GWAS, conservation fetched asynchronously

---

## Priority Score

| Component | Weight | Description |
|-----------|--------|-------------|
| Context++ Δ | ×100 | Allelic score change magnitude |
| Site type | ×5 | 8mer=4, 7mer-m8/1a=3, 6mer=2 |
| Binding change | +50 | If site gained, lost, or type changed |
| \|ΔMFE\| | ×3 | RNAhybrid thermodynamic shift |
| p-value bonus | +10/+3 | Per allele: +10 if p<0.05, +3 if p<0.20 |

---

## Tests

```bash
# Run all tests
python -m pytest test_full_flow.py test_rs10318.py test_ui.py -v

# Or individually
python test_full_flow.py
python test_rs10318.py
```

> **Note**: Tests require the server to be running (`python app.py` in a separate terminal).

---

## Reference Genome

TargetSNAP uses **GRCh37/hg19** coordinates, matching the TargetScan 7 reference files.

---

## License

[MIT License](LICENSE)

---

## Citation

If you use TargetSNAP in your research, please cite:

> Ahmad, S.S. (2026). TargetSNAP: An Integrated Local Web Tool for Predicting SNP-Mediated Disruption of miRNA Binding Sites Using TargetScan 7 and RNAhybrid. *[Journal]*, *[Volume]*, *[Pages]*.
