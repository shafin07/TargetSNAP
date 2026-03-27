# TargetSNAP

**SNP-aware miRNA target prediction using TargetScan Context++ scores.**

TargetSNAP is a web-based bioinformatics tool that predicts how a single-nucleotide polymorphism (SNP) in a gene's 3′UTR affects miRNA binding.  It runs [TargetScan 7.0](http://www.targetscan.org) locally on both the reference and alternate allele UTR sequences and classifies each miRNA target site as **Gain of Function (GOF)** or **Loss of Function (LOF)**.

---

## Features

| Feature | Description |
|---------|-------------|
| **SNP lookup** | Enter a dbSNP rs ID to retrieve chromosomal position, alleles, and overlapping genes from Ensembl GRCh37. |
| **Transcript selection** | View all transcript isoforms whose 3′UTR overlaps the SNP. |
| **Strand-aware allele handling** | For genes on the minus strand the alleles are complemented and the SNP offset is counted from the UTR end. |
| **TargetScan integration** | Runs TargetScan 7.0 Perl scripts to compute Context++ scores for both ref and alt UTR sequences. |
| **GOF / LOF classification** | miRNA families are classified based on whether binding is gained, lost, or unchanged between alleles. |
| **Interactive results** | Click any miRNA to expand full site-level details in a dropdown. |

---

## Quick start

### Prerequisites

* **Python 3.9+**
* **Perl 5** (for TargetScan scripts)
* The actual TargetScan 7.0 Perl scripts (see [Setup](#setup) below)

### Setup

```bash
# 1. Clone the repository
git clone https://github.com/shafin07/TargetSNAP.git
cd TargetSNAP

# 2. Create a virtual environment and install dependencies
python -m venv venv
source venv/bin/activate   # Windows: venv\Scripts\activate
pip install -r requirements.txt

# 3. Download TargetScan 7.0 scripts
#    Visit http://www.targetscan.org/cgi-bin/targetscan/data_download.cgi?db=vert_70
#    and replace the placeholder scripts in scripts/ with the real ones:
#      scripts/targetscan_70.pl
#      scripts/targetscan_70_BL_PCT.pl

# 4. Download the full miRNA family file
#    Replace data/miR_Family_Info.txt with the complete file from TargetScan.

# 5. Run the application
python app.py
```

Open **http://localhost:5000** in your browser.

---

## Usage

1. **Enter an rs ID** (e.g. `rs11552978`) and click **Search**.
2. The tool queries Ensembl GRCh37 for the SNP and lists genes whose 3′UTR overlaps.
3. **Click a gene/transcript** to run TargetScan on the reference and alternate UTR.
4. Results are grouped into **Loss of Function**, **Gain of Function**, and **Unchanged**.
5. **Click any miRNA family** to expand the full Context++ site-level details.

---

## How it works

### Strand handling

| Gene strand | Allele handling | SNP offset formula |
|:-----------:|-----------------|-------------------|
| **+** | Alleles used as reported | `snp_pos − utr_start` |
| **−** | Alleles are complemented (e.g. G>A → C>T) | `utr_end − snp_pos` |

### Classification

| Category | Condition |
|----------|-----------|
| **Loss of Function** | miRNA site present in ref but absent or weaker (higher score) in alt |
| **Gain of Function** | miRNA site absent or weaker in ref but present or stronger (lower score) in alt |
| **Unchanged** | No meaningful change in binding affinity |

---

## Project structure

```
TargetSNAP/
├── app.py                        # Flask web application
├── config.py                     # Paths and configuration
├── requirements.txt              # Python dependencies
├── targetscan/
│   ├── __init__.py
│   ├── snp_lookup.py             # Ensembl GRCh37 SNP & gene lookup
│   ├── utr_retrieval.py          # UCSC hg19 DAS UTR sequence fetch
│   ├── strand_handler.py         # Strand-aware allele & offset logic
│   ├── targetscan_runner.py      # TargetScan Perl script wrapper
│   └── analysis.py               # GOF / LOF classification
├── scripts/
│   ├── targetscan_70.pl          # TargetScan prediction (placeholder)
│   └── targetscan_70_BL_PCT.pl   # TargetScan BL/PCT (placeholder)
├── data/
│   └── miR_Family_Info.txt       # miRNA seed family info (sample)
├── templates/
│   └── index.html                # Main web page
├── static/
│   ├── css/style.css
│   └── js/main.js
└── tests/
    └── test_app.py               # Unit tests
```

---

## Data sources

| Data | Source |
|------|--------|
| SNP information | [Ensembl GRCh37 REST API](https://grch37.rest.ensembl.org) |
| UTR sequences | [UCSC Genome Browser DAS (hg19)](https://genome.ucsc.edu) |
| miRNA families | [TargetScan 7.0](http://www.targetscan.org) |

---

## Running tests

```bash
python -m pytest tests/ -v
```

---

## License

This project is provided for academic and research use.  TargetScan scripts are subject to their own [license terms](http://www.targetscan.org).
