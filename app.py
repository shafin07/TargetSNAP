"""
app.py – TargetSNAP Flask web application.

Routes
------
GET  /          → index page (rs ID input form)
POST /analyze   → run TargetSNAP analysis, return results page
"""

from __future__ import annotations

from markupsafe import Markup, escape

from flask import Flask, render_template, request

from targetscan.analysis import analyze_snp, WINDOW, FILTER_DISTANCE

app = Flask(__name__)

# ---------------------------------------------------------------------------
# Jinja2 helper: annotated sequence display
# ---------------------------------------------------------------------------

DISPLAY_HALF = 45  # nucleotides to show on each side of the SNP


def _annotate_seq(seq: str, snp_idx: int, affected_matches: list[dict]) -> Markup:
    """Return a Markup string with per-character <span> highlights.

    - SNP position     → ``highlight-snp``
    - Seed-match range → ``highlight-match`` (or ``highlight-snp-match`` if
                          the SNP falls inside the seed)
    """
    start = max(0, snp_idx - DISPLAY_HALF)
    end = min(len(seq), snp_idx + DISPLAY_HALF + 1)

    # Build a set of positions covered by affected seed matches
    match_pos: set[int] = set()
    for m in affected_matches:
        for i in range(max(m["start"], start), min(m["end"], end)):
            match_pos.add(i)

    parts: list[str] = []
    if start > 0:
        parts.append('<span class="text-muted">…</span>')

    for i in range(start, end):
        char = str(escape(seq[i]))
        is_snp = i == snp_idx
        in_match = i in match_pos

        if is_snp and in_match:
            parts.append(f'<span class="highlight-snp-match">{char}</span>')
        elif is_snp:
            parts.append(f'<span class="highlight-snp">{char}</span>')
        elif in_match:
            parts.append(f'<span class="highlight-match">{char}</span>')
        else:
            parts.append(char)

    if end < len(seq):
        parts.append('<span class="text-muted">…</span>')

    return Markup("".join(parts))


@app.template_filter("annotate_ref")
def annotate_ref(analysis: dict) -> Markup:
    """Annotate the reference sequence, highlighting LOF sites."""
    return _annotate_seq(
        analysis["ref_seq"], analysis["snp_idx"], analysis["lof"]
    )


@app.template_filter("annotate_alt")
def annotate_alt(analysis: dict) -> Markup:
    """Annotate the alternate sequence, highlighting GOF sites."""
    return _annotate_seq(
        analysis["alt_seq"], analysis["snp_idx"], analysis["gof"]
    )

# ---------------------------------------------------------------------------
# Routes
# ---------------------------------------------------------------------------


@app.route("/")
def index():
    return render_template("index.html")


@app.route("/analyze", methods=["POST"])
def analyze():
    rsid = request.form.get("rsid", "").strip()
    if not rsid:
        return render_template("index.html", error="Please enter an rs ID.")
    result = analyze_snp(rsid)
    return render_template(
        "results.html", result=result, window=WINDOW, filter_distance=FILTER_DISTANCE
    )


if __name__ == "__main__":
    import os
    debug = os.environ.get("FLASK_DEBUG", "0") == "1"
    app.run(debug=debug, host="0.0.0.0", port=5000)
