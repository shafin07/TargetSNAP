import requests, json, time

payload = {
    'gene_id': 'ENSG00000166923',
    'gene_name': 'GREM1',
    'transcript_id': 'ENST00000300177',
    'chromosome': '15',
    'position': 33025979,
    'strand': '+',
    'ref_allele': 'C',
    'alt_allele': 'T',
    'rs_id': 'rs10318'
}
s = time.time()
r = requests.post('http://127.0.0.1:5000/api/compare-alleles', json=payload, timeout=900)
elapsed = time.time() - s
d = r.json()
print(f'compare-alleles: {elapsed:.1f}s, status: {r.status_code}')
if 'error' in d:
    print(f'Error: {d["error"]}')
else:
    summary = d.get('summary', {})
    print(f'Summary: {json.dumps(summary)}')
    comp = d.get('comparison', [])
    print(f'Total comparisons: {len(comp)}')
    lof = d.get('loss_of_function', [])
    gof = d.get('gain_of_function', [])
    neutral = d.get('neutral', [])
    print(f'LOF: {len(lof)}, GOF: {len(gof)}, Neutral: {len(neutral)}')
    
    # Check for the expected miRNAs
    all_effects = lof + gof + neutral
    for mirna_name in ['hsa-miR-331-5p', 'hsa-miR-4678']:
        found = [e for e in all_effects if e.get('mirna_id') == mirna_name]
        if found:
            e = found[0]
            print(f'{mirna_name}: effect={e.get("effect_type","?")} ref_score={e.get("ref_score","?")} mut_score={e.get("mut_score","?")} site={e.get("site_type","")}')
        else:
            print(f'{mirna_name}: NOT FOUND in results')

    # Show first few LOF and GOF
    for label, items in [('LOF', lof[:3]), ('GOF', gof[:3])]:
        for item in items:
            print(f'  {label}: {item.get("mirna_id","")} ref={item.get("ref_score","")} mut={item.get("mut_score","")} site={item.get("site_type","")}')
