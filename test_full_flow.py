import requests, json, time

base = 'http://127.0.0.1:5000'

# Step 1: Search (POST to /api/search-rs)
print('=== Step 1: Search ===')
t0 = time.time()
r = requests.post(f'{base}/api/search-rs', json={'rs_id': 'rs10318'})
print(f'  Status: {r.status_code}, Time: {time.time()-t0:.1f}s')
data = r.json()
genes = data.get('genes', [])
print(f'  RS: {data.get("rs_id", "??")}')
print(f'  Genes: {[g["gene_name"] for g in genes]}')

if not genes:
    print(f'NO GENES - stopping. Response: {json.dumps(data, indent=2)[:500]}')
    exit()

gene = genes[0]
print(f'  Gene: {gene["gene_name"]} ({gene["gene_id"]})')

# Step 2: Get transcripts (POST to /api/gene-transcripts)
print('\n=== Step 2: Transcripts ===')
t0 = time.time()
r = requests.post(f'{base}/api/gene-transcripts', json={'gene_id': gene['gene_id']})
print(f'  Status: {r.status_code}, Time: {time.time()-t0:.1f}s')
tdata = r.json()
transcripts = tdata.get('transcripts', [])
print(f'  Transcripts: {len(transcripts)}')
for t in transcripts[:5]:
    print(f'    {t["transcript_id"]}: {t.get("transcript_name","")} UTR={t.get("utr_length","NA")}bp')

if not transcripts:
    print('NO TRANSCRIPTS - stopping')
    exit()

# Step 3: Compare alleles
transcript = transcripts[0]
print(f'\n=== Step 3: Compare alleles ({transcript["transcript_id"]}) ===')
t0 = time.time()
r = requests.post(f'{base}/api/compare-alleles', json={
    'rs_id': data['rs_id'],
    'gene_id': gene['gene_id'],
    'transcript_id': transcript['transcript_id']
})
elapsed = time.time()-t0
print(f'  Status: {r.status_code}, Time: {elapsed:.1f}s')
print(f'  Response size: {len(r.content)} bytes ({len(r.content)//1024}KB)')
if r.status_code == 200:
    cdata = r.json()
    lof = cdata.get('loss_of_function', [])
    gof = cdata.get('gain_of_function', [])
    neutral = cdata.get('neutral', [])
    print(f'  LOF: {len(lof)}, GOF: {len(gof)}, Neutral: {len(neutral)}')
    print(f'  total_ref_targets: {cdata.get("total_ref_targets")}')
    print(f'  total_mut_targets: {cdata.get("total_mut_targets")}')
    for m in lof:
        print(f'    LOF: {m.get("mirna_id")} ref={m.get("ref_score")} mut={m.get("mut_score")} site={m.get("site_type")}')
    for m in gof:
        print(f'    GOF: {m.get("mirna_id")} ref={m.get("ref_score")} mut={m.get("mut_score")} site={m.get("site_type")}')
    perf = cdata.get('performance_ms', {})
    print(f'\n  Performance: {json.dumps(perf, indent=4)}')
else:
    print(f'  Error: {r.text[:500]}')
