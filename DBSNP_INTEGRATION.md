# dbSNP Integration - TargetSNAP Web

## What Changed

The web interface now integrates with **real dbSNP data** from NCBI and Ensembl APIs!

### How It Works

1. **SNP Search**: When you enter an rs ID, it queries Ensembl REST API
2. **Gene Mapping**: Returns real genes affected by that SNP
3. **Transcripts**: Fetches actual transcript information from Ensembl
4. **Sequences**: Retrieves real mRNA/cDNA sequences for analysis

### Benefits

✅ Real SNP data (not just test data)  
✅ Any rs ID from dbSNP will work  
✅ Actual gene and transcript information  
✅ Real sequences for accurate TargetScan predictions  

## Setup

### Install requests library

```powershell
pip install requests
```

Or install all web requirements:

```powershell
pip install -r requirements_web.txt
```

## Restart the Server

1. Stop the current server (Ctrl+C in terminal)
2. Run again:
   ```powershell
   python app.py
   ```
3. Open browser to `http://localhost:5000`

## Try Real SNPs

Now you can search for ANY real dbSNP rs ID:

- `rs11552978` - Will now work! (if it exists in dbSNP)
- `rs1333049` - Common variant
- `rs662` - PON1 gene variant
- `rs429358` - APOE gene variant

## How It Works (Technical)

### API Chain

```
User enters rs ID
    ↓
Queries Ensembl VEP API (Variant Effect Predictor)
    ↓
Returns variant consequences and gene info
    ↓
Fetches transcript data from Ensembl lookup API
    ↓
Retrieves cDNA sequences from Ensembl sequence API
    ↓
Runs local TargetScan predictions
```

### Strand Handling

For dbSNP results:
- **+ strand**: Position calculated from start
- **- strand**: Position calculated from end with reverse complement

### Data Caching

Results are cached locally to speed up repeated searches and reduce API calls.

## What if dbSNP API is Down?

The system has fallback to local test data:
- `rs12345678` (BRCA2)
- `rs61733396` (CDKN2A)

These local SNPs will always work even if APIs are unavailable.

## Error Handling

If a specific rs ID not found:
- Shows error message
- Tries local database
- Falls back gracefully

## Limitations

- API calls may take 2-5 seconds per search
- Very rare SNPs might not have complete gene annotations
- Some SNPs might have multiple gene hits

## Next Steps

1. Restart server with `python app.py`
2. Try searching for real SNPs!
3. Check the console for API response times

Questions? The system logs API queries to help troubleshoot.
