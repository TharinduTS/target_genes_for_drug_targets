# Target genes for drug targets

# 1) Mapping genes to drug targets

## 1-I Introduction

I am trying to see whether there is a connection between the success of drugs (whether they made into testing phases 2,3,4 etc) and the cell specificity of their drug targets.

I started by downloading the list of drugs with their MOA (Mechanism of Action) for humans from 'https://www.ebi.ac.uk/chembl/explore/drug_mechanisms/STATE_ID:YB5SeM-jQgB2RiRcCSv8sA%3D%3D'

The download looks like following

```tsv
"Parent Molecule ChEMBL ID"     "Parent Molecule Name"  "Parent Molecule Type"  "Max Phase"     "First Approval"        "USAN Stem"     "Smiles"        "Mechanism of Action"   "Mechanism Comment"       "Selectivity Comment"   "Target ChEMBL ID"      "Target Name"   "Action Type"   "Target Type"   "Target Organism"       "Binding Site Name"     "Binding Site Comment"    "References"    "Synonyms"      "Parent Molecule ChEMBL ID"
"CHEMBL360055"  "GALLAMINE"     "Small molecule"        "4"     "1982"  ""      "CC[N+](CC)(CC)CCOc1cccc(OCC[N+](CC)(CC)CC)c1OCC[N+](CC)(CC)CC" "Muscle-type nicotinic acetylcholine receptor antagonist" "non-depolarizing"      ""      "CHEMBL2362997" "Muscle-type nicotinic acetylcholine receptor"  "ANTAGONIST"    "PROTEIN COMPLEX GROUP" "Homo sapiens"  ""        ""      "{'CHEMBL1200993': [{'ref_id': '9780702034718 PP. 164', 'ref_url': 'http://www.isbnsearch.org/isbn/9780702034718', 'ref_type': 'ISBN'}]}"       "Gallamine triethiodide (INN)|"   "CHEMBL360055"
"CHEMBL1259"    "METOCURINE"    "Small molecule"        "4"     "1982"  "'-ium'"        "COc1ccc2cc1Oc1cc3c(cc1OC)CC[N+](C)(C)[C@H]3Cc1ccc(cc1)Oc1c(OC)c(OC)cc3c1[C@@H](C2)[N+](C)(C)CC3" "Muscle-type nicotinic acetylcholine receptor antagonist"       ""      ""      "CHEMBL2362997" "Muscle-type nicotinic acetylcholine receptor"  "ANTAGONIST"    "PROTEIN COMPLEX GROUP"   "Homo sapiens"  ""      ""      "{'CHEMBL1739': [{'ref_id': '18633030', 'ref_url': 'http://europepmc.org/abstract/MED/18633030', 'ref_type': 'PubMed'}, {'ref_id': 'Dimethyltubocurarine', 'ref_url': 'http://en.wikipedia.org/wiki/Dimethyltubocurarine', 'ref_type': 'Wikipedia'}]}"    "Metocurine iodide (USAN, USP)|Dimethyltubocurarinium chloride (INN)|"    "CHEMBL1259"
```

*and renamed the file to drug_targets.tsv*

Then I wanted to find genes each of these drugs affect.

for that, I can use the "Target ChEMBL ID" and search it in their target database like 'https://www.ebi.ac.uk/chembl/explore/target/CHEMBL2362997' , changing Target ChEMBL ID. It has the genes those drugs affect under the section '#GeneCrossRefs'

*Note that this is not the same as "Parent Molecule ChEMBL ID"*

Because I had thousands of targets, I wanted to automate it and wrote the script 'chembl_target_scraper.py'

This automatically adds genes affected by each gene in front of each gene seperated by commas.

## 1-II Script

chembl_target_scraper.py

```py

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import asyncio
import json
import os
import re
from typing import Dict, List, Any, Optional

import aiohttp
import pandas as pd
from tqdm import tqdm


# =================
# I/O: table utils
# =================

def read_table(path: str) -> pd.DataFrame:
    lower = path.lower()
    if lower.endswith((".tsv", ".txt")):
        return pd.read_csv(path, sep="\t")
    if lower.endswith(".csv"):
        return pd.read_csv(path)
    if lower.endswith(".xlsx"):
        return pd.read_excel(path, engine="openpyxl")
    if lower.endswith(".xls"):
        return pd.read_excel(path, engine="xlrd")
    raise ValueError("Unsupported input format. Use .tsv/.csv/.xlsx/.xls")


def write_table(df: pd.DataFrame, path: str) -> None:
    lower = path.lower()
    if lower.endswith((".tsv", ".txt")):
        df.to_csv(path, sep="\t", index=False)
    elif lower.endswith(".csv"):
        df.to_csv(path, index=False)
    elif lower.endswith(".xlsx"):
        df.to_excel(path, index=False, engine="openpyxl")
    else:
        raise ValueError("Unsupported output format. Use .tsv/.csv/.xlsx")


# ===============================
# API JSON → cross-ref extraction
# ===============================

# Normalize different source labels to canonical keys
SRC_NORMALIZE = {
    # Genes
    "HGNC": "HGNC",
    # Proteins
    "UNIPROT": "UniProt",
    "UNIPROTKB": "UniProt",
    "UNIPROT/UNIPROTKB": "UniProt",
    # Ensembl
    "ENSEMBL": "Ensembl",
    "ENSEMBLGENE": "Ensembl",
    "ENSEMBLGENEID": "Ensembl",
    # Domains
    "PFAM": "Pfam",
    "PFAMA": "Pfam",       # sometimes appears as "Pfam-A"
    "PFAM-A": "Pfam",
    "INTERPRO": "InterPro",
}

CANON_KEYS = ["HGNC", "UniProt", "Ensembl", "Pfam", "InterPro"]

def _dedupe_keep_order(items: List[str]) -> List[str]:
    seen = set()
    out = []
    for it in items:
        if it not in seen:
            out.append(it)
            seen.add(it)
    return out

def _clean_val(val: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]", "", val or "")

def _get_target_object(api_json: Any) -> Optional[Dict[str, Any]]:
    """
    The /target/<ID>.json endpoint typically returns a single target object.
    Be tolerant to slight variations: {'target': {...}} or {'targets': [...]}.
    """
    if isinstance(api_json, dict):
        if "target_components" in api_json:
            return api_json
        if "target" in api_json and isinstance(api_json["target"], dict):
            return api_json["target"]
        if "targets" in api_json and isinstance(api_json["targets"], list) and api_json["targets"]:
            return api_json["targets"][0]
    return None

def extract_xrefs_from_api_json(api_json: Any) -> Dict[str, List[str]]:
    """
    Parse the API JSON for a single target and collect HGNC/UniProt/Ensembl/Pfam/InterPro
    from target_components → target_component_xrefs.
    NOTE (fix): accept both 'xref_src' and 'xref_src_db'; for HGNC prefer 'xref_name' (gene symbol).
    """
    refs: Dict[str, List[str]] = {k: [] for k in CANON_KEYS}
    tgt = _get_target_object(api_json)
    if not tgt:
        return refs

    components = tgt.get("target_components") or []
    for comp in components:
        xrefs = comp.get("target_component_xrefs") or []
        for x in xrefs:
            # Some payloads use 'xref_src', others 'xref_src_db'
            src_raw = (x.get("xref_src") or x.get("xref_src_db") or "").strip()
            if not src_raw:
                continue

            norm_key = SRC_NORMALIZE.get(src_raw.replace(" ", "").replace("-", "").upper())
            if norm_key not in refs:
                continue

            # Value strategy:
            # - HGNC → prefer gene symbol in xref_name; fallback to xref_id.
            # - Others → prefer accession-like xref_id; fallback to xref_name.
            if norm_key == "HGNC":
                val_raw = (x.get("xref_name") or x.get("xref_id") or "").strip()
            else:
                val_raw = (x.get("xref_id") or x.get("xref_name") or "").strip()

            val = _clean_val(val_raw)
            if val:
                refs[norm_key].append(val)

    for k in refs:
        refs[k] = _dedupe_keep_order(refs[k])
    return refs


# ===================
# Async HTTP fetching
# ===================

def build_api_url(target_id: str) -> str:
    # Official REST endpoint (JSON): https://www.ebi.ac.uk/chembl/api/data/target/<ID>.json
    return f"https://www.ebi.ac.uk/chembl/api/data/target/{target_id}.json"

async def fetch_one(
    session: aiohttp.ClientSession,
    target_id: str,
    timeout: float,
    retries: int,
    delay: float
) -> Optional[Dict[str, List[str]]]:
    url = build_api_url(target_id)
    for attempt in range(retries + 1):
        try:
            async with session.get(url, timeout=timeout) as resp:
                if resp.status == 200:
                    data = await resp.json(content_type=None)
                    refs = extract_xrefs_from_api_json(data)
                    if delay > 0:
                        await asyncio.sleep(delay)
                    return refs
                elif resp.status in (404, 410):
                    return None
                else:
                    await asyncio.sleep(0.5 * (attempt + 1))
        except Exception:
            await asyncio.sleep(0.5 * (attempt + 1))
    return None

async def scrape_all(
    target_ids: List[str],
    user_agent: str,
    timeout: float,
    retries: int,
    parallel: int,
    delay: float,
    cache: Dict[str, Dict[str, List[str]]]
) -> Dict[str, Optional[Dict[str, List[str]]]]:
    results: Dict[str, Optional[Dict[str, List[str]]]] = {}
    headers = {"User-Agent": user_agent, "Accept": "application/json"}
    connector = aiohttp.TCPConnector(limit=parallel)

    async with aiohttp.ClientSession(headers=headers, connector=connector) as session:
        sem = asyncio.Semaphore(parallel)

        async def worker(tid: str):
            if tid in cache:
                return tid, cache[tid]
            async with sem:
                refs = await fetch_one(session, tid, timeout, retries, delay)
                if refs is not None:
                    cache[tid] = refs
                return tid, refs

        tasks = [worker(tid) for tid in target_ids]
        for fut in tqdm(asyncio.as_completed(tasks), total=len(tasks), desc="Fetching from ChEMBL API"):
            tid, refs = await fut
            results[tid] = refs
    return results


# =====
# Main
# =====

def main():
    parser = argparse.ArgumentParser(
        description="ChEMBL Target cross-refs (HGNC/UniProt/Ensembl/Pfam/InterPro) via REST API (async, cached)"
    )
    parser.add_argument("-i", "--input", required=True, help="Input table (.tsv/.csv/.xlsx/.xls)")
    parser.add_argument("-o", "--output", required=True, help="Output table (.tsv/.csv/.xlsx)")
    parser.add_argument("-c", "--id-column", required=True, help="Column with Target ChEMBL IDs")
    parser.add_argument("--output-prefix", default="ChEMBL", help="Prefix for new columns (default: ChEMBL)")
    parser.add_argument("--timeout", type=float, default=12.0, help="HTTP timeout per request (s)")
    parser.add_argument("--retries", type=int, default=2, help="Retry attempts per target")
    parser.add_argument("--parallel", type=int, default=12, help="Max concurrent requests")
    parser.add_argument("--delay", type=float, default=0.0, help="Delay after each response per worker (s)")
    parser.add_argument("--user-agent", default="Mozilla/5.0", help="Custom User-Agent")
    parser.add_argument("--cache", default="chembl_cache.json", help="Path to JSON cache file")
    parser.add_argument("--verbose", action="store_true", help="Verbose logging")
    args = parser.parse_args()

    # Load table
    df = read_table(args.input)

    # Normalize headers (strip spaces)
    df.rename(columns={c: c.strip() for c in df.columns}, inplace=True)

    # Tolerant match for id column
    if args.id_column not in df.columns:
        candidates = [c for c in df.columns if c.strip().lower() == args.id_column.strip().lower()]
        if candidates:
            args.id_column = candidates[0]
        else:
            raise KeyError(f"Column '{args.id_column}' not found. Available: {list(df.columns)}")

    # Prepare output columns
    for key in CANON_KEYS:
        col = f"{args.output_prefix}_{key}"
        if col not in df.columns:
            df[col] = ""

    # Unique, clean target IDs
    target_ids = (
        df[args.id_column]
        .dropna()
        .astype(str)
        .map(str.strip)
        .replace({"": None, "nan": None, "None": None})
        .dropna()
        .unique()
        .tolist()
    )

    # Load cache
    cache: Dict[str, Dict[str, List[str]]] = {}
    if os.path.exists(args.cache):
        try:
            with open(args.cache, "r", encoding="utf-8") as f:
                cache = json.load(f)
        except Exception:
            cache = {}

    # Fetch in parallel (API)
    results = asyncio.run(
        scrape_all(
            target_ids=target_ids,
            user_agent=args.user_agent,
            timeout=args.timeout,
            retries=args.retries,
            parallel=args.parallel,
            delay=args.delay,
            cache=cache,
        )
    )

    # Save/refresh cache
    try:
        with open(args.cache, "w", encoding="utf-8") as f:
            json.dump(cache, f, indent=2, ensure_ascii=False)
    except Exception:
        pass

    # Fill rows
    prefix = args.output_prefix
    for idx, row in df.iterrows():
        tid = str(row.get(args.id_column, "")).strip()
        if not tid:
            continue
        refs = results.get(tid)
        if not refs:
            continue
        for key in CANON_KEYS:
            df.at[idx, f"{prefix}_{key}"] = ", ".join(refs.get(key, []))

    # Write output
    write_table(df, args.output)

    if args.verbose:
        print(f"Done. Output written to: {args.output}")

if __name__ == "__main__":
    main()
```

## 1-III CLI help

```txt

chembl_target_scraper.py
------------------------

Fetch HGNC gene symbols and other cross‑references (UniProt, Ensembl, Pfam,
InterPro) for each ChEMBL Target ID in a table. Uses the official ChEMBL
REST API and runs requests in parallel with caching.

USAGE:
  python chembl_target_scraper.py -i INPUT -o OUTPUT -c COLUMN [options]

REQUIRED ARGUMENTS:
  -i, --input FILE          Input table (.tsv, .csv, .xlsx, .xls)
  -o, --output FILE         Output table (.tsv, .csv, .xlsx)
  -c, --id-column NAME      Column containing Target ChEMBL IDs
                            (e.g. "CHEMBL2362997")

OPTIONAL ARGUMENTS:
  --output-prefix STR       Prefix for new columns (default: "ChEMBL")
  --parallel N              Max concurrent requests (default: 12)
  --timeout SEC             Request timeout (default: 12)
  --retries N               Retry attempts per target (default: 2)
  --delay SEC               Delay after each request (default: 0)
  --cache FILE              JSON cache file (default: chembl_cache.json)
  --user-agent STRING       Custom User-Agent header
  --refresh                 Ignore cache and refetch all targets
  --verbose                 Print progress details

OUTPUT COLUMNS ADDED:
  <prefix>_HGNC
  <prefix>_UniProt
  <prefix>_Ensembl
  <prefix>_Pfam
  <prefix>_InterPro

EXAMPLE:
  python chembl_target_scraper.py \
      -i targets.tsv \
      -o targets_with_genes.tsv \
      -c "Target ChEMBL ID" \
      --parallel 16 \
      --refresh \
      --verbose
```
## 1-IV Run command

I ran it like following.

*Note that I always start by removing cache as otherwise it does not really update the whole result*

```bash
rm chembl_cache.json
python chembl_target_scraper.py \
  -i drug_targets.tsv \
  -o drug_targets_with_genes.tsv \
  -c "Target ChEMBL ID" \
  --output-prefix "ChEMBL" \
  --parallel 16 \
  --timeout 12 \
  --retries 3 \
  --verbose
```
## 1-V Select columns

Then I selected only the needed columns for claity

```bash
less drug_targets_with_genes.tsv | cut -f 1,2,3,4,5,11,12,13,21 > final_drug_info.tsv
```





