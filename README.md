# target_genes_for_drug_targets

# Introduction

I downloaded a list of drug targets for humans from 'https://www.ebi.ac.uk/chembl/explore/drug_mechanisms/STATE_ID:vZE8QPzTY1YchWXcushDsw%3D%3D' and extracted a unique list of target names (like 'Muscle-type nicotinic acetylcholine receptor'). I am using following script to extract genes that are targeted by these from ChEMBL.

But first I installed needed package inside my environment

# Preparation
```
module load gcc arrow
module load python
python -m venv ~/envs/scanpy
source ~/envs/scanpy/bin/activate
```

```
pip install pandas openpyxl requests
```

copied the excel file that looks like following to the same directory

```excel
Target Name
Muscle-type nicotinic acetylcholine receptor
Tyrosine-protein kinase ABL1
Beta-1 adrenergic receptor
Heparan sulfate
Sphingosine 1-phosphate receptor 1
Monoamine transporter
<img width="129" height="141" alt="image" src="https://github.com/user-attachments/assets/7149bd4d-4860-4840-bcc9-aefcf4a683c5" />
```
and then ran the script


# CLI help


# chembl_target_genes_from_names.py — CLI Help

Map free‑text **target names** (e.g., “Muscle‑type nicotinic acetylcholine receptor”) to curated **ChEMBL targets** and extract **HGNC gene symbols** + **UniProt accessions** using the official ChEMBL REST API (`/target`, `/target_component`).

---

## 🔧 Installation

```bash
pip install pandas openpyxl requests
```

Usage
```
python chembl_target_genes_from_names.py \
  --input <file.xlsx> \
  --name-col <column> \
  --out <output.xlsx> \
  [--sheet <sheet_name>] \
  [--allow-nonhuman] \
  [--requests-per-sec <N>]
```
Required arguments
```
--input FILE
Path to input Excel (.xlsx) containing target names.


--name-col COLUMN
Column name holding the target names.


--out FILE
Output Excel (.xlsx) to write results.


Optional arguments


--sheet NAME
Worksheet name (default: first sheet).


--allow-nonhuman
Include non‑human targets (default: human‑only).


--requests-per-sec N
Throttle API calls (default: 5 requests/sec).


-h, --help
Show help and exit.
```
Example
```

python chembl_target_genes_from_names.py \
  --input target_list_for_genes.xlsx \
  --sheet Sheet1 \
  --name-col "Target Name" \
  --out genes_by_target.xlsx
```
 Output
```
The script produces an Excel file with one row per matched target–gene pair:

input_target_name — the original name from your file
matched_target_pref_name — curated ChEMBL target preferred name
target_chembl_id — the ChEMBL target ID
target_type — e.g., SINGLE PROTEIN, PROTEIN COMPLEX, FAMILY, NON‑PROTEIN
organism — typically Homo sapiens (unless --allow-nonhuman)
gene_symbol — HGNC gene symbol (multiple rows for complexes/subunits)
uniprot_accession — UniProt accession for the component
evidence — direct API URLs used for the mapping


Notes

Protein complexes (e.g., nicotinic receptors, integrins) return multiple gene symbols—one row per subunit.
Families/selectivity groups and non‑protein targets (e.g., DNA, heparan sulfate) will have no single gene symbol; you can filter by target_type after export.
```
*Politeness & Reliability*

The script rate‑limits requests (default 5 req/s) and adds simple retry/backoff for transient HTTP errors.
Increase/decrease with --requests-per-sec depending on your network and API responsiveness.

# Run Command

```
python chembl_target_genes_from_names.py \
  --input target_list_for_genes.xlsx \
  --sheet Sheet1 \
  --name-col "Target Name" \
  --out genes_by_target.xlsx \
  --requests-per-sec 50
```

# Script 

chembl_target_genes_from_names.py

```py

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Map free-text target names to ChEMBL target gene symbols (HGNC) and UniProt accessions.

Workflow per input name:
  1) Search ChEMBL targets by curated preferred name (case-insensitive) using /chembl/api/data/target
  2) Fallback: broader full-text search with 'search' if needed
  3) Retrieve components of each matched target via /chembl/api/data/target_component
  4) Extract gene symbols (from component_synonyms) and UniProt accessions
  5) Output an Excel with one row per (input_name × matched_target × gene)

ChEMBL REST API docs & explorer:
- https://www.ebi.ac.uk/chembl/api/data/docs (live resource list and filter syntax)  [used: target, target_component]  # [1](https://chembl.blogspot.com/2023/03/what-is-max-phase-in-chembl.html)
- https://chembl.gitbook.io/chembl-interface-documentation/web-services/chembl-data-web-services  # [2](https://www.ebi.ac.uk/chembl/explore/target/CHEMBL3350223)
"""

import argparse
import sys
import time
from typing import Dict, List, Optional, Tuple

import pandas as pd
import requests

# ----------------------------
# Config
# ----------------------------
BASE = "https://www.ebi.ac.uk/chembl/api/data"
REQS_PER_SEC = 5.0          # polite default; adjust with --requests-per-sec
TIMEOUT = 30                # seconds
MAX_PAGE = 1000             # API pagination page size
MAX_RETRIES = 5             # transient HTTP retries
BACKOFF_SECS = 0.8          # base backoff

# HTTP session
session = requests.Session()
session.headers.update({"User-Agent": "chembl-target-gene-mapper/1.1 (+research use)"})


def _throttle():
    time.sleep(1.0 / max(REQS_PER_SEC, 0.1))


def _request_with_retry(url: str, params: Dict = None) -> requests.Response:
    """GET with basic retry/backoff for transient errors."""
    params = params or {}
    attempt = 0
    while True:
        _throttle()
        try:
            r = session.get(url, params=params, timeout=TIMEOUT)
            if r.status_code in (429, 500, 502, 503, 504):
                # backoff and retry
                attempt += 1
                if attempt > MAX_RETRIES:
                    r.raise_for_status()
                time.sleep(BACKOFF_SECS * (2 ** (attempt - 1)))
                continue
            r.raise_for_status()
            return r
        except requests.RequestException as e:
            attempt += 1
            if attempt > MAX_RETRIES:
                raise
            time.sleep(BACKOFF_SECS * (2 ** (attempt - 1)))


def paged_get(url: str, params: Dict = None) -> List[Dict]:
    """
    GET with ChEMBL REST-style pagination (limit/offset).
    Tries to infer the list payload key (e.g., 'targets', 'target_components').
    """
    out = []
    offset = 0
    while True:
        q = dict(params or {})
        q.update({"limit": MAX_PAGE, "offset": offset})
        r = _request_with_retry(url, q)
        data = r.json()

        # Determine list payload
        payload = None
        for key in ("targets", "target_components", "mechanisms", "molecules", "activities"):
            if isinstance(data, dict) and key in data:
                payload = data[key]
                break

        if payload is None:
            if isinstance(data, list):
                payload = data
            elif isinstance(data, dict):
                for v in data.values():
                    if isinstance(v, list):
                        payload = v
                        break

        if not payload:
            break

        out.extend(payload)
        if len(payload) < MAX_PAGE:
            break
        offset += MAX_PAGE

    return out


def search_targets_by_pref_name_icontains(name: str) -> List[Dict]:
    """Primary lookup: curated preferred name (case-insensitive contains)."""
    url = f"{BASE}/target.json"
    params = {"pref_name__icontains": name}
    return paged_get(url, params)


def search_targets_by_fulltext(name: str) -> List[Dict]:
    """Fallback: broader full-text search if pref_name query returned no matches."""
    url = f"{BASE}/target.json"
    params = {"search": name}
    return paged_get(url, params)


def get_target_components(target_chembl_id: str) -> List[Dict]:
    """Retrieve UniProt + synonyms (incl. HGNC gene symbols) for a target."""
    url = f"{BASE}/target_component.json"
    params = {"target_chembl_id": target_chembl_id}
    return paged_get(url, params)


def pick_best_targets(
    candidates: List[Dict],
    query_name: str,
    restrict_human: bool = True
) -> List[Dict]:
    """
    Heuristic ranking/filtering:
      1) Keep Homo sapiens when possible (unless --allow-nonhuman)
      2) Exact (case-insensitive) match on pref_name preferred
      3) Otherwise, return all remaining candidates
    """
    if not candidates:
        return []

    if restrict_human:
        # Keep human (and also allow None organism for higher-level types)
        c2 = [t for t in candidates if (t.get("organism") == "Homo sapiens" or t.get("organism") is None)]
        if c2:
            candidates = c2

    # Exact match on pref_name (case-insensitive)
    exact = [t for t in candidates if str(t.get("pref_name", "")).strip().lower() == query_name.strip().lower()]
    if exact:
        return exact

    return candidates


def extract_gene_symbols_from_components(components: List[Dict]) -> List[Tuple[Optional[str], Optional[str]]]:
    """
    From target components, extract [(gene_symbol, uniprot_accession), ...]
    Gene symbol is read from component_synonyms with syn_type containing 'GENE_SYMBOL'/'HGNC SYMBOL'.
    """
    out: List[Tuple[Optional[str], Optional[str]]] = []
    for comp in components:
        acc = comp.get("accession")  # UniProt
        syns = comp.get("component_synonyms") or []
        gene_syms = []

        for s in syns:
            syn = s.get("component_synonym")
            typ = (s.get("syn_type") or "").upper()
            if "GENE_SYMBOL" in typ or "HGNC" in typ:
                if syn:
                    gene_syms.append(syn)

        # conservative fallback if nothing labeled as gene symbol
        if not gene_syms:
            for s in syns:
                syn = s.get("component_synonym", "")
                if syn.isupper() and 2 <= len(syn) <= 12 and syn.replace("-", "").isalpha():
                    gene_syms.append(syn)

        if not gene_syms:
            gene_syms = [None]

        for g in gene_syms:
            out.append((g, acc))

    # De-duplicate
    seen = set()
    uniq: List[Tuple[Optional[str], Optional[str]]] = []
    for g, a in out:
        key = (g or "", a or "")
        if key not in seen:
            uniq.append((g, a))
            seen.add(key)
    return uniq


def process_targets(
    df: pd.DataFrame,
    name_col: str,
    restrict_human: bool = True
) -> pd.DataFrame:
    records = []

    for idx, row in df.iterrows():
        raw_name = str(row[name_col]).strip()
        if not raw_name or raw_name.lower() in {"nan", "none"}:
            continue

        # 1) narrow search on curated preferred name
        candidates = search_targets_by_pref_name_icontains(raw_name)
        # 2) fallback to a broader full-text search
        if not candidates:
            candidates = search_targets_by_fulltext(raw_name)

        best = pick_best_targets(candidates, raw_name, restrict_human=restrict_human)

        if not best:
            records.append({
                "input_target_name": raw_name,
                "matched_target_pref_name": None,
                "target_chembl_id": None,
                "target_type": None,
                "organism": None,
                "gene_symbol": None,
                "uniprot_accession": None,
                "evidence": None
            })
            continue

        for t in best:
            tchembl = t.get("target_chembl_id") or t.get("chembl_id")
            pref = t.get("pref_name")
            ttype = t.get("target_type")
            org = t.get("organism")

            comps = get_target_components(tchembl) if tchembl else []
            genes = extract_gene_symbols_from_components(comps) if comps else [(None, None)]

            # Evidence URLs for traceability
            target_url = f"{BASE}/target.json?target_chembl_id={tchembl}" if tchembl else None
            comp_url = f"{BASE}/target_component.json?target_chembl_id={tchembl}" if tchembl else None
            ev = ";".join([u for u in [target_url, comp_url] if u])

            for gene_symbol, acc in genes:
                records.append({
                    "input_target_name": raw_name,
                    "matched_target_pref_name": pref,
                    "target_chembl_id": tchembl,
                    "target_type": ttype,
                    "organism": org,
                    "gene_symbol": gene_symbol,
                    "uniprot_accession": acc,
                    "evidence": ev
                })

    out_df = pd.DataFrame.from_records(records)
    # Sort for readability
    sort_cols = ["input_target_name", "matched_target_pref_name", "target_chembl_id", "gene_symbol"]
    have_cols = [c for c in sort_cols if c in out_df.columns]
    if have_cols:
        out_df = out_df.sort_values(have_cols, kind="stable").reset_index(drop=True)
    return out_df


def main():
    global REQS_PER_SEC  # <-- Fix: declare global at the top before any references

    p = argparse.ArgumentParser(description="Map target names to gene symbols via ChEMBL.")
    p.add_argument("--input", required=True, help="Path to input Excel (.xlsx)")
    p.add_argument("--sheet", default=None, help="Worksheet name (default: first sheet)")
    p.add_argument("--name-col", required=True, help="Column name containing target names")
    p.add_argument("--out", required=True, help="Output Excel filename (.xlsx)")
    p.add_argument("--allow-nonhuman", action="store_true", help="Do not restrict to Homo sapiens")
    p.add_argument("--requests-per-sec", type=float, default=REQS_PER_SEC,
                   help="Throttle rate for API calls (default 5 req/s)")
    args = p.parse_args()

    # Apply user-specified rate
    if args.requests_per_sec and args.requests_per_sec > 0:
        REQS_PER_SEC = args.requests_per_sec

    # Read Excel
    try:
        df = pd.read_excel(args.input, sheet_name=args.sheet, engine="openpyxl")
    except Exception as e:
        print(f"Failed to read Excel: {e}", file=sys.stderr)
        sys.exit(1)

    if args.name_col not in df.columns:
        print(f"Column '{args.name_col}' not found in {args.input}. Columns: {list(df.columns)}", file=sys.stderr)
        sys.exit(1)

    # Process
    out_df = process_targets(df, args.name_col, restrict_human=(not args.allow_nonhuman))

    # Write Excel
    try:
        with pd.ExcelWriter(args.out, engine="openpyxl") as xw:
            out_df.to_excel(xw, index=False, sheet_name="genes_by_target")
    except Exception as e:
        print(f"Failed to write output Excel: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"Done. Wrote {len(out_df)} rows to {args.out}")


if __name__ == "__main__":
    main()
```






