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

Then I selected only the needed columns for clarity

```bash
less drug_targets_with_genes.tsv | cut -f 1,2,3,4,5,11,12,13,21 > final_drug_info.tsv
```
# 2) Adding side effects to the database

## 2-I Introduction

In addition to the drug targets and affected genes above, I found the drug warnings table from 'https://www.ebi.ac.uk/chembl/explore/drug_warnings'.

Here I am going to combine those warnings to the same database

I started this by downloading the table and renaming it as 'warnings_table.tsv'

it looks like following

```tsv
"Molecule ChEMBL ID"    "Parent Molecule ChEMBL ID"     "Parent Molecule Name"  "Warning Type"  "Warning Class"       "Description"   "Country"       "First Withdrawn Year"  "EFO ID"        "EFO Term"            "EFO ID for Warning Class"      "References"
"CHEMBL16073"   "CHEMBL16073"   "PHENACETIN"    "Withdrawn"     "hepatotoxicity"        "nephropathy; liver and kidney toxicity; carcinogenicity; haemolytic anaemia and methaemoglobinaemia" "['Chile', 'New Zealand', 'United States', 'Nigeria', 'Oman', 'Bahrain', 'Philippines', 'Rwanda', 'Malaysia', 'Thailand', 'United Arab Emirates', 'Hong Kong', 'Ireland', 'Netherlands', 'Sweden', 'Greece', 'Bangladesh', 'Cyprus', 'Yemen', 'Mauritius', 'Turkey', 'Nepal', 'Finland', 'Panama', 'Denmark', 'Romania', 'Austria', 'Ethiopia', 'Japan', 'Israel', 'United Kingdom', 'India', 'Italy', 'Suriname', 'Brazil', 'Norway', 'Germany']"    "[1965, 1965, 1965, 1965, 1965, 1965, 1965, 1965, 1965, 1965, 1965, 1965, 1965, 1965, 1965, 1965, 1965, 1965, 1965, 1965, 1965, 1965, 1965, 1965, 1965, 1965, 1965, 1965, 1965, 1965, 1965, 1965, 1965, 1965, 1965, 1965, 1965]"    "EFO:0011052"   "hepatotoxicity"        "EFO:0011052"   "Type: WHO RefID: CL12.pdf URL: https://www.un.org/esa/coordination/CL12.pdf|Type: WHO RefID: CL12.pdf URL: https://www.un.org/esa/coordination/CL12.pdf"
"CHEMBL4594503" "CHEMBL4594503" "TALQUETAMAB"   "Black Box Warning"     "immune system toxicity"              "None"  "['United States']"     "[None]"        "None"  "None"  "EFO:0011053"   ""
"CHEMBL608"     "CHEMBL608"     "PROBUCOL"      "Withdrawn"     "metabolic toxicity"    "Reduction in serum HDL-C levels and partly because of possible QT interval prolongation or ventricular arrhythmias"        "['United States']"     "[1995]"        "EFO:0004612"   "high density lipoprotein cholesterol measurement"    "EFO:0011054"   "Type: USGPO RefID: FR-1995-10-20/html/95-26053.htm URL: https://www.govinfo.gov/content/pkg/FR-1995-10-20/html/95-26053.htm|Type: PubMed RefID: 19457483 URL: http://europepmc.org/article/MED/19457483"
```
I am merging essential data from this dataframe to my main dataframe with my script merge_tsv_by_keys.py from 

## 2-II Script

```url
https://github.com/TharinduTS/cell_type_enrichment_v2?tab=readme-ov-file#1-ii-merge-script
```

## 2-III Run command

```
python ./merge_tsv_by_keys.py \
  --left final_drug_info.tsv \
  --right warnings_table.tsv \
  --left-keys "Parent Molecule ChEMBL ID,Parent Molecule Name" \
  --right-keys "Parent Molecule ChEMBL ID,Parent Molecule Name" \
  --right-cols "Warning Type,Warning Class,Description,Country,First Withdrawn Year,EFO ID,EFO Term" \
  --out  drug_info_with_genes_and_warnings.tsv
```
### Remove rows without gene info

There were some rows that had drug mechanisms but no related genes. I am going to drop rows without any genetic info

This can be done with following drop_empty_rows_plus.py

#### Script 

```py

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Enhanced CLI utility to drop rows with empty/missing values in specified column(s), with:
  • Clear messages showing initial rows, dropped counts, and final rows
  • Audit file saving the rows that were removed
  • Inverse mode to keep ONLY rows where selected columns are empty
  • Optional pre- and post- filters by regex or list membership
  • Batch processing of multiple input files

Usage examples (single file):
  # Basic: drop rows where ChEMBL_HGNC is empty (ANY mode)
  python drop_empty_rows_plus.py -i input.csv -o output.csv -c ChEMBL_HGNC

  # Drop only if ALL specified columns are empty; also write audit
  python drop_empty_rows_plus.py -i data.tsv -o cleaned.tsv -c ChEMBL_HGNC "Target Name" --mode all --audit-out dropped.tsv

  # Keep ONLY rows where chosen columns are empty (inverse)
  python drop_empty_rows_plus.py -i data.xlsx -o kept.xlsx -c ChEMBL_HGNC --inverse

  # With pre-filter: keep rows where Target Name matches regex and Action Type is in list
  python drop_empty_rows_plus.py -i input.csv -o output.csv -c ChEMBL_HGNC \
    --pre-keep-regex "Target Name=nicotinic.*receptor" \
    --pre-keep-in "Action Type=ANTAGONIST,AGONIST" --regex-ignore-case

  # With post-filter: after dropping, keep rows where Country is in set
  python drop_empty_rows_plus.py -i input.csv -o output.csv -c ChEMBL_HGNC \
    --post-keep-in "Country=United States,Canada"

Batch examples:
  # Process many CSV files at once, write outputs to a directory
  python drop_empty_rows_plus.py -i data/*.csv -c ChEMBL_HGNC --out-dir out/ --suffix _clean

  # Also write audit files for each input
  python drop_empty_rows_plus.py -i data/*.csv -c ChEMBL_HGNC --out-dir out/ --suffix _clean \
    --audit-dir audits/ --audit-suffix _dropped
"""

import argparse
import os
import sys
import re
import glob
from typing import List, Tuple

import pandas as pd

# ---------------------------- Defaults ---------------------------------------
DEFAULT_NA_VALUES = [
    "", " ", "  ",  # blank/whitespace (we also strip and convert to NA)
    "na", "n/a", "na/na",
    "null", "none", "nil",
    "nan", "—", "-", "--"
]

# ------------------------- IO helpers ----------------------------------------

def detect_format(path: str) -> str:
    ext = os.path.splitext(path.lower())[1]
    if ext in [".csv", ".tsv", ".txt"]:
        return "delimited"
    if ext in [".xlsx", ".xls"]:
        return "excel"
    return "delimited"


def read_table(
    path: str,
    sep: str = None,
    sheet: str = None,
    na_values: List[str] = None,
) -> Tuple[pd.DataFrame, str]:
    fmt = detect_format(path)
    if na_values is None:
        na_values = DEFAULT_NA_VALUES

    if fmt == "excel":
        engine = "openpyxl" if path.lower().endswith(".xlsx") else "xlrd"
        df = pd.read_excel(path, sheet_name=sheet, engine=engine)
    else:
        # Use engine='python' and sep=None to auto-detect if sep not provided
        if sep is None:
            df = pd.read_csv(
                path,
                sep=None,
                engine="python",
                na_values=na_values,
                keep_default_na=True,
            )
        else:
            df = pd.read_csv(
                path,
                sep=sep,
                na_values=na_values,
                keep_default_na=True,
            )
    return df, fmt


def write_table(df: pd.DataFrame, path: str):
    os.makedirs(os.path.dirname(path), exist_ok=True) if os.path.dirname(path) else None
    ext = os.path.splitext(path.lower())[1]
    if ext in [".xlsx", ".xls"]:
        # Prefer writing .xlsx; pandas does not support writing legacy .xls without extra deps
        if ext == ".xls":
            # Force to .xlsx by changing extension
            path = os.path.splitext(path)[0] + ".xlsx"
        df.to_excel(path, index=False, engine="openpyxl")
    else:
        sep = "\t" if ext == ".tsv" else ","
        df.to_csv(path, index=False, sep=sep)


# --------------------- Column & cleaning helpers -----------------------------
def normalize_columns_for_ci(df: pd.DataFrame, target_cols: List[str], case_insensitive: bool) -> List[str]:
    if not case_insensitive:
        # Validate presence
        missing = [c for c in target_cols if c not in df.columns]
        if missing:
            raise KeyError(f"Columns not found: {missing}. Available: {list(df.columns)}")
        return target_cols

    def norm(s: str) -> str:
        return s.strip().lower()

    actual_map = {norm(c): c for c in df.columns}
    resolved, missing = [], []
    for col in target_cols:
        key = norm(col)
        if key in actual_map:
            resolved.append(actual_map[key])
        else:
            missing.append(col)
    if missing:
        raise KeyError(f"(CI) Columns not found: {missing}. Available: {list(df.columns)}")
    return resolved


def strip_whitespace_and_empty_to_na(df: pd.DataFrame) -> pd.DataFrame:
    obj_cols = df.select_dtypes(include=["object"]).columns
    if len(obj_cols) > 0:
        df[obj_cols] = df[obj_cols].apply(lambda col: col.astype("string").str.strip())
        df[obj_cols] = df[obj_cols].replace({"": pd.NA})
    return df


# ------------------------ Filtering helpers ----------------------------------
def _parse_key_to_list(item: str) -> Tuple[str, List[str]]:
    """Parse 'Column=val1,val2,val3' -> (Column, [val1, val2, val3])"""
    if "=" not in item:
        raise ValueError(f"Expected KEY=V1,V2,... format, got: {item}")
    key, values = item.split("=", 1)
    values = [v.strip() for v in values.split(",") if v.strip() != ""]
    return key.strip(), values


def _parse_key_to_regex(item: str) -> Tuple[str, str]:
    """Parse 'Column=pattern' -> (Column, pattern)"""
    if "=" not in item:
        raise ValueError(f"Expected KEY=REGEX format, got: {item}")
    key, pattern = item.split("=", 1)
    return key.strip(), pattern.strip()


def apply_filters(
    df: pd.DataFrame,
    keep_regex: List[str],
    drop_regex: List[str],
    keep_in: List[str],
    drop_in: List[str],
    regex_ignore_case: bool,
    ci_columns: bool,
) -> pd.DataFrame:
    """
    Apply keep/drop filters.
    - For keep filters: rows must satisfy ALL keep conditions (AND across conditions).
    - For drop filters: row is removed if it matches ANY drop condition (OR across conditions).
    """
    working = df

    # Column resolution helper
    def resolve_col(name: str) -> str:
        if not ci_columns:
            return name
        # case-insensitive lookup
        lower_map = {c.lower(): c for c in working.columns}
        key = name.lower()
        if key not in lower_map:
            raise KeyError(f"Filter column not found: {name}. Available: {list(working.columns)}")
        return lower_map[key]

    # Build mask for keep (start as all True, AND with each condition)
    keep_mask = pd.Series(True, index=working.index)
    flags = re.IGNORECASE if regex_ignore_case else 0

    for expr in keep_regex or []:
        col, pattern = _parse_key_to_regex(expr)
        col = resolve_col(col)
        # Ensure string dtype for regex checks; NaNs -> ""
        s = working[col].astype("string").fillna("")
        keep_mask &= s.str.contains(pattern, regex=True, flags=flags)

    for expr in keep_in or []:
        col, values = _parse_key_to_list(expr)
        col = resolve_col(col)
        keep_mask &= working[col].isin(values)

    if (keep_regex and len(keep_regex) > 0) or (keep_in and len(keep_in) > 0):
        working = working[keep_mask]

    # Build mask for drop (start as all False, OR with each condition)
    if (drop_regex and len(drop_regex) > 0) or (drop_in and len(drop_in) > 0):
        drop_mask = pd.Series(False, index=working.index)
        for expr in drop_regex or []:
            col, pattern = _parse_key_to_regex(expr)
            col = resolve_col(col)
            s = working[col].astype("string").fillna("")
            drop_mask |= s.str.contains(pattern, regex=True, flags=flags)
        for expr in drop_in or []:
            col, values = _parse_key_to_list(expr)
            col = resolve_col(col)
            drop_mask |= working[col].isin(values)
        working = working[~drop_mask]

    return working


# --------------------------- Core logic --------------------------------------
def compute_empty_masks(df: pd.DataFrame, columns: List[str], mode: str) -> Tuple[pd.Series, pd.DataFrame]:
    """Return (mask_empty_row, isna_by_col) where:
    - mask_empty_row: True for rows considered empty based on mode ('any'/'all')
    - isna_by_col: DataFrame of booleans for empties per selected column
    """
    isna = df[columns].isna()
    if mode == "any":
        mask_empty = isna.any(axis=1)
    else:
        mask_empty = isna.all(axis=1)
    return mask_empty, isna


# ------------------------------ CLI ------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(
        description="Drop (or keep only) rows with empty/missing values in specified column(s). Supports audit, filters, and batch.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("-i", "--input", nargs="+", required=True,
                   help="Input file(s). You can pass multiple paths or use shell globs like data/*.csv")
    p.add_argument("-o", "--output", default=None,
                   help="Output file (single-input mode only). For batch, use --out-dir instead.")
    p.add_argument("--out-dir", default=None, help="Directory to write outputs for batch mode.")
    p.add_argument("--suffix", default="_clean", help="Suffix to append before extension for batch outputs.")

    p.add_argument("-c", "--columns", nargs="+", required=True,
                   help="Column name(s) to check for empty values.")
    p.add_argument("--mode", choices=["any", "all"], default="any",
                   help="'any': row is considered empty if ANY chosen column is empty; 'all': only if ALL are empty.")
    p.add_argument("--inverse", action="store_true",
                   help="Keep ONLY rows that are considered empty per --mode (instead of dropping them).")

    p.add_argument("--sep", default=None,
                   help="Delimiter for delimited files (CSV/TSV). If omitted, autodetect is used.")
    p.add_argument("--sheet", default=None, help="Excel sheet name (when reading Excel).")
    p.add_argument("--ci", action="store_true", help="Case-insensitive column name matching.")

    p.add_argument("--extra-na", nargs="*", default=[],
                   help="Extra strings to treat as NA in addition to defaults.")

    # Filters (pre and post)
    p.add_argument("--pre-keep-regex", action="append", default=[],
                   help="Pre-filter: keep rows where COLUMN matches REGEX. Format: 'Column=pattern'. Repeatable.")
    p.add_argument("--pre-drop-regex", action="append", default=[],
                   help="Pre-filter: drop rows where COLUMN matches REGEX. Format: 'Column=pattern'. Repeatable.")
    p.add_argument("--pre-keep-in", action="append", default=[],
                   help="Pre-filter: keep rows where COLUMN is in list. Format: 'Column=v1,v2'. Repeatable.")
    p.add_argument("--pre-drop-in", action="append", default=[],
                   help="Pre-filter: drop rows where COLUMN is in list. Format: 'Column=v1,v2'. Repeatable.")

    p.add_argument("--post-keep-regex", action="append", default=[],
                   help="Post-filter: keep rows where COLUMN matches REGEX. Format: 'Column=pattern'. Repeatable.")
    p.add_argument("--post-drop-regex", action="append", default=[],
                   help="Post-filter: drop rows where COLUMN matches REGEX. Format: 'Column=pattern'. Repeatable.")
    p.add_argument("--post-keep-in", action="append", default=[],
                   help="Post-filter: keep rows where COLUMN is in list. Format: 'Column=v1,v2'. Repeatable.")
    p.add_argument("--post-drop-in", action="append", default=[],
                   help="Post-filter: drop rows where COLUMN is in list. Format: 'Column=v1,v2'. Repeatable.")

    p.add_argument("--regex-ignore-case", action="store_true",
                   help="Make regex filters case-insensitive.")

    # Audit outputs
    p.add_argument("--audit-out", default=None,
                   help="Path to write audit file (single-input mode). Contains rows removed from the main output.")
    p.add_argument("--audit-dir", default=None, help="Directory to write per-file audits in batch mode.")
    p.add_argument("--audit-suffix", default="_dropped", help="Suffix for per-file audit outputs in batch mode.")

    p.add_argument("--list-columns", action="store_true", default=False,
                   help="List columns from the first input file and exit.")

    return p.parse_args()


def build_output_path(input_path: str, out_dir: str, suffix: str) -> str:
    base = os.path.basename(input_path)
    stem, ext = os.path.splitext(base)
    # Preserve original extension; if '.xls', we'll convert to .xlsx in writer
    out_name = f"{stem}{suffix}{ext if ext else '.csv'}"
    return os.path.join(out_dir, out_name) if out_dir else out_name


def build_audit_path(input_path: str, audit_dir: str, audit_suffix: str) -> str:
    base = os.path.basename(input_path)
    stem, ext = os.path.splitext(base)
    out_name = f"{stem}{audit_suffix}{ext if ext else '.csv'}"
    return os.path.join(audit_dir, out_name) if audit_dir else out_name


# ------------------------------ Runner ---------------------------------------
def process_one(
    in_path: str,
    out_path: str,
    audit_path: str,
    args,
) -> Tuple[int, int, int]:
    # Read
    try:
        na_values = DEFAULT_NA_VALUES + [s for s in (args.extra_na or [])]
        df, _fmt = read_table(in_path, sep=args.sep, sheet=args.sheet, na_values=na_values)
    except Exception as e:
        print(f"[ERROR] Failed to read '{in_path}': {e}", file=sys.stderr)
        return (0, 0, -1)

    # Optionally list columns and exit early (single-file diagnostic)
    if getattr(args, 'list_columns', False):
        print(f"[COLUMNS] {in_path}: {list(df.columns)}")
        return (0, 0, 0)

    # Normalize target columns
    try:
        target_cols = normalize_columns_for_ci(df, args.columns, args.ci)
    except KeyError as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        return (0, 0, -2)

    # Strip whitespace to ensure blanks become NA
    df = strip_whitespace_and_empty_to_na(df)

    # Pre-filters
    n0 = len(df)
    if any([args.pre_keep_regex, args.pre_drop_regex, args.pre_keep_in, args.pre_drop_in]):
        df = apply_filters(
            df,
            keep_regex=args.pre_keep_regex,
            drop_regex=args.pre_drop_regex,
            keep_in=args.pre_keep_in,
            drop_in=args.pre_drop_in,
            regex_ignore_case=args.regex_ignore_case,
            ci_columns=args.ci,
        )
    n_after_pre = len(df)

    # Compute empty masks and split into kept/removed according to inverse flag
    mask_empty, isna_by_col = compute_empty_masks(df, target_cols, args.mode)

    # Create a helper series listing which of the target columns are empty per row
    empty_cols_list = isna_by_col.apply(lambda row: ",".join([c for c, empty in row.items() if empty]), axis=1)

    if args.inverse:
        # Keep only rows that are empty per mode
        df_kept = df[mask_empty].copy()
        df_removed = df[~mask_empty].copy()
        drop_label = "not_empty"
    else:
        # Drop rows that are empty per mode
        df_kept = df[~mask_empty].copy()
        df_removed = df[mask_empty].copy()
        drop_label = "empty"

    # Attach audit metadata
    if len(df_removed) > 0:
        df_removed = df_removed.assign(__drop_reason=drop_label)
        # Add which target columns were empty (based on original isna_by_col)
        df_removed["__empty_columns"] = empty_cols_list.loc[df_removed.index].values

    # Post-filters on the kept set
    if any([args.post_keep_regex, args.post_drop_regex, args.post_keep_in, args.post_drop_in]):
        df_kept = apply_filters(
            df_kept,
            keep_regex=args.post_keep_regex,
            drop_regex=args.post_drop_regex,
            keep_in=args.post_keep_in,
            drop_in=args.post_drop_in,
            regex_ignore_case=args.regex_ignore_case,
            ci_columns=args.ci,
        )

    n_final = len(df_kept)
    n_dropped = n0 - n_final  # relative to initial dataset

    # Write outputs
    try:
        write_table(df_kept, out_path)
    except Exception as e:
        print(f"[ERROR] Failed to write output '{out_path}': {e}", file=sys.stderr)
        return (n0, n_final, -3)

    if audit_path:
        try:
            write_table(df_removed, audit_path)
        except Exception as e:
            print(f"[ERROR] Failed to write audit '{audit_path}': {e}", file=sys.stderr)
            # do not hard fail the whole run

    # Messages (what you asked for)
    print(f"[OK] {os.path.basename(in_path)} | rows: start={n0}, after_pre={n_after_pre}, final={n_final}, dropped={n_dropped}")
    print(f"[OK] -> output: {out_path}" + (f" | audit: {audit_path}" if audit_path else ""))

    return (n0, n_final, n_dropped)


def main():
    args = parse_args()

    # Expand potential globs
    inputs: List[str] = []
    for pattern in args.input:
        matches = glob.glob(pattern)
        if len(matches) == 0 and os.path.isfile(pattern):
            matches = [pattern]
        inputs.extend(matches)
    # Deduplicate while preserving order
    seen = set()
    unique_inputs = []
    for pth in inputs:
        if pth not in seen:
            seen.add(pth)
            unique_inputs.append(pth)

    if not unique_inputs:
        print("[ERROR] No input files found.", file=sys.stderr)
        sys.exit(1)

    # Single vs batch output handling
    single_mode = (len(unique_inputs) == 1)

    if single_mode:
        if args.output is None and args.out_dir is None:
            print("[ERROR] Provide -o/--output for single file, or --out-dir.", file=sys.stderr)
            sys.exit(2)
        in_path = unique_inputs[0]
        out_path = args.output if args.output else build_output_path(in_path, args.out_dir, args.suffix)
        audit_path = args.audit_out if args.audit_out else (build_audit_path(in_path, args.audit_dir, args.audit_suffix) if args.audit_dir else None)
        n0, n_final, n_dropped = process_one(in_path, out_path, audit_path, args)
        if n_dropped < 0:
            sys.exit(3)
        print(f"[SUMMARY] processed=1 | total_start={n0} | total_final={n_final} | total_dropped={n_dropped}")
    else:
        # Batch mode
        if args.out_dir is None:
            print("[ERROR] Multiple inputs detected. Please provide --out-dir for outputs.", file=sys.stderr)
            sys.exit(4)
        os.makedirs(args.out_dir, exist_ok=True)
        if args.audit_dir:
            os.makedirs(args.audit_dir, exist_ok=True)

        total_start = 0
        total_final = 0
        total_dropped = 0
        processed = 0
        for in_path in unique_inputs:
            out_path = build_output_path(in_path, args.out_dir, args.suffix)
            audit_path = None
            if args.audit_dir:
                audit_path = build_audit_path(in_path, args.audit_dir, args.audit_suffix)
            n0, nfin, ndrop = process_one(in_path, out_path, audit_path, args)
            if ndrop == -1 or ndrop == -2 or ndrop == -3:
                # Error already logged; continue with others
                continue
            total_start += n0
            total_final += nfin
            total_dropped += ndrop
            processed += 1
        print(f"[SUMMARY] processed={processed} files | total_start={total_start} | total_final={total_final} | total_dropped={total_dropped}")


if __name__ == "__main__":
    main()
```

#### CLI help

```txt
# drop_empty_rows_plus.py — Command Line Usage

Clean tabular data by dropping (or keeping only) rows with empty values
in one or more specified columns. Supports audit logs, regex/list filters,
and batch processing across many files.

------------------------------------------------------------
USAGE
------------------------------------------------------------
  python drop_empty_rows_plus.py -i INPUT [...] -c COLS [...] [OPTIONS]

Single file:
  python drop_empty_rows_plus.py -i data.csv -o clean.csv -c ChEMBL_HGNC

Batch mode (multiple files or globs):
  python drop_empty_rows_plus.py -i data/*.csv -c ChEMBL_HGNC --out-dir out/

------------------------------------------------------------
REQUIRED ARGUMENTS
------------------------------------------------------------
-i, --input FILE(S)
    One or more input files. Supports paths and glob patterns.

-c, --columns COL [COL ...]
    Column(s) used to determine if a row is empty.

------------------------------------------------------------
OUTPUT OPTIONS
------------------------------------------------------------
-o, --output FILE
    Output file (single-input mode).

--out-dir DIR
    Directory for batch outputs. Required if multiple inputs.

--suffix TEXT
    Suffix added before extension for batch outputs (default: _clean).

------------------------------------------------------------
EMPTY-ROW HANDLING
------------------------------------------------------------
--mode {any, all}
    any  = row considered empty if ANY chosen column is empty (default)
    all  = row considered empty only if ALL chosen columns are empty

--inverse
    Keep ONLY rows considered empty per --mode (reverse behavior).

------------------------------------------------------------
AUDIT LOGGING
------------------------------------------------------------
--audit-out FILE
    Write dropped rows to this file (single-input mode).

--audit-dir DIR
    Directory to write per-file audits (batch mode).

--audit-suffix TEXT
    Suffix for audit files (default: _dropped).

------------------------------------------------------------
FILTERING (PRE and POST)
------------------------------------------------------------
--pre-keep-regex    "COL=REGEX"
--pre-drop-regex    "COL=REGEX"
--pre-keep-in       "COL=V1,V2,..."
--pre-drop-in       "COL=V1,V2,..."

--post-keep-regex   "COL=REGEX"
--post-drop-regex   "COL=REGEX"
--post-keep-in      "COL=V1,V2,..."
--post-drop-in      "COL=V1,V2,..."

--regex-ignore-case
    Make regex matching case-insensitive.

------------------------------------------------------------
DATA CLEANING & PARSING
------------------------------------------------------------
--sep SEP
    Input delimiter (CSV/TSV). Auto-detected if omitted.

--sheet NAME
    Sheet name when reading Excel files.

--ci
    Case-insensitive column name matching.

--extra-na STR [STR ...]
    Additional strings treated as missing values (e.g., UNKNOWN TBD).

------------------------------------------------------------
UTILITY
------------------------------------------------------------
--list-columns
    Print column names from the first input file and exit.

------------------------------------------------------------
EXAMPLES
------------------------------------------------------------
# Drop rows with empty ChEMBL_HGNC
python drop_empty_rows_plus.py -i input.csv -o clean.csv -c ChEMBL_HGNC

# Drop rows where ANY of the listed columns are empty
python drop_empty_rows_plus.py -i data.tsv -o clean.tsv \
  -c ChEMBL_HGNC "Target Name" "Action Type"

# Keep ONLY rows where ChEMBL_HGNC is empty
python drop_empty_rows_plus.py -i input.csv -o empty_only.csv -c ChEMBL_HGNC --inverse

# Batch process multiple files
python drop_empty_rows_plus.py -i data/*.csv -c ChEMBL_HGNC \
  --out-dir out/ --suffix _clean --audit-dir audits/

------------------------------------------------------------
```


#### Run command

```bash
python drop_empty_rows_plus.py -i drug_info_with_genes_and_warnings.tsv -o drug_info_with_genes_and_warnings_fully_cleaned.tsv -c ChEMBL_HGNC
```

# 3) Adding enrichment data - Cell type Enrichment

## 3-I Introduction

Then I wanted to add enrichment data here 

I am doing it with 'gene_enrichment_join.py'

## 3-II Script 

gene_enrichment_join.py
```py

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Join gene lists from Table 1 with per-gene metrics from Table 2 and emit an
enriched list like: <value><pair-sep><gene><gene-key-sep><key_combo>{item-sep}...

Usage (example):
    python gene_enrichment_join.py \
        --table1 chembl_table1.tsv \
        --table2 enrichment_table2.tsv \
        --output chembl_with_enrichment.tsv \
        --t1-col ChEMBL_HGNC \
        --t1-sep "," \
        --t2-gene-col "Gene name" \
        --t2-key-col "Cell type" \
        --t2-value-col "log2_enrichment_penalized" \
        --gene-key-sep "|" \
        --pair-sep ":" \
        --item-sep " & " \
        --newcol-name ChEMBL_HGNC_enrichment

Notes:
- If --t2-gene-col is omitted, the script tries to auto-detect from common names:
  ["Gene name", "Gene", "gene", "HGNC symbol", "HGNC", "symbol"]
- Add multiple --t2-key-col flags if you want to combine several keys:
  e.g., --t2-key-col "Cell type" --t2-key-col "Cell type group"
"""

import argparse
import logging
import sys
from typing import List, Dict, Tuple, Optional

import pandas as pd


def detect_t2_gene_col(df2: pd.DataFrame, user_choice: Optional[str]) -> str:
    if user_choice and user_choice in df2.columns:
        return user_choice
    candidates = ["Gene name", "Gene", "gene", "HGNC symbol", "HGNC", "symbol"]
    for c in candidates:
        if c in df2.columns:
            logging.info(f"Auto-detected Table 2 gene column: {c}")
            return c
    raise ValueError(
        "Could not detect the gene column in Table 2. "
        "Please provide --t2-gene-col explicitly."
    )


def build_gene_index(
    df2: pd.DataFrame,
    gene_col: str,
    key_cols: List[str],
    value_col: str,
    case_sensitive: bool,
    key_join_sep: str,
    dedupe: bool,
) -> Dict[str, List[Tuple[str, str]]]:
    """
    Returns a dict:
        gene_norm -> list of (key_combo_str, value_str)
    """
    for c in [gene_col, value_col] + key_cols:
        if c not in df2.columns:
            raise KeyError(f"Column '{c}' not found in Table 2.")

    # Prepare normalized gene column
    gene_series = df2[gene_col].astype(str)
    if not case_sensitive:
        gene_series = gene_series.str.upper()

    # Build key combo
    key_combo = df2[key_cols].astype(str).apply(
        lambda r: key_join_sep.join(r.tolist()), axis=1
    )

    # Format value as string (preserve as-is)
    value_series = df2[value_col].apply(lambda x: "" if pd.isna(x) else str(x))

    tmp = pd.DataFrame({
        "_gene_norm": gene_series,
        "_key": key_combo,
        "_val": value_series
    })

    gene_map: Dict[str, List[Tuple[str, str]]] = {}

    if dedupe:
        # Drop exact duplicates for (gene, key, value)
        tmp = tmp.drop_duplicates(subset=["_gene_norm", "_key", "_val"])

    for _, row in tmp.iterrows():
        g = row["_gene_norm"]
        k = row["_key"]
        v = row["_val"]
        gene_map.setdefault(g, []).append((k, v))
    return gene_map


def parse_gene_list(
    cell_value: object,
    sep: str,
    null_tokens: List[str],
    case_sensitive: bool
) -> List[Tuple[str, str]]:
    """
    Returns list of tuples: (gene_original_token, gene_norm_for_lookup)
    """
    if pd.isna(cell_value):
        return []
    text = str(cell_value)
    if text.strip() == "":
        return []
    parts = [p.strip() for p in text.split(sep)]
    out = []
    null_set = {t.strip() for t in null_tokens}
    for token in parts:
        if token == "" or token in null_set:
            continue
        gene_norm = token if case_sensitive else token.upper()
        out.append((token, gene_norm))
    return out


def enrich_row(
    gene_tokens: List[Tuple[str, str]],
    gene_map: Dict[str, List[Tuple[str, str]]],
    gene_key_sep: str,
    pair_sep: str,
    item_sep: str,
    keep_empty: bool,
) -> str:
    """
    Build the enriched string for a single Table 1 row in the new order:
        <value>{pair_sep}<gene>{gene_key_sep}<key_combo>{item_sep}...
    """
    items: List[str] = []
    for gene_original, gene_norm in gene_tokens:
        hits = gene_map.get(gene_norm, [])
        if not hits:
            if keep_empty:
                # Keep the gene token without any value (makes presence explicit)
                items.append(gene_original)
            continue
        for key_combo, val in hits:
            # CHANGED ORDER: value first, then gene + key combo
            items.append(f"{val}{pair_sep}{gene_original}{gene_key_sep}{key_combo}")
    return item_sep.join(items)


def main():
    parser = argparse.ArgumentParser(
        description="Expand genes in Table 1 by joining to Table 2 and attaching values (value first)."
    )
    parser.add_argument("--table1", required=True, help="Path to Table 1 (TSV).")
    parser.add_argument("--table2", required=True, help="Path to Table 2 (TSV).")
    parser.add_argument("--output", required=True, help="Output TSV file path.")

    parser.add_argument("--t1-col", default="ChEMBL_HGNC",
                        help="Column in Table 1 containing the gene list.")
    parser.add_argument("--t1-sep", default=",",
                        help="Separator used in the Table 1 gene list.")
    parser.add_argument("--t1-null-values", default="NA,N/A,No NA",
                        help="Comma-separated placeholders to ignore in Table 1 gene list.")

    parser.add_argument("--t2-gene-col", default=None,
                        help="Gene column in Table 2. Auto-detected if omitted.")
    parser.add_argument("--t2-key-col", action="append", required=True,
                        help="Repeat for each key column in Table 2 (e.g., 'Cell type').")
    parser.add_argument("--t2-value-col", required=True,
                        help="Value column from Table 2 to attach (e.g., 'log2_enrichment_penalized').")

    parser.add_argument("--gene-key-sep", default="|",
                        help="Separator between gene and key combo in output.")
    parser.add_argument("--key-join-sep", default="::",
                        help="Separator joining multiple Table 2 key columns.")
    parser.add_argument("--pair-sep", default="=",
                        help="Separator between value and (gene+key) in output.")
    parser.add_argument("--item-sep", default="; ",
                        help="Separator used between multiple items in the output cell.")

    parser.add_argument("--overwrite", action="store_true",
                        help="Overwrite the Table 1 gene column with enriched text.")
    parser.add_argument("--newcol-name", default=None,
                        help="Name of the new column (if not overwriting). Defaults to <t1-col>_enrichment.")

    parser.add_argument("--case-sensitive", action="store_true",
                        help="Use case-sensitive gene matching (default: case-insensitive).")
    parser.add_argument("--keep-empty", action="store_true",
                        help="If a gene has no matches, keep the gene token (without value).")
    parser.add_argument("--no-dedupe", dest="dedupe", action="store_false",
                        help="Do not deduplicate identical (gene,key,value) entries.")
    parser.set_defaults(dedupe=True)

    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO, format="%(levelname)s: %(message)s", stream=sys.stderr
    )

    # Load input
    try:
        df1 = pd.read_csv(args.table1, sep="\t", dtype=str, keep_default_na=False)
        df2 = pd.read_csv(args.table2, sep="\t", dtype=str, keep_default_na=False)
    except Exception as e:
        logging.error(f"Failed to read inputs: {e}")
        sys.exit(1)

    if args.t1_col not in df1.columns:
        logging.error(f"Column '{args.t1_col}' not found in Table 1.")
        sys.exit(1)

    try:
        t2_gene_col = detect_t2_gene_col(df2, args.t2_gene_col)
    except Exception as e:
        logging.error(str(e))
        sys.exit(1)

    # Build index from Table 2
    try:
        gene_map = build_gene_index(
            df2=df2,
            gene_col=t2_gene_col,
            key_cols=args.t2_key_col,
            value_col=args.t2_value_col,
            case_sensitive=args.case_sensitive,
            key_join_sep=args.key_join_sep,
            dedupe=args.dedupe,
        )
    except Exception as e:
        logging.error(f"Failed to build gene index: {e}")
        sys.exit(1)

    # Prepare null tokens for Table 1 gene list
    null_tokens = [t.strip() for t in args.t1_null_values.split(",") if t.strip()]

    # Compute enriched column
    enriched_col_name = (
        args.t1_col if args.overwrite else (args.newcol_name or f"{args.t1_col}_enrichment")
    )

    def process_cell(cell):
        gene_tokens = parse_gene_list(
            cell_value=cell,
            sep=args.t1_sep,
            null_tokens=null_tokens,
            case_sensitive=args.case_sensitive,
        )
        return enrich_row(
            gene_tokens=gene_tokens,
            gene_map=gene_map,
            gene_key_sep=args.gene_key_sep,
            pair_sep=args.pair_sep,
            item_sep=args.item_sep,
            keep_empty=args.keep_empty,
        )

    try:
        df1[enriched_col_name] = df1[args.t1_col].apply(process_cell)
    except Exception as e:
        logging.error(f"Failed while enriching Table 1: {e}")
        sys.exit(1)

    # Write output
    try:
        df1.to_csv(args.output, sep="\t", index=False)
        logging.info(f"Wrote output to: {args.output}")
    except Exception as e:
        logging.error(f"Failed to write output: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
```
## 3-III CLI help

```txt
Gene Enrichment Join — Command Line Help
=======================================

This script expands gene lists in a Table 1 TSV file by matching each gene
against Table 2 and attaching per-gene values such as enrichment scores,
combined with one or more key columns (e.g., cell type).

Basic Usage
-----------
    python gene_enrichment_join.py \
        --table1 <table1.tsv> \
        --table2 <table2.tsv> \
        --output <output.tsv> \
        --t1-col "ChEMBL_HGNC" \
        --t2-gene-col "Gene name" \
        --t2-key-col "Cell type" \
        --t2-value-col "log2_enrichment_penalized"

Required Arguments
------------------
--table1 PATH              Input Table 1 (TSV format)
--table2 PATH              Input Table 2 (TSV format)
--output PATH              Output TSV file path
--t2-key-col COL           Key column(s) from Table 2 (repeatable)
--t2-value-col COL         Value column from Table 2 to attach

Optional Arguments
------------------
--t1-col COL               Gene list column in Table 1 (default: ChEMBL_HGNC)
--t1-sep SEP               Separator for genes in Table 1 (default: ",")
--t1-null-values LIST      Comma-separated tokens to ignore (default: "NA,N/A,No NA")

--t2-gene-col COL          Gene column in Table 2 (auto-detected if omitted)

--gene-key-sep SEP         Separator between gene and key combo (default: "|")
--key-join-sep SEP         Separator joining multiple key columns (default: "::")
--pair-sep SEP             Separator between <gene|key> and value (default: "=")
--item-sep SEP             Separator between enriched items (default: "; ")

--overwrite                Replace original gene list column in Table 1
--newcol-name NAME         Create a new output column (default: <t1-col>_enrichment)

--case-sensitive           Match genes case-sensitively (default: off)
--keep-empty               Include genes with no matches as empty entries
--no-dedupe                Do not deduplicate repeated (gene, key, value) entries

Description
-----------
For each gene listed in Table 1, the script finds all matching rows in Table 2,
combines the gene with one or more key columns, attaches the chosen value
column, and outputs a joined string such as:

    GENE|CellType=Value; GENE|CellType2=Value2

The result is written either to a new column or overwrites the existing gene
column as specified.
```

## 3-IV Run command

```bash
python gene_enrichment_join.py \
  --table1 drug_info_with_genes_and_warnings_fully_cleaned.tsv \
  --table2 Final_data_with_tissue_expression_data.tsv \
  --output chembl_with_enrichment.tsv \
  --t1-col "ChEMBL_HGNC" \
  --t1-sep "," \
  --t2-gene-col "Gene name" \
  --t2-key-col "Cell type" \
  --t2-value-col "log2_enrichment_penalized" \
  --gene-key-sep "-" \
  --pair-sep ":" \
  --item-sep " & " \
  --newcol-name "ChEMBL_HGNC_enrichment"
```
## 3-V selecting maximum enrichment

Because I need a value to prioritize when I plot, I am selecting the maximum enrichment value for each row.

### Script

extract_extreme_enrichment.py
```py

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Extract the max/min/mean enrichment from a "<value><pair-sep><label>" list column.

NEW FORMAT (value first):
    "6.18:CHRNA1-myosatellite cells & 5.44:CHRNA1-myonuclei & ..."

- item_sep splits items (e.g., " & ")
- pair_sep splits value from label (from the LEFT now, since value comes first)

Example (your case):
  python extract_extreme_enrichment.py \
    --input chembl_with_enrichment.tsv \
    --output chembl_with_enrichment_scored.tsv \
    --source-col "ChEMBL_HGNC_enrichment" \
    --item-sep " & " \
    --gene-key-sep "-" \
    --pair-sep ":" \
    --agg max \
    --label-col "Cell_type_with_max_enrichment" \
    --value-col "Penalized_enrichment_value" \
    --html-unescape
"""

import argparse
import html
import math
import sys
from typing import List, Tuple, Optional

import pandas as pd


def parse_items(
    text: str,
    item_sep: str,
    pair_sep: str,
    html_unescape: bool,
) -> List[Tuple[str, float]]:
    """
    Parse a cell like:
        "6.18:CHRNA1-myosatellite cells & 5.44:CHRNA1-myonuclei"
    into list of (label_str, value_float).

    NEW ORDER: value first, then label.
    We split from the LEFT (split(pair_sep, 1)) to capture the value correctly.
    """
    if text is None or str(text).strip() == "":
        return []

    s = str(text)
    if html_unescape:
        s = html.unescape(s)

    items = []
    for raw in s.split(item_sep):
        token = raw.strip()
        if not token:
            continue
        if pair_sep not in token:
            # Skip malformed item
            continue

        # NEW: value first, then label
        val_str, label = token.split(pair_sep, 1)
        val_str = val_str.strip()
        label = label.strip()

        try:
            v = float(val_str)
            if math.isfinite(v):
                items.append((label, v))
        except Exception:
            # Skip if value can't be parsed
            continue
    return items


def aggregate_items(
    items: List[Tuple[str, float]],
    agg: str,
    mean_label_policy: str,
    joiner_for_join_policy: str = "; "
) -> Tuple[str, Optional[float]]:
    """
    Returns (label_for_output, value_for_output).

    - For 'max': label is the label with the maximum value.
    - For 'min': label is the label with the minimum value.
    - For 'mean':
        value is the mean of values; label depends on mean_label_policy:
          * 'empty' -> ""
          * 'max'   -> label of max value
          * 'min'   -> label of min value
          * 'join'  -> join all labels (may be long)
    """
    if not items:
        return ("", None)

    if agg == "max":
        label, v = max(items, key=lambda x: x[1])
        return (label, v)

    if agg == "min":
        label, v = min(items, key=lambda x: x[1])
        return (label, v)

    if agg == "mean":
        values = [v for _, v in items]
        mean_val = sum(values) / len(values) if values else None
        if mean_label_policy == "empty" or mean_val is None:
            return ("", mean_val)
        elif mean_label_policy == "max":
            label, _ = max(items, key=lambda x: x[1])
            return (label, mean_val)
        elif mean_label_policy == "min":
            label, _ = min(items, key=lambda x: x[1])
            return (label, mean_val)
        elif mean_label_policy == "join":
            labels = [lbl for lbl, _ in items]
            return (joiner_for_join_policy.join(labels), mean_val)
        else:
            return ("", mean_val)

    # Unknown agg -> default to max
    label, v = max(items, key=lambda x: x[1])
    return (label, v)


def main():
    ap = argparse.ArgumentParser(
        description="Extract max/min/mean from a '<value><pair-sep><label>' list column and write two new columns."
    )
    ap.add_argument("--input", required=True, help="Input TSV file path.")
    ap.add_argument("--output", required=True, help="Output TSV file path.")
    ap.add_argument("--source-col", required=True,
                    help="Column containing items like 'Value:GENE-CellType' separated by a delimiter.")
    ap.add_argument("--item-sep", default="; ",
                    help="Separator between items (default: '; ')")
    ap.add_argument("--gene-key-sep", default="|",
                    help="(Informational) separator used between gene and key; not required for parsing.")
    ap.add_argument("--pair-sep", default="=",
                    help="Separator between value and label (default: '=')")
    ap.add_argument("--agg", choices=["max", "min", "mean"], default="max",
                    help="Aggregation to compute over values (default: max).")
    ap.add_argument("--label-col", default="Cell_type_with_max_enrichment",
                    help="Output column name for the selected label (default: 'Cell_type_with_max_enrichment').")
    ap.add_argument("--value-col", default="Penalized_enrichment_value",
                    help="Output column name for the numeric value (default: 'Penalized_enrichment_value').")
    ap.add_argument("--value-round", type=int, default=None,
                    help="Round numeric value to N decimals (optional).")
    ap.add_argument("--html-unescape", action="store_true",
                    help="HTML-unescape the source text (e.g., '&amp;' -> '&').")
    ap.add_argument("--mean-label-policy", choices=["empty", "max", "min", "join"],
                    default="max",
                    help="Label policy when --agg=mean (default: max).")
    ap.add_argument("--joiner-for-mean-label", default="; ",
                    help="Joiner used when --mean-label-policy=join (default: '; ')")

    args = ap.parse_args()

    # Read
    try:
        df = pd.read_csv(args.input, sep="\t", dtype=str, keep_default_na=False)
    except Exception as e:
        sys.stderr.write(f"ERROR reading input: {e}\n")
        sys.exit(1)

    if args.source_col not in df.columns:
        sys.stderr.write(f"ERROR: source column '{args.source_col}' not found.\n")
        sys.exit(1)

    # Compute per-row extraction
    labels_out = []
    values_out = []

    for _, row in df.iterrows():
        items = parse_items(
            text=row[args.source_col],
            item_sep=args.item_sep,
            pair_sep=args.pair_seP if hasattr(args, "pair_seP") else args.pair_sep,  # guard just in case
            html_unescape=args.html_unescape,
        )
        label, val = aggregate_items(
            items=items,
            agg=args.agg,
            mean_label_policy=args.mean_label_policy,
            joiner_for_join_policy=args.joiner_for_mean_label,
        )
        if val is not None and args.value_round is not None:
            try:
                val = round(float(val), args.value_round)
            except Exception:
                pass
        labels_out.append(label)
        values_out.append("" if val is None else val)

    df[args.label_col] = labels_out
    df[args.value_col] = values_out

    # Write
    try:
        df.to_csv(args.output, sep="\t", index=False)
    except Exception as e:
        sys.stderr.write(f"ERROR writing output: {e}\n")
        sys.exit(1)


if __name__ == "__main__":
    main()

```

### CLI help 

```txt
Extract Extreme Enrichment — Command Line Help
==============================================

This script processes a column containing items of the form:
    GENE-CellType:Value & GENE-OtherCellType:Value2 & ...
and extracts the maximum, minimum, or mean enrichment value per row.
Two new columns are added to the output table: one for the selected label
and one for the numeric value.

Basic Usage
-----------
    python extract_extreme_enrichment.py \
        --input chembl_with_enrichment.tsv \
        --output chembl_with_enrichment_scored.tsv \
        --source-col "ChEMBL_HGNC_enrichment" \
        --item-sep " & " \
        --gene-key-sep "-" \
        --pair-sep ":" \
        --agg max \
        --label-col "Cell_type_with_max_enrichment" \
        --value-col "Penalized_enrichment_value" \
        --html-unescape

Required Arguments
------------------
--input PATH                Input TSV file
--output PATH               Output TSV file
--source-col COL            Column containing "label:value" items
--t2? no but source-col? yes

Optional Arguments
------------------
--item-sep SEP              Separator between items (default: "; ")
--pair-sep SEP              Separator between label and value (default: "=")
--gene-key-sep SEP          Informational only; not required for parsing
--agg {max,min,mean}        Aggregation function (default: max)
--label-col NAME            Output column name for the selected label
--value-col NAME            Output column name for the numeric value
--value-round N             Round numeric value to N decimals
--html-unescape             Convert HTML escapes (e.g., "&amp;" → "&")

Mean-Aggregation Options
------------------------
--mean-label-policy {empty,max,min,join}
    Label choice when --agg=mean (default: max)
--joiner-for-mean-label STR
    Separator when joining labels under policy "join" (default: "; ")

Description
-----------
For each row, the script parses entries like:
    GENE-CellType:6.18 & GENE-OtherType:5.44
and returns the label/value pair with the max, min, or mean enrichment.
Output is written to two new columns in the TSV file.
```
### Run command

```bash
python extract_extreme_enrichment.py \
  --input chembl_with_enrichment.tsv \
  --output chembl_with_enrichment_scored.tsv \
  --source-col "ChEMBL_HGNC_enrichment" \
  --item-sep " & " \
  --gene-key-sep "-" \
  --pair-sep ":" \
  --agg max \
  --label-col "Cell_type_with_max_enrichment" \
  --value-col "Penalized_enrichment_value" \
  --html-unescape \
  --value-round 6
```


# 4) Fill any missing values


Even from the filtered data, only some drugs have warnings. Therefore I am filling the rows with no warnings with custom values for ease of data handling with fill_empty_values.py

## 4-I Script

fill_empty_values.py
```
#!/usr/bin/env python3
"""
fill_empty_values.py — Fill empty/NaN cells per-column with provided defaults.

Examples:
  python fill_empty_values.py \
    -i chembl_with_enrichment_scored.tsv \
    -o chembl_with_enrichment_scored_empties_filled.tsv \
    --fill "Parent Molecule ChEMBL ID=none" \
    --fill "Parent Molecule Name=No_name" \
    --fill "Parent Molecule Type=none" \
    --fill "Target ChEMBL ID=none" \
    --fill "Target Name=none" \
    --fill "Action Type=none" \
    --fill "ChEMBL_HGNC=none" \
    --fill "Warning Type=No Warnings found" \
    --fill "Warning Class=none" \
    --fill "Description=none" \
    --fill "Country=none" \
    --fill "First Withdrawn Year=none" \
    --fill "EFO ID=none" \
    --fill "EFO Term=none" \
    --fill "ChEMBL_HGNC_enrichment=none" \
    --fill "Cell_type_with_max_enrichment=none" \
    --fill "Max Phase=none" \
    --fill "First Approval=none" \
    --fill "Penalized_enrichment_value=none"
"""

import argparse
import sys
import ast
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Fill empty/NaN cells with per-column defaults safely (dtype-aware).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("-i", "--input", required=True, help="Path to TSV/CSV input file.")
    p.add_argument("-o", "--output", required=True, help="Path to output file.")
    p.add_argument(
        "--sep", default="\t", help=r"Field separator (e.g., '\t' for TSV, ',' for CSV)."
    )
    p.add_argument(
        "--encoding", default="utf-8", help="Input/output file encoding."
    )
    p.add_argument(
        "--na-values",
        nargs="*",
        default=["", " ", "  ", "NA", "N/A", "na", "n/a", "NULL", "null", "-", "--"],
        help="Tokens considered as NA when reading.",
    )
    p.add_argument(
        "--keep-default-na",
        action="store_true",
        help="Keep pandas default NA tokens in addition to --na-values.",
    )
    p.add_argument(
        "--no-strip-whitespace",
        dest="strip_whitespace",
        action="store_false",
        help="Do NOT strip whitespace before treating empties.",
    )
    p.set_defaults(strip_whitespace=True)

    p.add_argument(
        "--fill",
        action="append",
        default=[],
        help='Repeated flag: --fill "Column=Value". Example: --fill "Max Phase=0"',
    )
    p.add_argument(
        "--chunksize",
        type=int,
        default=None,
        help="Optional chunked processing for very large files.",
    )
    return p.parse_args()


def parse_fill_pairs(fill_args: List[str]) -> Dict[str, str]:
    """
    Parse ["Col=Val", "Other=Val2"] → { "Col": "Val", "Other": "Val2" }.
    Keeps values as strings; we infer numeric types later per-column.
    """
    mapping: Dict[str, str] = {}
    for item in fill_args:
        if "=" not in item:
            raise ValueError(f'Invalid --fill entry (needs "="): {item}')
        col, val = item.split("=", 1)
        col = col.strip()
        val = val.strip()
        if not col:
            raise ValueError(f"Empty column name in --fill: {item}")
        mapping[col] = val
    return mapping


def normalize_object_columns(df: pd.DataFrame, strip_whitespace: bool) -> pd.DataFrame:
    """
    Strip whitespace (optional) and convert empty strings to NaN for object/string-like columns.
    """
    def _clean(s: pd.Series) -> pd.Series:
        if s.dtype == "object" or pd.api.types.is_string_dtype(s):
            s = s.astype("string")
            if strip_whitespace:
                s = s.str.strip()
            s = s.replace("", pd.NA)
        return s
    return df.apply(_clean)


def looks_numeric(value: str) -> bool:
    """
    Heuristic: does a string look like an int/float/bool? (We handle bool separately.)
    """
    v = value.strip().lower()
    if v in {"true", "false"}:
        return True
    # Try literal_eval for ints/floats
    try:
        lit = ast.literal_eval(value)
        return isinstance(lit, (int, float, bool))
    except Exception:
        return False


def cast_series_and_value_for_fill(
    ser: pd.Series, fill_value_str: str
) -> Tuple[pd.Series, object]:
    """
    Decide whether to treat column as numeric or string for filling.
    - If fill value looks numeric and the column can be interpreted as numeric
      (ignoring current NAs), cast column to numeric and cast fill to numeric too.
    - Else, ensure column is string dtype and keep fill as a string.

    Returns: (possibly-casted series, casted fill value)
    """
    # Try to coerce series to numeric (preserving NA)
    ser_numeric = pd.to_numeric(ser, errors="coerce")
    can_be_numeric = not ser_numeric.dropna().empty  # at least one non-NA numeric

    # Try to interpret fill value as literal (int/float/bool) if it looks numeric
    fill_value_casted = fill_value_str
    if looks_numeric(fill_value_str):
        try:
            lit = ast.literal_eval(fill_value_str)
            fill_value_casted = lit
        except Exception:
            pass  # keep as string if literal_eval fails

    # Decide path
    if isinstance(fill_value_casted, (int, float, bool)) and can_be_numeric:
        # Prefer pandas nullable dtypes for integers to allow NA
        if isinstance(fill_value_casted, bool):
            # Boolean path (nullable boolean dtype)
            ser_out = ser_numeric.astype("Float64")  # intermediate to preserve NA
            # Convert to boolean only if values are 0/1/NA; otherwise keep as Float64
            uniq_non_na = pd.unique(ser_out.dropna())
            if set(uniq_non_na).issubset({0.0, 1.0}):
                ser_out = ser_out.astype("Int64").astype("boolean")
                fill_value_casted = bool(fill_value_casted)
                return ser_out, fill_value_casted
            else:
                # If not binary, fill as numeric (Float64)
                ser_out = ser_numeric.astype("Float64")
                fill_value_casted = 1.0 if fill_value_casted else 0.0
                return ser_out, fill_value_casted

        # Distinguish int-like vs float-like based on existing non-NA data
        # If all non-NA are integers, use nullable Int64; else Float64
        non_na = ser_numeric.dropna()
        if not non_na.empty and np.all(np.equal(np.mod(non_na, 1), 0)):
            ser_out = ser_numeric.astype("Int64")
            if isinstance(fill_value_casted, float) and fill_value_casted.is_integer():
                fill_value_casted = int(fill_value_casted)
            if isinstance(fill_value_casted, int):
                return ser_out, fill_value_casted
            # If fill is float but column looks int-ish, try to cast fill to int
            try:
                fill_value_casted = int(fill_value_casted)
                return ser_out, fill_value_casted
            except Exception:
                # Fall back to float if we can't (rare)
                ser_out = ser_numeric.astype("Float64")
                fill_value_casted = float(fill_value_casted)
                return ser_out, fill_value_casted
        else:
            ser_out = ser_numeric.astype("Float64")
            fill_value_casted = float(fill_value_casted)
            return ser_out, fill_value_casted

    # String path: ensure string dtype
    ser_out = ser.astype("string")
    fill_value_casted = str(fill_value_casted)
    return ser_out, fill_value_casted


def apply_fills(df: pd.DataFrame, fills: Dict[str, str]) -> pd.DataFrame:
    """
    Apply per-column fills safely, converting dtypes as needed to avoid
    'incompatible dtype' warnings/errors.
    """
    for col, fill_str in fills.items():
        if col not in df.columns:
            raise KeyError(f"Column not found: {col}")

        # Treat empties as NA already; build the mask
        mask = df[col].isna()

        if not mask.any():
            continue  # nothing to fill

        # Decide dtype strategy and cast both series and fill value
        ser_casted, fill_casted = cast_series_and_value_for_fill(df[col], fill_str)

        # Assign the casted series back (to keep dtype), then fill only masked rows
        df[col] = ser_casted
        df.loc[mask, col] = fill_casted

    return df


def process_in_memory(
    path: str,
    sep: str,
    na_values: List[str],
    keep_default_na: bool,
    strip_whitespace: bool,
    fills: Dict[str, str],
    encoding: str,
) -> pd.DataFrame:
    df = pd.read_csv(
        path,
        sep=sep,
        dtype=str,  # read as strings first
        na_values=na_values,
        keep_default_na=keep_default_na,
        encoding=encoding,
    )
    df = normalize_object_columns(df, strip_whitespace)
    df = apply_fills(df, fills)
    return df


def process_in_chunks(
    path: str,
    sep: str,
    na_values: List[str],
    keep_default_na: bool,
    strip_whitespace: bool,
    fills: Dict[str, str],
    encoding: str,
    chunksize: int,
    output: str,
):
    first = True
    for chunk in pd.read_csv(
        path,
        sep=sep,
        dtype=str,
        na_values=na_values,
        keep_default_na=keep_default_na,
        encoding=encoding,
        chunksize=chunksize,
    ):
        chunk = normalize_object_columns(chunk, strip_whitespace)
        chunk = apply_fills(chunk, fills)
        chunk.to_csv(
            output, sep=sep, index=False, header=first, mode="w" if first else "a", encoding=encoding
        )
        first = False


def main():
    args = parse_args()
    fills = parse_fill_pairs(args.fill)

    try:
        if args.chunksize:
            process_in_chunks(
                path=args.input,
                sep=args.sep,
                na_values=args.na_values,
                keep_default_na=args.keep_default_na,
                strip_whitespace=args.strip_whitespace,
                fills=fills,
                encoding=args.encoding,
                chunksize=args.chunksize,
                output=args.output,
            )
        else:
            df = process_in_memory(
                path=args.input,
                sep=args.sep,
                na_values=args.na_values,
                keep_default_na=args.keep_default_na,
                strip_whitespace=args.strip_whitespace,
                fills=fills,
                encoding=args.encoding,
            )
            df.to_csv(args.output, sep=args.sep, index=False, encoding=args.encoding)

        print(f"Saved filled file to: {args.output}")
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()

```

#### CLI help

```txt
fill_empty_values.py — Command Line Tool
----------------------------------------

Fill empty values (NaN or blank/whitespace cells) in selected columns of CSV, TSV,
or Excel files. Supports custom fill values, multiple columns, case-insensitive
matching, and outputs a summary showing how many cells were filled.

Usage:
  python fill_empty_values.py -i INPUT -o OUTPUT [options]

Required Arguments:
  -i, --input FILE           Path to input file (.csv, .tsv, .xlsx, .xls)
  -o, --output FILE          Path to output file (.csv, .tsv, .xlsx)

Optional Arguments:
  --sheet NAME               Sheet name/index when reading Excel files
  -f, --fill COL=VALUE       Fill empty cells in a specific column.
                             Repeat for multiple columns.
                             Example: --fill "Warning Type=No Warnings found"
  -c, --columns COL [COL...] Specify one or more columns to fill using --value.
  -v, --value VALUE          Value applied to all columns listed with --columns
  --case-insensitive         Match column names ignoring case
  -h, --help                 Show this help message

Behavior:
  • Empty = NaN, None, or blank/whitespace-only string
  • Supports CSV, TSV, XLSX, XLS (read-only)
  • Writes CSV, TSV, or XLSX output
  • Prints summary:
        - Total number of rows
        - Number of cells filled per column

Examples:
  Fill one column:
    python fill_empty_values.py -i data.csv -o out.csv \
      --fill "Warning Type=No Warnings found"

  Fill two columns with their own values:
    python fill_empty_values.py -i drugs.xlsx -o out.xlsx \
      --fill "Warning Type=No Warnings found" \
      --fill "Description=Not Available"

  Fill multiple columns using the same value:
    python fill_empty_values.py -i data.tsv -o out.tsv \
      --columns "Warning Type" "Description" --value "N/A"

----------------------------------------
```

#### Run Command

```bash
python fill_empty_values.py \
  -i chembl_with_enrichment_scored.tsv \
  -o chembl_with_enrichment_scored_empties_filled.tsv \
  --fill "Parent Molecule ChEMBL ID=none" \
  --fill "Parent Molecule Name=No_name" \
  --fill "Parent Molecule Type=none" \
  --fill "Target ChEMBL ID=none" \
  --fill "Target Name=none" \
  --fill "Action Type=none" \
  --fill "ChEMBL_HGNC=none" \
  --fill "Warning Type=No Warnings found" \
  --fill "Warning Class=none" \
  --fill "Description=none" \
  --fill "Country=none" \
  --fill "First Withdrawn Year=none" \
  --fill "EFO ID=none" \
  --fill "EFO Term=none" \
  --fill "ChEMBL_HGNC_enrichment=none" \
  --fill "Cell_type_with_max_enrichment=none" \
  --fill "Max Phase=0" \ 
  --fill "First Approval=0" \
  --fill "Penalized_enrichment_value=0"
```

# 5) Making interactive plots

## Run command

```bash
python universal_plot_maker_plus.py \
  --file chembl_with_enrichment_scored_empties_filled.tsv \
  --out chembl_with_enrichment_scored.html \
  --plot-type bar \
  --x-choices "Parent Molecule ChEMBL ID | Parent Molecule Name | Target ChEMBL ID | Target Name" \
  --y-choices "Penalized_enrichment_value" \
  --default-x "Target Name" \
  --default-y "Penalized_enrichment_value" \
  --color-col "Max Phase" \
  --color-choices "Max Phase|Parent Molecule Name|Parent Molecule Type|First Approval|Target Name|Action Type|Warning Type|Warning Class" \
  --filter-cols "Max Phase|Parent Molecule Name|Parent Molecule Type|First Approval|Target Name|Action Type|Warning Type|Warning Class" \
  --search-cols "Max Phase|Parent Molecule Name|Parent Molecule Type|First Approval|Target Name|Action Type|Warning Type|Warning Class" \
  --details "Parent Molecule ChEMBL ID|Parent Molecule Name|Parent Molecule Type|Max Phase|First Approval|Target ChEMBL ID|Target Name|Action Type|ChEMBL_HGNC|Warning Type|Warning Class|Description|Country|First Withdrawn Year|EFO ID|EFO Term|Cell_type_with_max_enrichment|Penalized_enrichment_value" \
  --title "Drug targets V1" \
  --dup-policy overlay \
  --sort-primary "Max Phase" \
  --sort-primary-order desc \
  --sort-secondary "Target Name" \
  --sort-secondary-order asc \
  --initial-zoom 100 \
  --self-contained \
  --lang en \
  --pt-enable \
  --pt-col "ChEMBL_HGNC_enrichment" \
  --pt-title "Enrichment for cell types" \
  --pt-x-label "cell type" \
  --pt-y-label "log2 Enrichment Penalized" \
  --pt-color "#2a9d8f" \
  --pt-height 360 \
  --pt-width auto \
  --pt-rotate -35 \
  --pt-container-id "present-tissues-plot" \
  --pt-enable --pt-mode flow \
  --pt-anchor "#rowDetails" --pt-position after \
  --pt-offset-x -300 --pt-offset-y -10
```
