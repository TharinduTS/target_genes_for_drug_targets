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
python drop_empty_rows_plus.py -i drug_info_with_genes_and_warnings.tsv -o drug_info_with_genes_and_warnings_empties_dropped.tsv -c ChEMBL_HGNC
```
### Fill empty values

Even from the filtered data, only some drugs have warnings. Therefore I am filling the rows with no warnings with custom values for ease of data handling with fill_empty_values.py

#### Script

fill_empty_values.py
```
#!/usr/bin/env python3
"""
fill_empty_values.py

CLI utility to fill empty (NaN/blank/whitespace-only) cells in selected column(s)
with user-specified value(s) for CSV or Excel files.

Usage examples:

1) Single column → single value
   python fill_empty_values.py -i input.csv -o output.csv \
       --fill "Warning Type=No Warnings found"

2) Multiple columns → each with its own value (repeat --fill)
   python fill_empty_values.py -i drugs.xlsx -o drugs_filled.xlsx --sheet "Sheet1" \
       --fill "Warning Type=No Warnings found" \
       --fill "Description=No description provided"

3) Apply the same value to several columns
   python fill_empty_values.py -i input.csv -o output.csv \
       --columns "Warning Type" "Description" \
       --value "N/A"

Notes:
- Empty is defined as: NaN/None OR a string that is empty after stripping whitespace.
- Column name matching is exact by default; use --case-insensitive to match ignoring case.
- Supports .csv, .xlsx (and legacy .xls for reading).

Author: Generated by M365 Copilot
"""

import argparse
import os
import sys
from typing import Dict, List, Optional, Tuple
import pandas as pd

# Pandas engine hints per environment requirements
READ_EXCEL_ENGINES = {
    ".xlsx": "openpyxl",
    ".xls": "xlrd",
}


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def parse_fill_pairs(pairs: List[str]) -> Dict[str, str]:
    """Parse a list of strings like ["Col=Value", "Other=Val"] into a dict.
    Allows values to contain '=' by splitting only on the first '='.
    Strips surrounding quotes and whitespace.
    """
    mapping: Dict[str, str] = {}
    for p in pairs or []:
        if "=" not in p:
            raise ValueError(
                f"Invalid --fill entry '{p}'. Expected format COL=VALUE."
            )
        col, val = p.split("=", 1)
        col = col.strip().strip('"').strip("'")
        val = val.strip().strip('"').strip("'")
        if not col:
            raise ValueError("Column name in --fill cannot be empty.")
        mapping[col] = val
    return mapping


def normalize_columns(df: pd.DataFrame, case_insensitive: bool) -> Tuple[pd.DataFrame, Dict[str, str]]:
    """If case_insensitive, create a mapping from lower-cased name -> actual name.
    Returns the (unchanged) df and a lookup dict.
    """
    if not case_insensitive:
        return df, {c: c for c in df.columns}
    lookup = {c.lower(): c for c in df.columns}
    return df, lookup


def resolve_column_name(user_col: str, lookup: Dict[str, str], case_insensitive: bool) -> Optional[str]:
    if case_insensitive:
        return lookup.get(user_col.lower())
    return lookup.get(user_col)


def is_empty_series(s: pd.Series) -> pd.Series:
    """Return a Boolean mask where cells are considered empty (NaN/None/blank)."""
    if s.dtype == object:
        # For strings: empty if NaN OR strip()==""
        return s.isna() | s.astype(str).str.strip().eq("")
    else:
        # For non-strings: empty if NaN
        return s.isna()


def apply_fills(df: pd.DataFrame, fills: Dict[str, str], columns: List[str], value: Optional[str], case_insensitive: bool) -> pd.DataFrame:
    df = df.copy()

    # Build column lookup
    _, lookup = normalize_columns(df, case_insensitive)

    # Mode A: --fill pairs take precedence if provided
    if fills:
        for user_col, fill_val in fills.items():
            actual_col = resolve_column_name(user_col, lookup, case_insensitive)
            if actual_col is None:
                raise KeyError(f"Column '{user_col}' not found in input.")
            mask = is_empty_series(df[actual_col])
            df.loc[mask, actual_col] = fill_val

    # Mode B: --columns with a single --value
    if columns:
        if value is None:
            raise ValueError("When using --columns, you must also provide --value.")
        for user_col in columns:
            actual_col = resolve_column_name(user_col, lookup, case_insensitive)
            if actual_col is None:
                raise KeyError(f"Column '{user_col}' not found in input.")
            mask = is_empty_series(df[actual_col])
            df.loc[mask, actual_col] = value

    return df


def read_table(path: str, sheet: Optional[str]) -> pd.DataFrame:
    ext = os.path.splitext(path)[1].lower()
    if ext == ".csv":
        return pd.read_csv(path)
    elif ext in READ_EXCEL_ENGINES:
        engine = READ_EXCEL_ENGINES[ext]
        return pd.read_excel(path, sheet_name=sheet or 0, engine=engine)
    else:
        raise ValueError(f"Unsupported input extension '{ext}'. Use .csv, .xlsx, or .xls")


def write_table(df: pd.DataFrame, path: str, sheet: Optional[str]):
    ext = os.path.splitext(path)[1].lower()
    if ext == ".csv":
        df.to_csv(path, index=False)
    elif ext in READ_EXCEL_ENGINES:
        # For writing, pandas chooses engine automatically for xlsx; we can be explicit when possible
        if ext == ".xlsx":
            df.to_excel(path, index=False, sheet_name=sheet or "Sheet1", engine="openpyxl")
        elif ext == ".xls":
            # pandas no longer writes .xls in many environments; provide a helpful error
            raise ValueError(
                "Writing .xls is not supported. Please use an .xlsx output file."
            )
    else:
        raise ValueError(f"Unsupported output extension '{ext}'. Use .csv or .xlsx")


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Fill empty cells in selected columns with provided values (CSV/Excel)."
    )
    p.add_argument("-i", "--input", required=True, help="Path to input CSV/XLSX/XLS file.")
    p.add_argument("-o", "--output", required=True, help="Path to output CSV/XLSX file.")
    p.add_argument("--sheet", help="Excel sheet name/index (reading). If omitted, first sheet is used.")
    p.add_argument(
        "-f",
        "--fill",
        action="append",
        help=(
            "Column-to-value mapping in the form COL=VALUE. Repeatable. "
            "Example: --fill 'Warning Type=No Warnings found'"
        ),
    )
    p.add_argument(
        "-c",
        "--columns",
        nargs="*",
        help="Column(s) to fill (use with --value for a single value applied to all).",
    )
    p.add_argument(
        "-v",
        "--value",
        help="Value to use with --columns for filling empty cells.",
    )
    p.add_argument(
        "--case-insensitive",
        action="store_true",
        help="Match column names ignoring case.",
    )
    return p


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    try:
        fills = parse_fill_pairs(args.fill) if args.fill else {}
        if not fills and not args.columns:
            raise ValueError(
                "You must provide either --fill COL=VALUE (one or more) or --columns ... --value ..."
            )

        df = read_table(args.input, args.sheet)
        df_out = apply_fills(
            df,
            fills=fills,
            columns=args.columns or [],
            value=args.value,
            case_insensitive=args.case_insensitive,
        )
        write_table(df_out, args.output, args.sheet)
        print(f"✅ Wrote output to {args.output}")
        return 0
    except Exception as e:
        eprint(f"Error: {e}")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
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
  -i drug_info_with_genes_and_warnings_empties_dropped.tsv \
  -o drug_info_with_genes_and_warnings_fully_cleaned.tsv \
  --fill "Warning Type=No Warnings found" \
  --fill "Warning Class=NA" \
  --fill "Description=NA" \
  --fill "Country=NA" \
  --fill "First Withdrawn Year=NA" \
  --fill "EFO ID=NA" \
  --fill "EFO Term=No NA"
```
# 3) Adding enrichment data
