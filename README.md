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

# Run command

```
python chembl_target_genes_from_names.py \
  --input target_list_for_genes.xlsx \
  --sheet Sheet1 \
  --name-col 'Target Name' \
  --out genes_by_target.xlsx
```
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





