# Enzyme Variant & Structure Analysis for Butyrate Kinase (PDB: 1SAZ)

This analysis connects **3D structure**, **secondary-structure geometry**, and **sequence variation across species** for a butyrate kinase family, using PDB entry **1SAZ** as the structural template.

The workflow is implemented in three Python scripts:

1. `script1_loops.py` – identify loop regions in 1SAZ and summarize B-factors (flexibility).
2. `script2_secondary_distances.py` – compute end-to-end distances for helices and strands based on CA coordinates.
3. `script3.py` – group sequences by species using a UniProt-style table and a CLUSTAL multiple sequence alignment (MSA), then compute amino-acid usage per alignment position for each species group.

---

## Data & Files Used

* **PDB structure:** `1SAZ.pdb` (butyrate kinase).
* **Script 1:** `script1_loops.py`
* **Script 2:** `script2_secondary_distances.py`
* **Script 3:** `script3.py`
* **Sequence metadata table (TSV):** `ungrouped_table.tsv`

  * Contains at least: `Entry Name` (or `Entry`) and `Organism`.
* **Multiple sequence alignment:** `aligned.txt` (CLUSTAL format)
* **Outputs:**

  * `loops.txt` – loop summary and B-factors (Script 1).
  * `secondary_distances.txt` – helix/strand distances (Script 2).
  * `grouped_table.tsv` – input table plus a `Group` (species) column (Script 3).
  * `enzyme_analysis_grouped.txt` – amino-acid usage per position per species group (Script 3).

All scripts use **Python 3** and only the standard library.

---

## How to Run the Scripts

From the command line, the typical usage is:

```bash
# Script 1: Loop identification + B-factors
python script1_loops.py 1SAZ.pdb > loops.txt

# Script 2: Secondary-structure distances
python script2_secondary_distances.py 1SAZ.pdb > secondary_distances.txt

# Script 3: Species groups + amino-acid usage
python script3.py \
  --table ungrouped_table.tsv \
  --msa aligned.txt \
  --out_table grouped_table.tsv \
  > enzyme_analysis_grouped.txt
```

---

## Script 1 – Loop Identification and B-Factor Analysis (`script1_loops.py`)

### Objective

Script 1 uses the 1SAZ butyrate kinase structure to:

* Identify **loop regions** (residues not in helices or sheets).
* Summarize **average B-factors** for each loop as a proxy for flexibility.
* List which residues belong to each loop region.

### Step 1: Parse ATOM Records

The script first parses `ATOM` lines from the PDB:

```python
residue_atoms = parse_pdb_atoms(pdb_path)
```

For each residue, it stores:

* `resname` – three-letter residue name (e.g., LYS, GLY).
* `B_factors` – list of B-factors for all atoms in that residue.

Residues are keyed as `(chain_id, residue_number)`. This creates a dictionary that can be used later to compute per-residue and per-loop B-factors.

### Step 2: Identify Secondary Structure Regions

Next, the script parses `HELIX` and `SHEET` records:

```python
sec_residues = parse_secondary_structure_regions(pdb_path)
```

This step collects every residue index that belongs to either an alpha helix or a beta strand. The result is a set of residues considered **secondary structure**, which lets the script classify everything else as a loop.

### Step 3: Find Loop Segments

Loop segments are defined as **continuous stretches of residues that are not in helices or sheets**:

```python
loops = find_loops(residue_atoms, sec_residues)
```

For each chain:

1. Residue indices are sorted.
2. The script walks through residues and starts a new loop every time it encounters one or more consecutive residues that are **not** in `sec_residues`.

Conceptually, if residues 40–55 are helix, 56–60 are not in any helix/sheet, and 61–70 are sheet, then residues 56–60 form one loop.

### Step 4: Compute Average B-Factors for Each Loop

For each loop:

```python
for loop in loops:
    b_values = []
    residue_labels = []
    for (chain, res_seq) in loop:
        info = residue_atoms.get((chain, res_seq))
        ...
```

* All B-factors from residues in the loop are concatenated.
* The mean B-factor is computed as the simple average over all residues in the loop:

  `Average B-factor = (Σᵢ Bᵢ) / N`, where `N` is the number of residues in the loop.

$$
\text{Average B-factor} = \frac{\sum_i B_i}{N}
$$


* Each residue is labeled (e.g., `LYS45A`) for easy mapping back onto the structure.

### Step 5: Output and Interpretation

Example text output:

```text
Loop 1: Chain A, residues 45–52
  Average B-factor: 19.23 Å^2
  Residues:
    GLY45A, SER46A, ASP47A, ...
```

In practice, these loop regions can be visualized in PyMOL by highlighting the reported residues (Figure 1), allowing us to connect **high B-factors** (more flexible loops) to potential functional roles in butyrate kinase (e.g., near the active site or substrate-binding region).

> **Figure 1.** Loop regions identified by `script1_loops.py` mapped onto butyrate kinase (1SAZ); loops are colored to highlight regions with higher average B-factors.

---

## Script 2 – Helix and Strand Distance Measurements (`script2_secondary_distances.py`)

### Objective

Script 2 analyzes the same 1SAZ structure to measure:

* The **distance between CA atoms** of the first and last residue of each alpha helix.
* The same distance for each beta strand.

This provides a compact geometric descriptor of secondary-structure elements in butyrate kinase.

### Step 1: Extract CA Coordinates

The script first reads all `ATOM` lines and keeps only the CA atoms:

```python
ca_coords = parse_ca_coordinates(pdb_path)
```

For each CA atom, it stores (x, y, z) coordinates keyed by `(chain_id, res_seq)`.

### Step 2: Parse Helix and Sheet Definitions

Helix and strand intervals are parsed from the PDB `HELIX` and `SHEET` records:

```python
helices, sheets = parse_helices_and_sheets(pdb_path)
```

Each item stores:

* An ID (helix index or sheet label).
* Chain ID.
* Start residue index.
* End residue index.

For example:

```text
Helix 1: Chain A, residues 40–55
Sheet A: Chain A, residues 100–105
```

### Step 3: Compute Start–End Distances

For each helix or strand, the script locates the CA coordinates of the start and end residues and computes the Euclidean distance:

```python
d = distance(ca_coords[start_key], ca_coords[end_key])
```


$$
d = \sqrt{(x_2 - x_1)^2 + (y_2 - y_1)^2 + (z_2 - z_1)^2}
$$

### Step 4: Output and Interpretation

Example output:

```text
Alpha helix distances:
  Helix 1 (Chain A 40–55): CA–CA distance = 22.45 Å
  Helix 2 (Chain A 60–70): CA–CA distance = 16.83 Å

Beta strand distances:
  Strand A (Chain A 100–105): CA–CA distance = 10.32 Å
  ...
```

These distances can be compared among helices and strands. Longer distances can imply extended helices/strands; shorter distances might indicate kinks or bends. Mapping these onto 1SAZ (Figure 2) provides structural context for where **rigid secondary-structure segments** occur relative to flexible loops identified in Script 1.

> **Figure 2.** Helices and beta strands in butyrate kinase (1SAZ) colored by start–end CA distance computed by `script2_secondary_distances.py`.

---

## Script 3 – Species Grouping and Amino-Acid Usage (`script3.py`)

### Objective

Script 3 integrates sequence and annotation data to examine how residues vary across species:

1. **Assigns a species-level group** (e.g., *Bacillus velezensis*) to each sequence using the `Organism` column in the TSV file.
2. Writes a new table with an added `Group` column (`grouped_table.tsv`).
3. Uses a CLUSTAL MSA to **calculate amino-acid usage at each alignment position for each group**, printing both **counts and percentages** for each residue.

### Command-Line Usage

```bash
python script3.py \
  --table ungrouped_table.tsv \
  --msa aligned.txt \
  --out_table grouped_table.tsv \
  > enzyme_analysis_grouped.txt
```

* `--table` – TSV with at least `Entry Name` (or `Entry`) and `Organism`.
* `--msa` – CLUSTAL alignment file for the butyrate kinase family.
* `--out_table` – output TSV that will include the new `Group` column.

### Step 1: Derive Species Group from Organism Names

The script reads each row of the TSV and uses a helper function to convert the full `Organism` string into a species-level group:

```python
species_group = extract_species(organism)
```

Key behaviors:

* Parenthetical information is removed:
  `"Bacillus velezensis (strain XYZ)" → "Bacillus velezensis"`.
* The first two tokens are used for the species name.
* If `Organism` is missing or empty, the group defaults to `"Unknown"`.

This species name is stored in a new `Group` column:

```python
row["Group"] = species_group
entry_to_group[entry_name] = species_group
```

The script then writes `grouped_table.tsv`, which mirrors the original table but with `Group` added.

### Step 2: Parse the CLUSTAL Alignment

The script reads the alignment:

```python
seq_dict = parse_clustal(args.msa)
```

* Skips the CLUSTAL header, blank lines, and consensus lines.
* Reconstructs full aligned sequences for each identifier.
* Uses a dictionary mapping `sequence_id → aligned_sequence`.

Example:

```text
tr|A0AAW7KBZ3|A0AAW7KBZ3_ENTFL → "------METVLVINPGSTSTKLA..."
```

### Step 3: Map MSA IDs to Groups

Sequence IDs in the CLUSTAL file are mapped to entry names and then to species groups:

```python
entry_name = clustal_id_to_entry_name(raw_id)
group = entry_to_group.get(entry_name, "Unknown")
group_to_seqs[group].append(aligned_seq)
```

* `clustal_id_to_entry_name` typically strips the prefix:

  * `tr|A0AAW7KBZ3|A0AAW7KBZ3_ENTFL` → `A0AAW7KBZ3_ENTFL`
* If no match is found in the table, the sequence is assigned to `Group = "Unknown"`.
  This is the source of the **"Group Unknown"** block in the output for sequences whose metadata was not present in the TSV or did not match the expected entry-name format.

The script prints a summary of how many sequences fall into each group (species).

### Step 4: Compute Amino-Acid Usage Per Position Per Group

The main analysis step is:

```python
compute_usage_per_position(group_to_seqs)
```

For each alignment position:

1. The script iterates over all groups (species).
2. For each group:

   * It looks at the residue at that position for **every sequence in that group**.
   * It counts residues using a `Counter`.
   * It computes `count / n * 100` to get percentages, where `n` is the number of sequences in that group.
   * Gaps (`-`) are labeled as `gap` and reported as well.

Example output (simplified):

```text
Position 181:
  Group Bacillus velezensis (n=2):
    K: 2/2 (100.0%)
  Group Unknown (n=19):
    K: 13/19 (68.4%)
    R: 6/19 (31.6%)
```

This confirms that the script:

* **Calculates amino acid usage at each alignment position for each group**, and
* **Prints both counts and percentages** for each residue.

Multiple positions show where particular residues are conserved or variable among species groups and among the `"Unknown"` group. Positions that are highly conserved in *Bacillus* but variable elsewhere may be functionally important, especially if they map close to loop regions or the active site in 1SAZ.

### Step 5: Interpretation for Butyrate Kinases

* Positions where **butyrate kinase sequences from the same species share the same residue** but other groups differ may point to **species-specific adaptations**.
* Positions that are conserved across almost all groups (e.g., 100% of sequences across species) likely indicate **core catalytic or structural residues**.
* Combining this with information from Script 1 (loop flexibility) and Script 2 (secondary-structure geometry) helps connect **sequence-level divergence** with possible structural and functional consequences in the butyrate kinase family.

> **Figure 3.** Example of position-specific amino-acid usage from `enzyme_analysis_grouped.txt` (e.g., Positions 181–190), highlighting residues where *Bacillus velezensis* differs from the broader “Unknown” group.

---

## Summary

Across the three scripts, the analysis pipeline for 1SAZ butyrate kinase:

* **Script 1** identifies flexible loop regions and quantifies their average B-factors.
* **Script 2** measures geometric properties of helices and strands using CA–CA distances.
* **Script 3** groups sequences by species and **calculates amino-acid usage at each MSA position for each group**, printing counts and percentages.

Together, these steps connect:

1. **Structure** (loop locations, helix/strand geometry),
2. **Dynamics** (via B-factors), and
3. **Evolutionary variation** (species-specific residue usage),

to provide a multi-layered view of the butyrate kinase enzyme family centered on the 1SAZ structure.
