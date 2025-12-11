#!/usr/bin/env python3
"""
Script 3: Enzyme variant analysis by group (species)

1. Takes a table of selected entries (TSV) and an MSA (CLUSTAL) as input.
2. Creates a 'Group' column for each entry (each group = species from the Organism column).
3. Calculates amino acid usage at each alignment position for each group.
4. Prints counts and percentages like:
   Amino acid usage per position per group:
   Position 1:
     Group Bacillus subtilis (n=7):
       M: 7/7 (100.0%)
       gap: 0/7 (0.0%)
   ...

Usage example:
    python script3.py --table ungrouped_table.tsv --msa aligned.txt --out_table grouped_table.tsv > enzyme_analysis_grouped.txt
"""

import argparse
import csv
import re
from collections import defaultdict, Counter


def parse_args():
    p = argparse.ArgumentParser(description="Script 3: AA usage per position per group (species).")
    p.add_argument(
        "--table",
        required=True,
        help="TSV table with at least 'Entry Name' and 'Organism' columns."
    )
    p.add_argument(
        "--msa",
        required=True,
        help="Multiple sequence alignment in CLUSTAL format."
    )
    p.add_argument(
        "--out_table",
        default="grouped_table.tsv",
        help="Output TSV with an added 'Group' column (default: grouped_table.tsv)."
    )
    return p.parse_args()


def extract_species(organism_str: str) -> str:
    """
    Take an Organism string like:
        'Clostridioides difficile (strain R20291)'
    and return a species-level group:
        'Clostridioides difficile'
    """
    if not organism_str:
        return "Unknown"

    # Remove any parentheses (strain info, etc.)
    cleaned = re.sub(r"\([^)]*\)", "", organism_str).strip()
    tokens = cleaned.split()
    if len(tokens) >= 2:
        return " ".join(tokens[:2])
    elif tokens:
        return tokens[0]
    else:
        return "Unknown"


def parse_clustal(msa_path: str):
    """
    Parse a CLUSTAL alignment and return:
        seq_dict: { raw_id -> full_aligned_sequence }
    """
    seq_dict = {}
    with open(msa_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            # Skip header, consensus, and blank lines
            if not line or line.startswith("CLUSTAL") or line[0].isspace():
                continue

            parts = line.split()
            if len(parts) < 2:
                continue

            seq_id = parts[0]
            frag = parts[1]

            seq_dict.setdefault(seq_id, "")
            seq_dict[seq_id] += frag

    return seq_dict


def clustal_id_to_entry_name(seq_id: str) -> str:
    """
    Convert a CLUSTAL sequence ID to UniProt entry name.

    Typical CLUSTAL IDs from UniProt look like:
        tr|A0A9D2L4V7|A0A9D2L4V7_9BACT
    We want:
        A0A9D2L4V7_9BACT
    """
    parts = seq_id.split("|")
    if len(parts) >= 3:
        return parts[2]
    return parts[-1]


def find_column_case_insensitive(fieldnames, target):
    """
    Return the actual column name in fieldnames that matches target (case-insensitive),
    or None if not found.
    """
    target_lower = target.lower()
    for name in fieldnames:
        if name.lower() == target_lower:
            return name
    return None


def load_table_and_add_groups(table_path: str):
    """
    Read the UniProt table (TSV) and:
      - add a 'Group' column = species extracted from 'Organism'
      - build mapping entry_name -> group

    Returns:
      rows: list[dict]
      entry_to_group: dict[entry_name] -> group (species)
      fieldnames: original header list (we'll add 'Group' if needed)
    """
    rows = []
    entry_to_group = {}

    with open(table_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fieldnames = reader.fieldnames
        if fieldnames is None:
            raise ValueError("Input table has no header row.")

        # Try to locate 'Entry Name' and 'Organism' flexibly
        entry_col = find_column_case_insensitive(fieldnames, "Entry Name") or \
                    find_column_case_insensitive(fieldnames, "Entry")
        org_col = find_column_case_insensitive(fieldnames, "Organism")

        if entry_col is None or org_col is None:
            raise ValueError(
                "Table must contain columns for Entry Name and Organism "
                "(case-insensitive)."
            )

        for row in reader:
            organism = row.get(org_col, "")
            species_group = extract_species(organism)
            row["Group"] = species_group
            rows.append(row)

            entry_name = row.get(entry_col)
            if entry_name:
                entry_to_group[entry_name] = species_group

    return rows, entry_to_group, fieldnames


def write_grouped_table(rows, fieldnames, out_path):
    """Write rows back out as TSV with 'Group' column included."""
    if "Group" not in fieldnames:
        fieldnames = list(fieldnames) + ["Group"]

    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def compute_usage_per_position(group_to_seqs):
    """
    Given:
        group_to_seqs: { group_name -> [aligned_seq1, aligned_seq2, ...] }

    Compute AA usage per position for each group and print counts + percentages.
    """
    # Make sure all sequences are the same length
    lengths = {len(seq) for seqs in group_to_seqs.values() for seq in seqs}
    if len(lengths) != 1:
        raise ValueError(f"Sequences are not all the same length: lengths = {lengths}")
    aln_length = lengths.pop()

    print("Amino acid usage per position per group:\n")

    # For each alignment position (1-based for printing)
    for pos_idx in range(aln_length):
        pos_num = pos_idx + 1
        print(f"Position {pos_num}:")

        # For each group, gather the residues at this position
        for group in sorted(group_to_seqs.keys()):
            seqs = group_to_seqs[group]
            n = len(seqs)
            residues = [seq[pos_idx] for seq in seqs]

            counts = Counter(residues)

            print(f"  Group {group} (n={n}):")

            # Sort keys so gaps come last, others alphabetical
            def sort_key(res):
                return (res == "-", res)

            for res in sorted(counts.keys(), key=sort_key):
                count = counts[res]
                aa_label = "gap" if res == "-" else res
                pct = (count / n) * 100.0
                print(f"    {aa_label}: {count}/{n} ({pct:.1f}%)")

        print()  # blank line between positions


def main():
    args = parse_args()

    # 1. Load table, add Group (species) column, and build entry_name -> group mapping
    rows, entry_to_group, fieldnames = load_table_and_add_groups(args.table)

    # 2. Parse the CLUSTAL MSA
    seq_dict = parse_clustal(args.msa)

    # 3. Map each MSA sequence ID to Entry Name, then to Group
    group_to_seqs = defaultdict(list)
    unmapped = []

    for raw_id, aligned_seq in seq_dict.items():
        entry_name = clustal_id_to_entry_name(raw_id)
        group = entry_to_group.get(entry_name, "Unknown")
        if group == "Unknown":
            unmapped.append(entry_name)
        group_to_seqs[group].append(aligned_seq)

    # 4. Write grouped table (for the assignment requirement: add Group column)
    write_grouped_table(rows, fieldnames, args.out_table)
    print(f"Grouped table written to: {args.out_table}\n")

    print("Groups (species) represented in the MSA:")
    for group, seqs in sorted(group_to_seqs.items(), key=lambda x: x[0]):
        print(f"  {group}: {len(seqs)} sequences")
    print()

    if unmapped:
        print("WARNING: Some MSA entries were not found in the table and were")
        print("assigned to group 'Unknown':")
        for eid in sorted(set(unmapped)):
            print(f"  {eid}")
        print()

    # 5. Calculate and print amino acid usage per position per group
    compute_usage_per_position(group_to_seqs)


if __name__ == "__main__":
    main()
