#!/usr/bin/env python3
"""
script2_secondary_distances.py

For each helix and beta sheet in a PDB file, compute the distance between
the alpha-carbon (CA) atoms of the initial and terminal residues.

Usage:
    python script2_secondary_distances.py 1SAZ.pdb
"""

import sys
import math

def parse_ca_coordinates(pdb_path):
    """
    Return dict[(chain_id, res_seq)] -> (x, y, z) for CA atoms.
    """
    ca_coords = {}
    with open(pdb_path) as f:
        for line in f:
            record = line[0:6].strip()
            if record != "ATOM":
                continue
            atom_name = line[12:16].strip()
            if atom_name != "CA":
                continue
            chain_id = line[21].strip() or "-"
            res_seq = int(line[22:26])
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            ca_coords[(chain_id, res_seq)] = (x, y, z)
    return ca_coords

def parse_helices_and_sheets(pdb_path):
    """
    Parse HELIX and SHEET records.

    Returns:
      helices: list of (helix_id, chain_id, start_res_seq, end_res_seq)
      sheets: list of (sheet_id, chain_id, start_res_seq, end_res_seq)
    """
    helices = []
    sheets = []
    with open(pdb_path) as f:
        for line in f:
            record = line[0:6].strip()
            if record == "HELIX":
                helix_id = line[7:10].strip()
                chain_id = line[19].strip() or "-"
                start_res = int(line[21:25])
                end_chain = line[31].strip() or chain_id
                end_res = int(line[33:37])
                helices.append((helix_id or str(len(helices)+1), chain_id, start_res, end_res))

            elif record == "SHEET":
                sheet_id = line[11:14].strip()
                chain_id = line[21].strip() or "-"
                start_res = int(line[22:26])
                end_chain = line[32].strip() or chain_id
                end_res = int(line[33:37])
                sheets.append((sheet_id or str(len(sheets)+1), chain_id, start_res, end_res))
    return helices, sheets

def distance(p1, p2):
    return math.sqrt(
        (p2[0] - p1[0]) ** 2 +
        (p2[1] - p1[1]) ** 2 +
        (p2[2] - p1[2]) ** 2
    )

def main():
    if len(sys.argv) != 2:
        print("Usage: python script2_secondary_distances.py <pdb_file>")
        sys.exit(1)

    pdb_path = sys.argv[1]
    ca_coords = parse_ca_coordinates(pdb_path)
    helices, sheets = parse_helices_and_sheets(pdb_path)

    if not helices:
        print("No HELIX records found in this PDB file.")
    if not sheets:
        print("No SHEET records found in this PDB file.")

    if not helices and not sheets:
        return

    print("Alpha helix distances:")
    for helix_id, chain, start_res, end_res in helices:
        start_key = (chain, start_res)
        end_key = (chain, end_res)
        if start_key not in ca_coords or end_key not in ca_coords:
            print(f"  Helix {helix_id} (Chain {chain} {start_res}-{end_res}): CA coordinates missing.")
            continue
        d = distance(ca_coords[start_key], ca_coords[end_key])
        print(f"  Helix {helix_id} (Chain {chain} {start_res}-{end_res}): "
              f"distance between CA atoms = {d:.2f} Å")

    print("\nBeta sheet strand distances:")
    for sheet_id, chain, start_res, end_res in sheets:
        start_key = (chain, start_res)
        end_key = (chain, end_res)
        if start_key not in ca_coords or end_key not in ca_coords:
            print(f"  Sheet {sheet_id} (Chain {chain} {start_res}-{end_res}): CA coordinates missing.")
            continue
        d = distance(ca_coords[start_key], ca_coords[end_key])
        print(f"  Sheet {sheet_id} (Chain {chain} {start_res}-{end_res}): "
              f"distance between CA atoms = {d:.2f} Å")

if __name__ == "__main__":
    main()
