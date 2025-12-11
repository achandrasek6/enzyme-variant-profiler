#!/usr/bin/env python3
"""
script1_loops.py

Identify loop regions in a PDB file and report:
- average B-factor per loop
- residue names and sequence numbers in each loop

Usage:
    python script1_loops.py 1SAZ.pdb
"""

import sys
from collections import defaultdict

def parse_pdb_atoms(pdb_path):
    """
    Return:
      residue_atoms: dict[(chain_id, res_seq)] -> {
          "resname": str,
          "B_factors": [float, ...]
      }
    """
    residue_atoms = {}
    with open(pdb_path) as f:
        for line in f:
            record = line[0:6].strip()
            if record != "ATOM":
                continue
            resname = line[17:20].strip()
            chain_id = line[21].strip() or "-"
            res_seq = int(line[22:26])
            try:
                b_factor = float(line[60:66])
            except ValueError:
                continue

            key = (chain_id, res_seq)
            if key not in residue_atoms:
                residue_atoms[key] = {"resname": resname, "B_factors": []}
            residue_atoms[key]["B_factors"].append(b_factor)
    return residue_atoms

def parse_secondary_structure_regions(pdb_path):
    """
    Parse HELIX and SHEET records and return a set of residues
    that are part of any secondary structure element.
    """
    sec_residues = set()
    with open(pdb_path) as f:
        for line in f:
            record = line[0:6].strip()
            if record == "HELIX":
                init_chain = line[19].strip() or "-"
                init_seq = int(line[21:25])
                end_chain = line[31].strip() or init_chain
                end_seq = int(line[33:37])
                if init_chain != end_chain:
                    chain_ids = [init_chain, end_chain]
                else:
                    chain_ids = [init_chain]
                for chain in chain_ids:
                    for res in range(init_seq, end_seq + 1):
                        sec_residues.add((chain, res))

            elif record == "SHEET":
                init_chain = line[21].strip() or "-"
                init_seq = int(line[22:26])
                end_chain = line[32].strip() or init_chain
                end_seq = int(line[33:37])
                if init_chain != end_chain:
                    chain_ids = [init_chain, end_chain]
                else:
                    chain_ids = [init_chain]
                for chain in chain_ids:
                    for res in range(init_seq, end_seq + 1):
                        sec_residues.add((chain, res))
    return sec_residues

def find_loops(residue_atoms, sec_residues):
    """
    Identify continuous stretches of residues that are NOT in sec_residues.
    Returns a list of loops, where each loop is a list of keys (chain, res_seq).
    """
    residues_by_chain = defaultdict(list)
    for (chain, res_seq) in residue_atoms.keys():
        residues_by_chain[chain].append(res_seq)

    loops = []
    for chain, res_list in residues_by_chain.items():
        res_list = sorted(set(res_list))
        current_loop = []
        prev_res = None

        for res in res_list:
            in_secondary = (chain, res) in sec_residues
            if in_secondary:
                if current_loop:
                    loops.append(current_loop)
                    current_loop = []
                prev_res = res
                continue

            # Residue is NOT in helix/sheet
            if prev_res is None or res != prev_res + 1 or (chain, prev_res) in sec_residues:
                if current_loop:
                    loops.append(current_loop)
                current_loop = [(chain, res)]
            else:
                current_loop.append((chain, res))

            prev_res = res

        if current_loop:
            loops.append(current_loop)

    return loops

def main():
    if len(sys.argv) != 2:
        print("Usage: python script1_loops.py <pdb_file>")
        sys.exit(1)

    pdb_path = sys.argv[1]
    residue_atoms = parse_pdb_atoms(pdb_path)
    sec_residues = parse_secondary_structure_regions(pdb_path)
    loops = find_loops(residue_atoms, sec_residues)

    if not loops:
        print("No loop regions found (all residues are in helices/sheets or no atoms parsed).")
        return

    for i, loop in enumerate(loops, start=1):
        b_values = []
        residue_labels = []
        for (chain, res_seq) in loop:
            info = residue_atoms.get((chain, res_seq))
            if not info:
                continue
            b_values.extend(info["B_factors"])
            residue_labels.append(f"{info['resname']}{res_seq}{chain}")

        if not b_values:
            avg_b = float("nan")
        else:
            avg_b = sum(b_values) / len(b_values)

        start_chain, start_res = loop[0]
        end_chain, end_res = loop[-1]
        print(f"Loop {i}: Chain {start_chain}, residues {start_res}–{end_res}")
        print(f"  Average B-factor: {avg_b:.2f} Å^2")
        print("  Residues:")
        print("   ", ", ".join(residue_labels))
        print()

if __name__ == "__main__":
    main()
