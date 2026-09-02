"""Cheminformatics utility functions for input validation and result visualisation."""

import os
import re

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.Draw import rdMolDraw2D

from config import WORKING_DIR


def check_reactiveN(substrate):
    """
    Check that the SMILES contains more than one reactive nitrogen.

    A nitrogen is considered reactive if it is bonded to fewer than three
    non-hydrogen atoms and has no triple bond to any neighbour.

    :param substrate: Canonicalised SMILES string.
    :return: True if the structure contains more than one reactive nitrogen,
             False otherwise.
    :rtype: bool
    """
    reactant = Chem.MolFromSmiles(substrate)
    n_reactive_nitrogen = 0

    for atom in reactant.GetAtoms():
        if atom.GetSymbol() == "N":
            neighbors = atom.GetNeighbors()
            n_heavy_atom = 0
            triple_bond = False
            for neighbor in neighbors:
                if neighbor.GetSymbol() != "H":
                    n_heavy_atom += 1
                bond = reactant.GetBondBetweenAtoms(neighbor.GetIdx(), atom.GetIdx())
                if bond.GetBondType() == Chem.BondType.TRIPLE:
                    triple_bond = True
            if n_heavy_atom < 3 and not triple_bond:
                n_reactive_nitrogen += 1

    return n_reactive_nitrogen > 1


def check_MW(substrate):
    """
    Check that the molecular weight of the substrate is below 500 g/mol.

    :param substrate: Canonicalised SMILES string.
    :return: True if MW < 500, False otherwise.
    :rtype: bool
    """
    mol = Chem.MolFromSmiles(substrate)
    return Descriptors.MolWt(mol) < 500


def check_no_salt(substrate):
    """
    Check that the SMILES does not represent a salt (i.e. contains no dot separator).

    :param substrate: Canonicalised SMILES string.
    :return: True if there is no salt, False otherwise.
    :rtype: bool
    """
    return "." not in substrate


def check_no_metal(substrate):
    """
    Check that the SMILES does not contain a metal atom.

    :param substrate: Canonicalised SMILES string.
    :return: True if no metal is present, False otherwise.
    :rtype: bool
    """
    metal_symbols = [
        "Li", "Na", "K", "Rb", "Cs", "Fr", "Be", "Mg", "Ca", "Sr", "Ba", "Ra",
        "Sc", "Y", "Ti", "Zr", "Hf", "V", "Nb", "Ta", "Cr", "Mo", "W", "Mn", "Tc",
        "Re", "Fe", "Ru", "Os", "Co", "Rh", "Ir", "Ni", "Pd", "Pt", "Cu", "Ag",
        "Au", "Zn", "Cd", "Hg", "Al", "Ga", "In", "Tl", "Sn", "Pb", "Bi",
    ]
    metal_pattern = "|".join(metal_symbols)
    return not re.compile(metal_pattern).search(substrate)


def check_charge(substrate):
    """
    Calculate the formal charge of the substrate molecule.

    :param substrate: Canonicalised SMILES string.
    :return: Formal charge of the molecule.
    :rtype: int
    """
    mol = Chem.MolFromSmiles(substrate)
    return Chem.GetFormalCharge(mol)


def draw_atom_label(mol, df, proba, run_name):
    """
    Draw the molecule with per-nitrogen N-alkylation probabilities annotated.

    Saves the image as ``N_alkylation_probability.png`` in the run directory.
    Note: currently only works correctly when there is no change in atom index
    during the deprotonation step.

    :param mol: RDKit molecule object (constructed from the original SMILES).
    :param df: Feature matrix DataFrame containing the ``nitrogen_index`` column
               (indices taken from the deprotonated structure).
    :param proba: Array of class probabilities from the model (shape: n_nitrogens x 2).
    :param run_name: Name of the run (used to locate the output directory).
    :return: None
    """
    directory_path = os.path.join(WORKING_DIR, run_name)

    for idx, row in df.iterrows():
        nitrogen_index = row["nitrogen_index"] - 1
        for atom in mol.GetAtoms():
            if atom.GetIdx() == nitrogen_index:
                atom.SetProp(
                    "atomNote",
                    f"{atom.GetIdx() + 1}: {proba[idx][1]:.2f}",
                )

    d = rdMolDraw2D.MolDraw2DCairo(800, 600)
    d.drawOptions().addAtomIndicies = False
    d.SetFontSize(400)
    d.DrawMolecule(mol)
    d.FinishDrawing()

    output_path = os.path.join(directory_path, "N_alkylation_probability.png")
    with open(output_path, "wb") as f:
        f.write(d.GetDrawingText())


def check_atom_index_change(run_name):
    """
    Check whether any nitrogen atom index changed during the reprotonation step.

    Parses the ``Reprot.out`` file produced by CREST and looks for the table
    of old-to-new atom index mappings. Prints a warning if any nitrogen atom
    index has changed.

    :param run_name: Name of the run (used to locate the run directory).
    :return: None
    """
    directory_path = os.path.join(WORKING_DIR, run_name)
    reprot_out_file = os.path.join(directory_path, "Reprot.out")
    found_table = False

    with open(reprot_out_file, "r") as mf:
        for line in mf:
            if "Input coordinate lines sorted:" in line:
                found_table = True
                table = []
                while True:
                    next_line = next(mf)
                    if next_line.strip():
                        table.append(next_line)
                    else:
                        break

                table_data = [row.split() for row in table]
                for row in table_data:
                    if len(row) > 2 and row[0] == "N" and row[1] != row[2]:
                        print("Warning: Nitrogen index has changed after reprotonation")
