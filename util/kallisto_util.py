"""Utility functions for running kallisto descriptor calculations."""

import os

from kallisto.molecule import Molecule
import kallisto.reader.strucreader as ksr

from config import WORKING_DIR


def get_kallisto_neutral(run_name):
    """
    Calculate Proximity Shells and Partial Charge (EEQ) for the neutral species.

    Uses the ``protonated.xyz`` file in the run directory.

    :param run_name: Name of the run (used to locate the run directory).
    :return: Tuple of (proximity_shells, partial_charges) for every atom.
    """
    directory_path = os.path.join(WORKING_DIR, run_name)
    xyz_path = os.path.join(directory_path, "protonated.xyz")

    with open(xyz_path, "r") as f:
        atoms = ksr.read(f)
        molecule = Molecule(symbols=atoms)

        # Proximity shells (size 2 and 3)
        prox_all = molecule.get_prox((2, 3))

        # EEQ partial charges for the neutral species
        eeq_all = molecule.get_eeq(charge=0)

    return prox_all, eeq_all


def get_kallisto_anion(run_name, charge):
    """
    Calculate Polarizability for the anionic species.

    Uses ``deprotonated.xyz`` (charge 0 input) or ``xtbopt.xyz`` (charge -1 input).

    :param run_name: Name of the run (used to locate the run directory).
    :param charge: Formal charge of the input molecule (0 or -1).
    :return: Polarizability for every atom in the anionic species.
    """
    directory_path = os.path.join(WORKING_DIR, run_name)

    if charge == 0:
        xyz_path = os.path.join(directory_path, "deprotonated.xyz")
    elif charge == -1:
        xyz_path = os.path.join(directory_path, "xtbopt.xyz")

    with open(xyz_path, "r") as f:
        atoms = ksr.read(f)
        molecule = Molecule(symbols=atoms)
        alp_all = molecule.get_alp(charge=-1)

    return alp_all
