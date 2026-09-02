"""Utility functions for running Corina 3D coordinate generation."""

import os
import subprocess

from config import WORKING_DIR, CORINA_EXEC


def corina(run_name):
    """
    Generate 3D coordinates from a SMILES file using Corina.

    Reads ``<run_name>.smi`` from the run directory and writes a mol2 file
    (``<run_name>.mol2``) to the same directory. The ``-d wh`` flag instructs
    Corina to add explicit hydrogen atoms.

    :param run_name: Name of the run (used to locate the run directory).
    :raises AssertionError: If the mol2 output file was not created or is empty.
    :return: None
    """
    directory_path = os.path.join(WORKING_DIR, run_name)
    in_path = os.path.join(directory_path, f"{run_name}.smi")
    out_mol2 = os.path.join(directory_path, f"{run_name}.mol2")

    if os.path.exists(in_path):
        corina_command = (
            f"{CORINA_EXEC} -i t=smiles {in_path} -o t=mol2 {out_mol2} -d wh"
        )
        subprocess.run(corina_command, shell=True)
    else:
        print("No SMILES file found")

    assert os.path.exists(out_mol2) and os.path.getsize(out_mol2) > 0, (
        "Warning: mol2 file was not generated or is empty."
    )
