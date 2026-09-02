"""Utility functions for running OpenBabel file format conversions."""

import os
import subprocess

from config import WORKING_DIR, OBABEL_EXEC


def get_obabel_coords(run_name):
    """
    Convert a Corina mol2 file to an xyz coordinate file using OpenBabel.

    Reads ``<run_name>.mol2`` and writes ``<run_name>.xyz`` in the run directory.

    :param run_name: Name of the run (used to locate the run directory).
    :raises AssertionError: If the xyz output file was not created.
    :return: None
    """
    directory_path = os.path.join(WORKING_DIR, run_name)
    mol2_path = os.path.join(directory_path, f"{run_name}.mol2")
    xyz_path = os.path.join(directory_path, f"{run_name}.xyz")

    obabel_command = f"{OBABEL_EXEC} -imol2 {mol2_path} -oxyz > {xyz_path}"
    subprocess.run(obabel_command, shell=True)

    assert os.path.exists(xyz_path), "Warning: xyz file was not generated."
