"""Utility functions for running xTB and CREST calculations."""

import os
import subprocess

from config import WORKING_DIR, XTB_EXEC, CREST_EXEC


def deprot_opt(run_name, charge):
    """
    Deprotonate neutral molecules or optimise the structure of anionic molecules.

    For neutral molecules (charge == 0) CREST is used to generate a deprotonated
    structure (``deprotonated.xyz``). For anionic molecules (charge == -1) xTB
    performs a geometry optimisation (``xtbopt.xyz``). Output is written to
    ``Deprot.out``.

    :param run_name: Name of the run (used to locate the run directory).
    :param charge: Formal charge of the input molecule (0 or -1).
    :return: None
    """
    directory_path = os.path.join(WORKING_DIR, run_name)
    xyz_path = os.path.join(directory_path, f"{run_name}.xyz")

    if charge == 0:
        crest_command = (
            f"{CREST_EXEC} {xyz_path} --deprotonate --alpb h2o --chrg 0 > Deprot.out"
        )
        process = subprocess.run(crest_command, shell=True, cwd=directory_path)
        assert process.returncode == 0, "CREST deprotonation failed"
    elif charge == -1:
        xtb_command = (
            f"{XTB_EXEC} {xyz_path} --opt --alpb h2o --chrg -1 > Deprot.out"
        )
        process = subprocess.run(xtb_command, shell=True, cwd=directory_path)
        assert process.returncode == 0, "xTB optimisation failed"


def deprot_molden(run_name, charge):
    """
    Generate a molden.input file for the optimised anionic structure.

    xTB is run with the ``--molden`` flag on the deprotonated or optimised xyz
    file. Output is written to ``xtb_Deprot.out``.

    :param run_name: Name of the run (used to locate the run directory).
    :param charge: Formal charge of the input molecule (0 or -1).
    :return: None
    """
    directory_path = os.path.join(WORKING_DIR, run_name)

    if charge == 0:
        xtb_command = (
            f"{XTB_EXEC} deprotonated.xyz --molden --alpb h2o --chrg -1 > xtb_Deprot.out"
        )
    elif charge == -1:
        xtb_command = (
            f"{XTB_EXEC} xtbopt.xyz --molden --alpb h2o --chrg -1 > xtb_Deprot.out"
        )

    subprocess.run(xtb_command, shell=True, cwd=directory_path)


def reprot_opt(run_name, charge):
    """
    Protonate all molecules to form a neutral species using CREST.

    The deprotonated or optimised xyz file is used as input. Output is written
    to ``Reprot.out``.

    :param run_name: Name of the run (used to locate the run directory).
    :param charge: Formal charge of the input molecule (0 or -1).
    :return: None
    """
    directory_path = os.path.join(WORKING_DIR, run_name)

    if charge == 0:
        xyz_path = os.path.join(directory_path, "deprotonated.xyz")
    elif charge == -1:
        xyz_path = os.path.join(directory_path, "xtbopt.xyz")

    crest_command = (
        f"{CREST_EXEC} {xyz_path} --protonate --alpb h2o --chrg -1 > Reprot.out"
    )
    subprocess.run(crest_command, shell=True, cwd=directory_path)


def reprot_molden(run_name):
    """
    Generate a molden.input file for the optimised neutral structure.

    xTB is run with the ``--molden`` flag on ``protonated.xyz``. Output is
    written to ``xtb_Reprot.out``.

    :param run_name: Name of the run (used to locate the run directory).
    :return: None
    """
    directory_path = os.path.join(WORKING_DIR, run_name)
    xtb_command = (
        f"{XTB_EXEC} protonated.xyz --molden --alpb h2o --chrg 0 > xtb_Reprot.out"
    )
    subprocess.run(xtb_command, shell=True, cwd=directory_path)
