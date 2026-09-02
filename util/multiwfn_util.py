### Utility functions to set up Multiwfn computations and extract results from output files.

import os
import pandas as pd
import subprocess
import re

from config import WORKING_DIR, MULTIWFN_EXEC, MULTIWFN_IP_IN


def deprot_IP(run_name):
    """
    Calculate the ionisation potential for the deprotonated species.

    Runs Multiwfn on the molden.input file in the run directory and writes
    the output to multiwfn_IP_Deprot.out.

    :param run_name: Name of the run (used to locate the run directory).
    :return: None
    """
    directory_path = os.path.join(WORKING_DIR, run_name)
    multiwfn_command = (
        f"{MULTIWFN_EXEC} < {MULTIWFN_IP_IN} > multiwfn_IP_Deprot.out"
    )
    subprocess.run(multiwfn_command, shell=True, cwd=directory_path)


def reprot_IP(run_name):
    """
    Calculate the ionisation potential for the protonated (neutral) species.

    Runs Multiwfn on the molden.input file in the run directory and writes
    the output to multiwfn_IP_Reprot.out.

    :param run_name: Name of the run (used to locate the run directory).
    :return: None
    """
    directory_path = os.path.join(WORKING_DIR, run_name)
    multiwfn_command = (
        f"{MULTIWFN_EXEC} < {MULTIWFN_IP_IN} > multiwfn_IP_Reprot.out"
    )
    subprocess.run(multiwfn_command, shell=True, cwd=directory_path)


def extract_IP_anion(run_name, charge):
    """
    Extract per-atom ionisation potential values for the deprotonated (anionic) species.

    Reads the Multiwfn output file and the corresponding xyz file, then extracts
    the minimal IP value for every nitrogen atom.

    :param run_name: Name of the run (used to locate the run directory).
    :param charge: Formal charge of the input molecule (0 or -1).
    :return: DataFrame with columns ``nitrogen_index`` and ``IP_anion``.
    """
    directory_path = os.path.join(WORKING_DIR, run_name)
    multiwfn_path = os.path.join(directory_path, "multiwfn_IP_Deprot.out")

    if charge == 0:
        xyz_path = os.path.join(directory_path, "deprotonated.xyz")
    elif charge == -1:
        xyz_path = os.path.join(directory_path, "xtbopt.xyz")

    # Number of atoms from xyz file
    with open(xyz_path, "r") as f:
        nat = int(f.readlines()[0].strip())

    # Extract ionisation potential table from Multiwfn output file
    try:
        datalines = []
        datalines3 = []
        with open(multiwfn_path, "r") as mf:
            for line in mf:
                if "Minimal, " in line:
                    for _ in range(nat + 1):
                        nextline = next(mf).strip()
                        if not nextline:
                            break
                        datalines.append(nextline)
                if "(The library is available at" in line:
                    for _ in range(nat + 1):
                        nextline = next(mf).strip()
                        if "The number of total inner-core electrons:" in nextline:
                            break
                        datalines3.append(nextline)

        df_multiwfn = pd.DataFrame([l.split() for l in datalines]).iloc[1:]
        df_multiwfn.columns = [
            "Atom#", "Area(Ang^2", "Min value", "Max value",
            "Average", "Variance", "Extra1", "Extra2",
        ]
        df_multiwfn = df_multiwfn.drop(["Extra1", "Extra2"], axis=1)

        df_element = pd.DataFrame([re.split(r"\(|\)", l) for l in datalines3])[[0, 1]]
        df_element.columns = ["Element", "Atom#"]
        df_element["Atom#"] = df_element["Atom#"].str.strip()
        df_element["Element"] = df_element["Element"].str.strip()

        df_result = pd.merge(df_multiwfn, df_element, on="Atom#", how="outer")
    except Exception:
        print("Error with Multiwfn output extraction or no file")

    try:
        data = [
            [row["Atom#"], row["Min value"]]
            for _, row in df_result.iterrows()
            if row["Element"] == "N"
        ]
    except Exception:
        print("Error with extracting nitrogen rows")

    IP_anion = pd.DataFrame(data, columns=["nitrogen_index", "IP_anion"])
    return IP_anion


def extract_IP_neutral(run_name):
    """
    Extract per-atom ionisation potential values for the protonated (neutral) species.

    Reads the Multiwfn output file and the protonated xyz file, then extracts
    the minimal IP value for every nitrogen atom.

    :param run_name: Name of the run (used to locate the run directory).
    :return: DataFrame with columns ``nitrogen_index`` and ``IP_neutral``.
    """
    directory_path = os.path.join(WORKING_DIR, run_name)
    multiwfn_path = os.path.join(directory_path, "multiwfn_IP_Reprot.out")
    xyz_path = os.path.join(directory_path, "protonated.xyz")

    # Number of atoms from xyz file
    with open(xyz_path, "r") as f:
        nat = int(f.readlines()[0].strip())

    # Extract ionisation potential table from Multiwfn output file
    try:
        datalines = []
        datalines3 = []
        with open(multiwfn_path, "r") as mf:
            for line in mf:
                if "Minimal, " in line:
                    for _ in range(nat + 1):
                        nextline = next(mf).strip()
                        if not nextline:
                            break
                        datalines.append(nextline)
                if "(The library is available at" in line:
                    for _ in range(nat + 1):
                        nextline = next(mf).strip()
                        if "The number of total inner-core electrons:" in nextline:
                            break
                        datalines3.append(nextline)

        df_multiwfn = pd.DataFrame([l.split() for l in datalines]).iloc[1:]
        df_multiwfn.columns = [
            "Atom#", "Area(Ang^2", "Min value", "Max value",
            "Average", "Variance", "Extra1", "Extra2",
        ]
        df_multiwfn = df_multiwfn.drop(["Extra1", "Extra2"], axis=1)

        df_element = pd.DataFrame([re.split(r"\(|\)", l) for l in datalines3])[[0, 1]]
        df_element.columns = ["Element", "Atom#"]
        df_element["Atom#"] = df_element["Atom#"].str.strip()
        df_element["Element"] = df_element["Element"].str.strip()

        df_result = pd.merge(df_multiwfn, df_element, on="Atom#", how="outer")
    except Exception:
        print("Error with Multiwfn output extraction or no file")

    try:
        data = [
            [row["Atom#"], row["Min value"]]
            for _, row in df_result.iterrows()
            if row["Element"] == "N"
        ]
    except Exception:
        print("Error with extracting nitrogen rows")

    IP_neutral = pd.DataFrame(data, columns=["nitrogen_index", "IP_neutral"])
    return IP_neutral
