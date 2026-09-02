"""Utility functions for generating molecular fingerprints."""

import numpy as np
import pandas as pd
from rdkit.Chem import AllChem as Chem


def generate_fingerprints(df, smiles_column):
    """
    Generate RDKit fingerprints (length 2048) for every row in a DataFrame.

    Each SMILES string in ``smiles_column`` is converted to an RDKit
    fingerprint and appended as bit-columns (``bit_0``, ``bit_1``, ...,
    ``bit_2047``) to the input DataFrame.

    :param df: DataFrame containing at least one SMILES column.
    :param smiles_column: Name of the column holding SMILES strings.
    :return: Copy of ``df`` with 2048 fingerprint bit-columns appended.
    :rtype: pandas.DataFrame
    """
    fingerprints = []

    for smiles in df[smiles_column]:
        mol = Chem.MolFromSmiles(smiles)
        fp = Chem.RDKFingerprint(mol)
        fingerprints.append([int(x) for x in fp.ToBitString()])

    fingerprint_df = pd.DataFrame(
        fingerprints, columns=[f"bit_{i}" for i in range(len(fingerprints[0]))]
    )
    return pd.concat([df, fingerprint_df], axis=1)
