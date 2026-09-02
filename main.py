# This is the main file for the N-alkylation/arylation predictor.
# It takes a nitrogen nucleophile molecule SMILES as input and
# gives the probabilities for each possible nitrogen reaction centre as output.

# TODO:
# The nitrogen indices are currently taken from the deprotonated structure. Check if the
# nitrogen index changes after reprotonation, so that it matches the indices from deprotonation.
# There may be a failure since there is a step that merges the dataframes from IP_anion and
# IP_neutral on 'nitrogen_index'. Suggestion: reorder the atoms in the xyz file at the
# beginning of the workflow.
# Stop the workflow if the check returns False (e.g. only 1 reactive nitrogen).
# Add checks if a step does not work and either reset or stop the workflow.


# 1 Imports
import os
import subprocess
import re

import pandas as pd
import numpy as np
import joblib

from rdkit import Chem
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import rdChemReactions as Reactions
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")
from rdkit.Chem import Descriptors

from kallisto.molecule import Molecule
import kallisto.reader.strucreader as ksr
import xgboost

from config import WORKING_DIR, MODEL_PATH
from util.cheminformatics import (
    check_reactiveN,
    check_MW,
    check_no_salt,
    check_no_metal,
    check_charge,
    draw_atom_label,
    check_atom_index_change,
)
from util.fp_util import generate_fingerprints
from util.corina_util import corina
from util.obabel_util import get_obabel_coords
from util.xtb_util import deprot_opt, deprot_molden, reprot_opt, reprot_molden
from util.multiwfn_util import deprot_IP, reprot_IP, extract_IP_anion, extract_IP_neutral
from util.kallisto_util import get_kallisto_anion, get_kallisto_neutral


# 2 Input from user
substrate = input("Molecule SMILES: ")
run_name = input("Name the molecule: ")

# 3 Canonicalise
substrate = Chem.MolToSmiles(Chem.MolFromSmiles(substrate))

# 4 Input checks
if not check_reactiveN(substrate):
    raise ValueError(
        "Warning: The substrate does not have more than one reactive nitrogen."
    )

if not check_MW(substrate):
    raise ValueError(
        "Warning: The MW of the substrate exceeds the limit for which this model was "
        "trained (MW > 500 g/mol)."
    )

if not check_no_salt(substrate):
    raise ValueError("Warning: Please submit a SMILES without a salt.")

if not check_no_metal(substrate):
    raise ValueError("Warning: The SMILES contains a metal.")

# 5 Calculate charge
charge = check_charge(substrate)
if charge not in [0, -1]:
    raise ValueError(
        f"Warning: The substrate charge is {charge}, but only 0 or -1 are supported."
    )

print("All checks complete")

# 6 Create run folder and write SMILES file for Corina
directory_path = os.path.join(WORKING_DIR, run_name)
os.makedirs(directory_path, exist_ok=True)
file_path = os.path.join(directory_path, f"{run_name}.smi")
with open(file_path, "w") as f:
    f.write(substrate)
print(f"Run folder created: {directory_path}")

# 7 Corina mol2 file generation
corina(run_name)
print("Corina generated mol2 file")

# 8 OpenBabel mol2 to xyz
get_obabel_coords(run_name)
print("OpenBabel converted mol2 file to xyz file")

# 9 CREST deprotonation and/or xTB optimisation
deprot_opt(run_name, charge)
print("CREST deprotonation complete")

# 10 Molden file generation for anion
deprot_molden(run_name, charge)

# 11 Ionisation potential for anion
deprot_IP(run_name)

# 12 Extract IP_anion
IP_anion = extract_IP_anion(run_name, charge)

# 13 CREST reprotonation
reprot_opt(run_name, charge)
print("CREST reprotonation complete")

# 14 Molden file generation for neutral (protonated)
reprot_molden(run_name)

# 15 Ionisation potential for neutral
reprot_IP(run_name)

# 16 Extract IP_neutral
IP_neutral = extract_IP_neutral(run_name)
print("Ionisation potential calculations complete")

# 17 Calculate kallisto descriptors:
#    - partial charge (neutral), proximity shells (neutral), polarizability (anion)
prox_all, eeq_all = get_kallisto_neutral(run_name)
alp_all = get_kallisto_anion(run_name, charge)

data = []
for atom_index, (prox, eeq) in enumerate(zip(prox_all, eeq_all)):
    data.append(
        {
            "nitrogen_index": atom_index + 1,
            "proximity_shells_neutral": prox,
            "partial_charge_neutral": eeq,
        }
    )
kallisto_neutral = pd.DataFrame(data).drop_duplicates(subset=["nitrogen_index"])

data = []
for atom_index, alp in enumerate(alp_all):
    data.append({"nitrogen_index": atom_index + 1, "polarizability_anion": alp})
kallisto_anion = pd.DataFrame(data).drop_duplicates(subset=["nitrogen_index"])

kallisto = kallisto_neutral.merge(kallisto_anion, on=["nitrogen_index"])
print("Kallisto calculations complete")

# 18 Check nitrogen indices after deprot/opt and reprot
check_atom_index_change(run_name)

# 19 Build feature matrix X
X_index = IP_anion.merge(IP_neutral[["nitrogen_index", "IP_neutral"]], on="nitrogen_index")
X_index["nitrogen_index"] = X_index["nitrogen_index"].astype(int)
X_index = X_index.merge(kallisto, on="nitrogen_index", how="left")
X_index["smiles"] = substrate
X_index = generate_fingerprints(X_index, "smiles")
X_index["IP_anion"] = X_index["IP_anion"].astype(float)
X_index["IP_neutral"] = X_index["IP_neutral"].astype(float)

features_csv = os.path.join(directory_path, "all_features.csv")
X_index.to_csv(features_csv, index=False)
print(f"Feature matrix saved to: {features_csv}")
print(X_index)

exclude_columns = ["nitrogen_index", "smiles"]
X = X_index.drop(exclude_columns, axis=1)
print(X)

# 20 Load model
model = joblib.load(MODEL_PATH)

# 21 Output prediction
pred = model.predict(X)
proba = model.predict_proba(X)

print(pred)
print(proba)

# 22 Visualisation
draw_atom_label(Chem.MolFromSmiles(substrate), X_index, np.round(proba, 2), run_name)
