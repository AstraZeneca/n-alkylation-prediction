![Maturity level-0](https://img.shields.io/badge/Maturity%20Level-ML--0-red)
# N-Alkylation / N-Arylation Regioselectivity Predictor

A machine-learning workflow for predicting the most likely nitrogen reaction centre
in N-alkylation and N-arylation reactions of poly-nitrogen nucleophiles.

Developed by Alice Oram, AstraZeneca.
With contributions from Christoph Bauer, Isabel Arranz and Per-Ola Norrby

---

## Overview

Given a molecule SMILES as input, the predictor:

1. Validates the input (reactive-nitrogen count, molecular weight, salt, metal, charge).
2. Generates 3D coordinates with **Corina**.
3. Converts coordinates to xyz format with **OpenBabel**.
4. Deprotonates (neutral molecules) or optimises (anionic molecules) the structure
   with **CREST** / **xTB**.
5. Calculates ionisation potentials for the anionic and neutral forms with **Multiwfn**.
6. Calculates atomic descriptors (proximity shells, partial charges, polarizability)
   with **kallisto**.
7. Combines the descriptors with RDKit fingerprints and feeds them into a trained
   **XGBoost** classifier.
8. Outputs the predicted reactive nitrogen and the probability for each nitrogen centre,
   and saves a labelled molecule image (`N_alkylation_probability.png`).

---

## Requirements

### Python packages

Install the Python dependencies with:

```bash
pip install -r requirements.txt
```

### External tools

The workflow requires the following external executables to be available on your
system (see [Configuration](#configuration) for how to point the code at them):

| Tool       | Version tested | Purpose                                    |
|------------|----------------|--------------------------------------------|
| Corina     | 3.60           | SMILES to 3D mol2 coordinate generation    |
| OpenBabel  | 3.1.1          | mol2 to xyz format conversion              |
| CREST      | latest         | Deprotonation / reprotonation              |
| xTB        | latest         | Geometry optimisation and molden files     |
| Multiwfn   | 3.8 (no-GUI)   | Ionisation potential calculation           |

---

## Installation

```bash
git clone https://github.com/<your-org>/n-alkylation-arylation.git
cd n-alkylation-arylation
pip install -r requirements.txt
```

---

## Configuration

All configurable paths are set in **`config.py`** in the repository root.
Edit this file before running the predictor:

```python
# Root directory where per-run sub-folders will be created
WORKING_DIR = "/path/to/your/working/directory"

# Paths to external tool executables
MULTIWFN_EXEC = "/path/to/Multiwfn_noGUI"
CREST_EXEC    = "/path/to/crest"
XTB_EXEC      = "xtb"          # or full path if not on PATH
CORINA_EXEC   = "/path/to/corina"
OBABEL_EXEC   = "obabel"       # or full path if not on PATH
```

Alternatively, each executable path can be overridden with an environment variable
(e.g. `export MULTIWFN_EXEC=/path/to/Multiwfn_noGUI`).

The trained XGBoost model (`model/XGBoost_fp.pkl`) and the Multiwfn input file
(`util/multiwfn_IP.in`) are bundled with the repository and require no additional
configuration.

---

## Usage

```bash
python main.py
```

The script will prompt for:

- **Molecule SMILES** - the input nitrogen nucleophile.
- **Molecule name** - used as the run folder name and file prefix.

Example input:

```
Molecule SMILES: Cc1ccc(NC(=O)c2ccccc2)cc1N
Name the molecule: example_molecule
```

Output files are written to `<WORKING_DIR>/<molecule_name>/`:

| File                            | Description                                          |
|---------------------------------|------------------------------------------------------|
| `<name>.smi`                    | Input SMILES file for Corina                         |
| `<name>.mol2`                   | 3D structure (Corina)                                |
| `<name>.xyz`                    | xyz coordinates (OpenBabel)                          |
| `deprotonated.xyz` / `xtbopt.xyz` | Optimised anionic structure                        |
| `protonated.xyz`                | Neutral (re-protonated) structure                    |
| `multiwfn_IP_Deprot.out`        | Multiwfn output for the anionic form                 |
| `multiwfn_IP_Reprot.out`        | Multiwfn output for the neutral form                 |
| `all_features.csv`              | Full feature matrix used for prediction              |
| `N_alkylation_probability.png`  | Labelled molecule image with per-nitrogen probabilities |

---

## Repository structure

```
n-alkylation-arylation/
├── config.py               # User-editable configuration (paths, executables)
├── main.py                 # Main workflow script
├── requirements.txt        # Python dependencies
├── model/
│   ├── XGBoost_fp.pkl      # Trained XGBoost classifier
│   └── README.md           # Model card
├── util/
│   ├── cheminformatics.py  # Input validation and visualisation
│   ├── corina_util.py      # Corina wrapper
│   ├── fp_util.py          # RDKit fingerprint generation
│   ├── kallisto_util.py    # Kallisto descriptor calculations
│   ├── multiwfn_IP.in      # Multiwfn input commands
│   ├── multiwfn_util.py    # Multiwfn wrapper and output parser
│   ├── obabel_util.py      # OpenBabel wrapper
│   └── xtb_util.py         # xTB / CREST wrapper
└── LICENSE
```

---

## Known limitations and TODOs

- Nitrogen indices are taken from the deprotonated structure. If the atom ordering
  changes during reprotonation the merge of `IP_anion` and `IP_neutral` on
  `nitrogen_index` may fail. The function `check_atom_index_change` prints a warning
  when this occurs; a full fix is planned.
- The workflow currently stops with an `AssertionError` / `ValueError` if any
  validation check fails. Graceful error recovery is planned.
- The model was trained on molecules with MW < 500 g/mol and charge 0 or -1.
  Predictions outside this domain should be treated with caution.

---

## License

This project is licensed under the Apache 2.0 License - see the [LICENSE](LICENSE) file for details.
