# Model

This directory contains the trained XGBoost model used for N-alkylation / N-arylation
regioselectivity prediction.

## File

| File              | Description                                                   |
|-------------------|---------------------------------------------------------------|
| `XGBoost_fp.pkl`  | Trained XGBoost classifier (serialised with joblib)           |

## Model details

- **Algorithm**: XGBoost gradient-boosted decision trees
- **Features**: Ionisation potentials (anion and neutral), kallisto atomic descriptors
  (proximity shells, EEQ partial charges, polarizability), and RDKit fingerprints (2048 bits)
- **Training domain**: Nitrogen nucleophiles with MW < 500 g/mol and formal charge 0 or -1
- **Output**: Binary classification (reactive / non-reactive nitrogen centre) with
  per-class probabilities

## Usage

The model is loaded automatically by `main.py` via the `MODEL_PATH` entry in `config.py`.
It should not be needed to interact with this file directly.
