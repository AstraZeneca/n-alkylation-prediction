"""
config.py - Configuration for the N-alkylation/arylation predictor.

Before running the predictor, set the paths below to match the locations
of the external tools and working directory on your system.
"""

import os

# ---------------------------------------------------------------------------
# Working directory
# ---------------------------------------------------------------------------
# Root directory where per-run sub-folders will be created.
# Default: a directory named "runs" next to this config file.
WORKING_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "runs")

# ---------------------------------------------------------------------------
# External tool executables
# ---------------------------------------------------------------------------
# Path to the Multiwfn (no-GUI) executable.
MULTIWFN_EXEC = os.environ.get("MULTIWFN_EXEC", "Multiwfn_noGUI")

# Path to the CREST executable.
CREST_EXEC = os.environ.get("CREST_EXEC", "crest")

# Path to the xTB executable.
XTB_EXEC = os.environ.get("XTB_EXEC", "xtb")

# Path to the Corina executable.
CORINA_EXEC = os.environ.get("CORINA_EXEC", "corina")

# Path to the OpenBabel executable.
OBABEL_EXEC = os.environ.get("OBABEL_EXEC", "obabel")

# ---------------------------------------------------------------------------
# Multiwfn input file
# ---------------------------------------------------------------------------
# Path to the Multiwfn input file used for ionisation potential calculations.
MULTIWFN_IP_IN = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "util", "multiwfn_IP.in"
)

# ---------------------------------------------------------------------------
# Model
# ---------------------------------------------------------------------------
# Path to the trained XGBoost model file.
MODEL_PATH = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "model", "XGBoost_fp.pkl"
)
