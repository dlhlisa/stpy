"""
This contains some helper functions.

"""

import sys 
import os
from pathlib import Path 
PROJECT_ROOT = Path(__file__).resolve().parent
print(f'Project root directory: {PROJECT_ROOT}')
sys.path.append(str(PROJECT_ROOT))
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

import requests
import logging

from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem import rdMolDescriptors, MACCSkeys 
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem import Draw
from sklearn.base import BaseEstimator, TransformerMixin


# Utility functions for molecule processing
def safe_canonicalsmi_from_smiles(smi):
    """Safely generate canonical SMILES from input SMILES string.

    Args:
        smi (string): SMILES string.

    Returns:
        string: Canonical SMILES string.

    Examples:
        >>> smiles = 'C1=CC=CC=C1OCOC'
        >>> canon_smi = safe_canonicalsmi_from_smiles(smiles)
        >>> print(canon_smi)
        COC1=CC=CC=C1O

        >>> a =['COCCCN', 'c1ccccc1OCOC', None, 'C1CCCCC1O', 'C1=CC=CC=C1', 'invalid_smiles'] 
        >>> df = pd.DataFrame({'smiles': a})
        >>> df['canonical_smi'] = df['smiles'].apply(safe_canonicalsmi_from_smiles)
        >>> print(df) 
        smiles canonical_smi
        0          COCCCN        COCCCN
        1    c1ccccc1OCOC  COCOc1ccccc1
        2            None          None
        3       C1CCCCC1O     OC1CCCCC1
        4     C1=CC=CC=C1      c1ccccc1
        5  invalid_smiles          None
    """
    try:
        canonical_smi = Chem.MolToSmiles(Chem.MolFromSmiles(smi))
        return canonical_smi
    except:
        return None


class MoleculeStandardizer:
    """
    A reusable, configurable RDKit molecule standardization pipeline.

    https://github.com/rdkit/rdkit/blob/master/Code/GraphMol/MolStandardize/TransformCatalog/normalizations.in

    Parameters
    ----------
    steps : list of str, optional
        Ordered list of standardization steps to apply.
        Supported steps:
            - "normalize"
            - "remove_fragments"
            - "reionize"
            - "tautomer_parent"
            - "stereo_parent"
            - "isotope_parent"
            - "charge_parent"
            - "super_parent"
            - "cleanup"
    largest_fragment : bool, optional
        Keep only the largest fragment after processing.
    numThreads : int, optional
        Number of threads for RDKit operations.
    logger : logging.Logger, optional
        Logger for reporting actions.

    Notes
    -----
    This class is designed for batch processing and reproducible pipelines.
    """

    def __init__(
        self,
        steps=None,
        largest_fragment=False,
        numThreads=4,
        logger=None,
    ):
        self.steps = steps or ["normalize", "remove_fragments", "reionize"]
        self.largest_fragment = largest_fragment
        self.numThreads = numThreads
        self.logger = logger or logging.getLogger(__name__)

        # Map step names to RDKit functions
        self.step_functions = {
            "normalize": rdMolStandardize.NormalizeInPlace,
            "remove_fragments": rdMolStandardize.RemoveFragmentsInPlace,
            "reionize": rdMolStandardize.ReionizeInPlace,
            "tautomer_parent": rdMolStandardize.TautomerParentInPlace,
            "stereo_parent": rdMolStandardize.StereoParentInPlace,
            "isotope_parent": rdMolStandardize.IsotopeParentInPlace,
            "charge_parent": rdMolStandardize.ChargeParentInPlace,
            "super_parent": rdMolStandardize.SuperParentInPlace,
            "cleanup": rdMolStandardize.CleanupInPlace,
        }

    # ------------------------------------------------------------
    # Core API
    # ------------------------------------------------------------
    def standardize_mols(self, mols):
        """Apply the configured standardization pipeline to RDKit Mol objects."""
        for step in self.steps:
            if step not in self.step_functions:
                raise ValueError(f"Unknown standardization step: {step}")

            self.logger.debug(f"Applying step: {step}")
            self.step_functions[step](mols, numThreads=self.numThreads)

        if self.largest_fragment:
            chooser = rdMolStandardize.LargestFragmentChooser()
            mols = [chooser.choose(m) for m in mols]
            self.logger.debug("Applied largest fragment selection")

        return mols

    # ------------------------------------------------------------
    # Public API for SMILES
    # ------------------------------------------------------------
    def standardize_smiles(self, smiles_list):
        """
        Standardize a list of SMILES strings.

        Returns
        -------
        list of str
            Standardized SMILES strings.
        """
        valid_mols = []
        invalid = []

        # Convert SMILES → Mol
        for s in smiles_list:
            mol = Chem.MolFromSmiles(s)
            if mol is None:
                invalid.append(s)
                self.logger.warning(f"Invalid SMILES skipped: {s}")
            else:
                valid_mols.append(mol)

        if invalid:
            self.logger.info(f"Skipped {len(invalid)} invalid SMILES")

        # Apply transformations
        before = len(valid_mols)
        valid_mols = self.standardize_mols(valid_mols)
        after = len(valid_mols)

        if after < before:
            self.logger.info(f"Removed {before - after} molecules during standardization")

        # Convert Mol → canonical SMILES
        return [Chem.MolToSmiles(m) for m in valid_mols]


# ----------------------------------------------------------------------
# Functional wrapper for convenience
# ----------------------------------------------------------------------
def molecule_standardize(
    smiles,
    steps=None,
    largest_fragment=False,
    numThreads=4,
    logger=None,
):
    """
    Convenience wrapper around MoleculeStandardizer.
    """
    std = MoleculeStandardizer(
        steps=steps,
        largest_fragment=largest_fragment,
        numThreads=numThreads,
        logger=logger,
    )
    return std.standardize_smiles(smiles)


# from rdkit import RDLogger
# from utils import MoleculeStandardizer
# RDLogger.DisableLog("rdApp.*")  # optional: silence RDKit warnings

# std = MoleculeStandardizer(
#     steps=["normalize", "remove_fragments", "reionize"],
#     largest_fragment=True
# )
# df = pd.DataFrame({ "compound_id": [1, 2, 3], "smiles": [ "CC.[Na+].[Cl-]", "[O-]c1[n+](C)cccc1", "invalid_smiles_here" ] })

# df["std_smiles"] = df["smiles"].apply(lambda s: std.standardize_smiles([s])[0] if Chem.MolFromSmiles(s) else None)

# -----------------------------------------------------
# Calculating fingerprint
def get_fingerprint(
    smiles,
    fp="morgan",
    radius=3,
    fpSize=2048,
    output="numpy",
    dtype="uint8",
    logger=None,
):
    """
    Memory‑optimized molecular fingerprint generator.

    Parameters
    ----------
    smiles : str
        Input SMILES string.
    fp : str
        Fingerprint type.
    output : str
        'numpy' (recommended) or 'vect'.
    dtype : str
        Data type for NumPy output: 'bool', 'uint8', or 'int8'.

    Returns
    -------
    np.ndarray or ExplicitBitVect or None
    """

    logger = logger or logging.getLogger(__name__)
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        logger.warning(f"Invalid SMILES skipped: {smiles}")
        return None

    fp = fp.lower()

    # ----------------------------
    # Special fingerprints
    # ----------------------------
    if fp == "maccs":
        bv = MACCSkeys.GenMACCSKeys(mol)
        if output == "vect":
            return bv
        arr = np.frombuffer(bv.ToBitString().encode(), "S1") == b"1"
        return arr.astype(dtype)

    if fp == "mqns":
        arr = np.array(rdMolDescriptors.MQNs_(mol))
        return arr

    # ----------------------------
    # Generator-based fingerprints
    # ----------------------------
    generators = {
        "morgan": lambda: rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=fpSize),
        "fmorgan": lambda: rdFingerprintGenerator.GetMorganGenerator(
            radius=radius,
            fpSize=fpSize,
            atomInvariantsGenerator=rdFingerprintGenerator.GetMorganFeatureAtomInvGen(),
        ),
        "rdkit": lambda: rdFingerprintGenerator.GetRDKitFPGenerator(fpSize=fpSize),
        "atompair": lambda: rdFingerprintGenerator.GetAtomPairGenerator(fpSize=fpSize),
        "toptor": lambda: rdFingerprintGenerator.GetTopologicalTorsionGenerator(fpSize=fpSize),
    }

    if fp not in generators:
        raise ValueError(f"Unknown fingerprint type: {fp}")

    gen = generators[fp]()

    if output == "vect":
        return gen.GetFingerprint(mol)

    # NumPy output (memory‑optimized)
    arr = gen.GetFingerprintAsNumPy(mol)
    return arr.astype(dtype)


def get_fingerprints_df(
    df,
    smiles_col="smiles",
    fp="morgan",
    radius=3,
    fpSize=2048,
    dtype="uint8",
):
    fps = [
        get_fingerprint(
            smi,
            fp=fp,
            radius=radius,
            fpSize=fpSize,
            output="numpy",
            dtype=dtype,
        )
        for smi in df[smiles_col]
    ]

    df[f"{fp}_fp"] = fps
    return df

# df = pd.DataFrame({"smiles": ["CCO", "c1ccccc1", "invalid"]})
# df = get_fingerprints_df(df, smiles_col="smiles", fp="morgan")
# print(df.head())


class FingerprintTransformer(BaseEstimator, TransformerMixin):
    def __init__(self, fp="morgan", radius=3, fpSize=2048, dtype="uint8"):
        self.fp = fp
        self.radius = radius
        self.fpSize = fpSize
        self.dtype = dtype

    def fit(self, X, y=None):
        return self

    def transform(self, X):
        fps = [
            get_fingerprint(
                smi,
                fp=self.fp,
                radius=self.radius,
                fpSize=self.fpSize,
                output="numpy",
                dtype=self.dtype,
            )
            for smi in X
        ]
        return np.stack(fps)


# from sklearn.pipeline import Pipeline
# from sklearn.ensemble import RandomForestClassifier
# pipe = Pipeline([
#     ("fp", FingerprintTransformer(fp="morgan")),
#     ("clf", RandomForestClassifier())
# ])
# X = df["smiles"]
# y = [0, 1, 0]
# pipe.fit(X, y)

def concat_fingerprints(
    smiles,
    fps=("morgan", "rdkit"),
    radius=3,
    fpSize=2048,
    dtype="uint8",
    logger=None,
):
    """
    Generate and concatenate multiple fingerprints for a single molecule.

    Parameters
    ----------
    smiles : str
        Input SMILES string.
    fps : tuple of str
        Fingerprint types to concatenate.
    radius : int
        Morgan radius (used where applicable).
    fpSize : int
        Bit vector size for fingerprints that support it.
    dtype : str
        Output dtype: 'bool', 'uint8', 'int8'.

    Returns
    -------
    np.ndarray or None
        Concatenated fingerprint vector.
    """

    logger = logger or logging.getLogger(__name__)
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        logger.warning(f"Invalid SMILES skipped: {smiles}")
        return None

    fp_list = []

    for fp in fps:
        arr = get_fingerprint(
            smiles,
            fp=fp,
            radius=radius,
            fpSize=fpSize,
            output="numpy",
            dtype=dtype,
            logger=logger,
        )

        if arr is None:
            logger.warning(f"Fingerprint {fp} failed for SMILES: {smiles}")
            return None

        fp_list.append(arr)

    return np.concatenate(fp_list)

# vec = concat_fingerprints(
#     "CCO",
#     fps=("morgan", "fmorgan", "maccs"),
#     fpSize=1024,
#     dtype="uint8"
# )
# print(vec.shape)


def concat_fingerprints_df(
    df,
    smiles_col="smiles",
    fps=("morgan", "rdkit"),
    radius=3,
    fpSize=2048,
    dtype="uint8",
    logger=None,
):
    """
    Compute concatenated fingerprints for all SMILES in a DataFrame.

    Returns
    -------
    np.ndarray
        Shape: (n_samples, total_fp_length)
    """

    logger = logger or logging.getLogger(__name__)

    fp_rows = []
    for smi in df[smiles_col]:
        fp = concat_fingerprints(
            smi,
            fps=fps,
            radius=radius,
            fpSize=fpSize,
            dtype=dtype,
            logger=logger,
        )
        fp_rows.append(fp)
    df[f"concat_fp"] = fp_rows
    # np.stack(fp_rows)
    return df


def fp_manipulate(df, fp1, fp2, mani='concat'):
    """
    Concatenate or sum two fingerprint columns in a DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        Input dataframe.
    fp1 : str
        First fingerprint column name.
    fp2 : str
        Second fingerprint column name.
    mani : str
        'concat' or 'sum'.

    Returns
    -------
    pd.DataFrame
        DataFrame with new column '{mani}_fp'.
    """

    # Validate columns
    if fp1 not in df.columns:
        raise KeyError(f"Fingerprint {fp1} not found in dataframe.")
    if fp2 not in df.columns:
        raise KeyError(f"Fingerprint {fp2} not found in dataframe.")

    # Drop rows with missing fingerprints
    before = len(df)
    df = df.dropna(subset=[fp1, fp2]).copy()
    removed = before - len(df)
    if removed:
        print(f"Removed {removed} rows due to missing fingerprints in {fp1} or {fp2}.")

    # Perform manipulation
    if mani == 'concat':
        df[f'{mani}_fp'] = df.apply(
            lambda row: np.concatenate(
                [np.asarray(row[fp1], dtype=np.uint8),
                 np.asarray(row[fp2], dtype=np.uint8)]
            ),
            axis=1
        )

    elif mani == 'sum':
        df[f'{mani}_fp'] = df.apply(
            lambda row: (
                np.asarray(row[fp1], dtype=np.uint8) +
                np.asarray(row[fp2], dtype=np.uint8)
            ).astype(np.uint8),
            axis=1
        )

    else:
        raise ValueError("mani must be 'concat' or 'sum'.")

    return df


def cas2smiles(cas_number='50-78-2'):
    """ 
    """
    # PubChem PUG-REST API endpoint
    api_url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{cas_number}/property/CanonicalSMILES/JSON'

    # Make the HTTP GET request
    response = requests.get(api_url)

    # Check if the request was successful
    if response.status_code == 200:
        # Parse the JSON response
        data = response.json()
        # Extract the SMILES from the response
        # print(data)
        smiles = data['PropertyTable']['Properties'][0]['ConnectivitySMILES']
        return smiles
    else:
        print(f"Error: {response.status_code}")
        return "NotFound"

# df['SMILES'] = df['CASRN'].apply(lambda x: cas2smiles(x))




