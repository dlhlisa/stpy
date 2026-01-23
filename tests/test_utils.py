import numpy as np
import pandas as pd
import pytest
from rdkit import Chem
from stpy.utils import (
    FingerprintTransformer,
    MoleculeStandardizer,
    cas2smiles,
    concat_fingerprints,
    concat_fingerprints_df,
    fp_manipulate,
    get_fingerprint,
    get_fingerprints_df,
    molecule_standardize,
    safe_canonicalsmi_from_smiles,
)


class TestSafeCanonicalSMILES:
    def test_valid_smiles(self):
        smiles = "C1=CC=CC=C1OCOC"
        result = safe_canonicalsmi_from_smiles(smiles)
        assert result == "COCOc1ccccc1"

    def test_invalid_smiles(self):
        result = safe_canonicalsmi_from_smiles("invalid_smiles")
        assert result is None

    def test_none_input(self):
        result = safe_canonicalsmi_from_smiles(None)
        assert result is None


class TestMoleculeStandardizer:
    def test_standardize_smiles(self):
        std = MoleculeStandardizer(steps=["normalize"])
        smiles_list = ["CC.[Na+].[Cl-]"]
        result = std.standardize_smiles(smiles_list)
        assert len(result) == 1
        assert Chem.MolFromSmiles(result[0]) is not None

    def test_invalid_smiles(self):
        std = MoleculeStandardizer()
        smiles_list = ["invalid", "CCO"]
        result = std.standardize_smiles(smiles_list)
        assert len(result) == 1
        assert result[0] == "CCO"


class TestMoleculeStandardizeFunction:
    def test_standardize_function(self):
        smiles = ["CC.[Na+].[Cl-]"]
        result = molecule_standardize(smiles)
        assert len(result) == 1
        assert Chem.MolFromSmiles(result[0]) is not None


class TestGetFingerprint:
    def test_morgan_fingerprint(self):
        smiles = "CCO"
        fp = get_fingerprint(smiles, fp="morgan", output="numpy")
        assert isinstance(fp, np.ndarray)
        assert fp.shape == (2048,)

    def test_invalid_smiles(self):
        fp = get_fingerprint("invalid", fp="morgan")
        assert fp is None

    def test_maccs_fingerprint(self):
        smiles = "CCO"
        fp = get_fingerprint(smiles, fp="maccs", output="numpy")
        assert isinstance(fp, np.ndarray)
        assert fp.shape == (167,)


class TestGetFingerprintsDF:
    def test_add_fingerprints(self):
        df = pd.DataFrame({"smiles": ["CCO", "c1ccccc1"]})
        result = get_fingerprints_df(df, smiles_col="smiles", fp="morgan")
        assert "morgan_fp" in result.columns
        assert len(result) == 2
        assert isinstance(result["morgan_fp"].iloc[0], np.ndarray)


class TestFingerprintTransformer:
    def test_fit_transform(self):
        transformer = FingerprintTransformer(fp="morgan")
        X = ["CCO", "c1ccccc1"]
        transformer.fit(X)
        result = transformer.transform(X)
        assert isinstance(result, np.ndarray)
        assert result.shape == (2, 2048)


class TestConcatFingerprints:
    def test_concat(self):
        smiles = "CCO"
        fp = concat_fingerprints(smiles, fps=("morgan", "rdkit"))
        assert isinstance(fp, np.ndarray)
        assert fp.shape == (2048 + 2048,)  # default fpSize

    def test_invalid_smiles(self):
        fp = concat_fingerprints("invalid", fps=("morgan",))
        assert fp is None


class TestConcatFingerprintsDF:
    def test_concat_df(self):
        df = pd.DataFrame({"smiles": ["CCO"]})
        result = concat_fingerprints_df(df, fps=("morgan", "rdkit"))
        assert "concat_fp" in result.columns
        assert len(result) == 1
        assert isinstance(result["concat_fp"].iloc[0], np.ndarray)


class TestFPManipulate:
    def test_concat_manipulate(self):
        df = pd.DataFrame(
            {
                "fp1": [
                    np.array([1, 0, 1], dtype=np.uint8),
                    np.array([0, 1, 0], dtype=np.uint8),
                ],
                "fp2": [
                    np.array([1, 1, 0], dtype=np.uint8),
                    np.array([1, 0, 1], dtype=np.uint8),
                ],
            }
        )
        result = fp_manipulate(df, "fp1", "fp2", mani="concat")
        assert "concat_fp" in result.columns
        assert len(result["concat_fp"].iloc[0]) == 6

    def test_sum_manipulate(self):
        df = pd.DataFrame(
            {
                "fp1": [
                    np.array([1, 0, 1], dtype=np.uint8),
                    np.array([0, 1, 0], dtype=np.uint8),
                ],
                "fp2": [
                    np.array([1, 1, 0], dtype=np.uint8),
                    np.array([1, 0, 1], dtype=np.uint8),
                ],
            }
        )
        result = fp_manipulate(df, "fp1", "fp2", mani="sum")
        assert "sum_fp" in result.columns
        assert result["sum_fp"].iloc[0][0] == 2  # 1+1

    def test_missing_column(self):
        df = pd.DataFrame({"fp1": [np.array([1])]})
        with pytest.raises(KeyError):
            fp_manipulate(df, "fp1", "fp2", mani="concat")


class TestCas2Smiles:
    def test_valid_cas(self):
        # Note: This test may fail if API is down or CAS is invalid
        # Using a known CAS for aspirin: 50-78-2
        result = cas2smiles("50-78-2")
        assert isinstance(result, str)
        assert result != "NotFound"

    def test_invalid_cas(self):
        result = cas2smiles("invalid-cas")
        assert result == "NotFound"
