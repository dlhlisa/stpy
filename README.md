# Small Tricks / Reusable Scripts

[![codecov](https://codecov.io/github/lihua-notes/stpy/graph/badge.svg?token=IFAUQPJU3U)](https://codecov.io/github/lihua-notes/stpy)
[![Documentation Status](https://readthedocs.org/projects/stpy-master/badge/?version=latest)](https://stpy-master.readthedocs.io/en/latest/?badge=latest)
[![PyPI version](https://img.shields.io/pypi/v/stpy.svg)](https://pypi.org/project/stpy/)
![CI](https://github.com/lihua-notes/stpy/actions/workflows/python-ci.yml/badge.svg)

Include general rules, tips, tricks, tools and others for my daily use.

## Getting Started with stpy

Welcome to the stpy package! This package is aimed to help in your data science and cheminformatics projects.

### Installation

You can install stpy using pip. Run the following command in your terminal:

```bash
pip install stpy
```

This will install the latest version of stpy along with its dependencies.

### Basic Usage

Once you have installed stpy, you can start using it in your Python scripts or Jupyter notebooks. Here is a simple example of how to use some of the utilities provided by stpy:

```python
from stpy.utils import canonicalize_smiles

smiles = "C1=CC=CC=C1"  # Benzene
canonical_smiles = canonicalize_smiles(smiles)
print(f"Canonical SMILES: {canonical_smiles}")
```

This code imports the canonicalize_smiles function from the `stpy.utils` module, canonicalizes a SMILES string for benzene, and prints the result.

For more detailed information on the various modules and functions available in `stpy`, please refer to the full documentation: [stpy Documentation](https://stpy-master.readthedocs.io/en/latest/getting_started.html). You will find comprehensive guides, API references, and examples to help you make the most of stpy.

### Support

If you encounter any issues or have questions about using stpy, feel free to reach out or open an issue for assistance.

We hope you find stpy useful for your projects! Happy coding!

#### For developers/maintainers

After updating your virtual environment, update the requirements.txt and env.yaml accordingly.
- `pip list --format=freeze > requirements.txt`
- `conda env export > env.yaml`

### Copyright

Copyright (c) 2026, Lihua Deng