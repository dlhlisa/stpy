Getting Started with stpy
=================

Welcome to the documentation for the `stpy` package! This guide will help you get started with installing and using `stpy` for your data science and cheminformatics projects.  

Installation
------------
You can install `stpy` using pip. Run the following command in your terminal:

.. code-block:: bash

   pip install stpy 

This will install the latest version of `stpy` along with its dependencies.

Basic Usage
-----------
Once you have installed `stpy`, you can start using it in your Python scripts or Jupyter notebooks. Here is a simple example of how to use some of the utilities provided by `stpy`:

.. code-block:: python  

   from stpy.utils import canonicalize_smiles

   smiles = "C1=CC=CC=C1"  # Benzene
   canonical_smiles = canonicalize_smiles(smiles)
   print(f"Canonical SMILES: {canonical_smiles}")

This code imports the `canonicalize_smiles` function from the `stpy.utils` module, canonicalizes a SMILES string for benzene, and prints the result.

Documentation
-------------
For more detailed information on the various modules and functions available in `stpy`, please refer to the full documentation. You can find it at:  
`stpy Documentation <https://your-documentation-link-here>`_
 
Here, you will find comprehensive guides, API references, and examples to help you make the most of `stpy`.

Support
-------
If you encounter any issues or have questions about using `stpy`, feel free to reach out via the GitHub repository:  
`stpy GitHub Repository <https://github.com/lihua-notes/stpy.github>`_ or open an issue for assistance.    

We hope you find `stpy` useful for your projects! Happy coding!
