# Configuration file for the Sphinx documentation builder.

import os
import sys
sys.path.insert(0, os.path.abspath('..'))

project = 'π-Stack Optimizer'
copyright = '2025, Arunima Ghosh, Susmita Barik, Roshan J. Singh, Sandeep K. Reddy'
author = 'Arunima Ghosh, Susmita Barik, Roshan J. Singh, Sandeep K. Reddy'
release = '1.0'

# -- General configuration ---------------------------------------------------

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',  # Support for Google/NumPy style docstrings
    'sphinx.ext.viewcode',
    'sphinx.ext.todo',
    'sphinx_design',
    'sphinx_copybutton',
    'sphinx_togglebutton',
    'sphinx.ext.mathjax',
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# -- Options for HTML output -------------------------------------------------

html_theme = 'sphinx_rtd_theme'
html_css_files = ["custom.css"]
html_static_path = ['_static']

html_context = {
    "display_github": True,
    "github_user": "sandeepgroup",
    "github_repo": "pi-stack-optimizer",
    "github_version": "main",
    "conf_py_path": "/docs/",
}

# -- Options for autodoc -----------------------------------------------------

autodoc_mock_imports = ['numpy', 'scipy', 'ase']  # Mock imports for RTD builds if needed
autoclass_content = 'both'  # Include both class docstring and __init__
