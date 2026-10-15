"""Sphinx configuration for the ThermOpt documentation."""

import os
import sys
from datetime import date

sys.path.insert(0, os.path.abspath("../src"))

project = "thermoptcobrapy"
author = "Pavan Kumar S"
copyright = f"{date.today().year}, BiSECt lab, IIT Madras"

try:
    from thermoptcobrapy import __version__ as release
except ImportError:  # pragma: no cover - docs build without the package
    release = "0.1.0.dev0"
version = release

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
html_title = f"{project} {release}"
html_theme_options = {
    "collapse_navigation": False,
    "navigation_depth": 3,
}

# Numpydoc-style docstrings, which is what the package uses.
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_use_param = True
napoleon_use_rtype = True

autodoc_member_order = "bysource"
autodoc_typehints = "description"
autodoc_default_options = {
    "show-inheritance": True,
}
autosummary_generate = True

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "cobra": ("https://cobrapy.readthedocs.io/en/latest/", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
}

nitpicky = False
