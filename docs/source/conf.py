# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "VRS-Python"
copyright = "VRS-Python Maintainers"  # TODO ???
author = "VRS-Python Maintainers"  # TODO ???
html_title = "vrs_python"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx_rtd_theme",
    "sphinx.ext.autodoc",
    "sphinx_autodoc_typehints",
    "sphinx.ext.linkcode",
    "sphinx_copybutton",
    "sphinx.ext.autosummary",
    "sphinx_github_changelog",
    "sphinx_tabs.tabs",
]

templates_path = ["_templates"]
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"

# -- autodoc things ----------------------------------------------------------
import os
import sys

sys.path.insert(0, os.path.abspath("../../"))
autodoc_preserve_defaults = True

# -- get version -------------------------------------------------------------
from importlib.metadata import version

release = version = version("ga4gh.vrs")

# -- linkcode ----------------------------------------------------------------
def linkcode_resolve(domain, info):
    if domain != "py":
        return None
    if not info["module"]:
        return None
    filename = info["module"].replace(".", "/")
    return f"https://github.com/ga4gh/vrs-python/blob/main/src/{filename}.py"


# -- code block style --------------------------------------------------------
pygments_style = "default"
pygements_dark_style = "monokai"
