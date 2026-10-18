# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'ecRadOff'
copyright = '2026, Paolo Andreozzi'
author = 'Paolo Andreozzi'
release = '0.1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration
templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_logo = '_static/logo_ecradoff.svg'
html_theme = 'sphinx_rtd_theme'
html_theme_options = {
    "prev_next_buttons_location": "bottom",
    "logo_only": True,
}
html_css_files = [
    "custom.css"
]

html_static_path = ['_static']

extensions = [
        "sphinx.ext.autodoc",
        "sphinx.ext.napoleon",
        "myst_parser"
]

myst_enable_extensions = [
    "dollarmath",
]

import os
import sys

sys.path.insert(0, os.path.abspath("../src"))
