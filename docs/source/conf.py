#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Configuration file for the Sphinx documentation builder."""

# -- Project information -----------------------------------------------------

project = 'SNe Constraints'
copyright = '2019, Daniel Perrefort'
author = 'Daniel Perrefort, Yike Zhang'
version, release = '', ''

# -- General configuration ---------------------------------------------------

extensions = [
    'sphinx.ext.autodoc',
    # 'sphinx.ext.viewcode',  # Can add back in after code is public
    'sphinx.ext.napoleon',
    'sphinxarg.ext'
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix(es) of source filenames. (source_suffix = ['.rst', '.md'])
source_suffix = '.rst'
exclude_patterns = []

# The master toctree document.
master_doc = 'index'

# -- Options for HTML output -------------------------------------------------

html_theme = 'sphinx_rtd_theme'
html_theme_options = {}
html_static_path = ['../_static']

html_context = {
    'css_files': [
        '_static/theme_overrides.css',  # override wide tables in RTD theme
    ],
}

# -- Options for HTMLHelp output ---------------------------------------------

# Output file base name for HTML help builder.
htmlhelp_basename = 'sdss_classification_doc'

# -- Options for LaTeX output ------------------------------------------------

latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    # 'papersize': 'letterpaper',

    # The font size ('10pt', '11pt' or '12pt').
    # 'pointsize': '10pt',

    # Additional stuff for the LaTeX preamble.
    # 'preamble': '',

    # Latex figure (float) alignment
    # 'figure_align': 'htbp',
}

latex_documents = [
    (master_doc, 'sdss_classification.tex', 'sdss\\fitting Documentation',
     'Daniel Perrefort, Yike Zhang', 'manual'),
]
