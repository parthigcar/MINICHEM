# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))
# import Sphinx

# -- Project information -----------------------------------------------------
# The master toctree document.
master_doc = 'index'
project = 'MINICHEM'
copyright = '2020, Parthkumar Patel, A. John Arul'
author = 'Parthkumar Patel, A. John Arul'

# The full version, including alpha/beta/rc tags
release = '0.0.1'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinx.ext.mathjax', 'sphinx_numfig', 'sphinx.ext.autosectionlabel', 'sphinx_rtd_theme']

mathjax_config = {
    "jax": ["input/TeX","output/HTML-CSS"],
    "displayAlign": "center"
}

numfig = True
# numref = True
# numfig_number_figures = True
numfig_figure_caption_prefix = "Fig."
# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
# The theme to use for HTML and HTML Help pages
import sphinx_rtd_theme
html_theme = 'sphinx_rtd_theme'
html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]
# html_static_path = ['_static']
# Register the theme as an extension to generate a sitemap.xml

# # The name for this set of Sphinx documents.  If None, it defaults to
# # "<project> v<release> documentation".
# html_title = "MINICHEM Documentation"

# A shorter title for the navigation bar.  Default is the same as html_title.
#html_short_title = None

# The name of an image file (within the static path) to use as favicon of the
# docs.  This file should be a Windows icon file (.ico) being 16x16 or 32x32
# pixels large.
#html_favicon = None

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
# html_static_path = ['_static']
# html_static_path = []
html_css_files = [
    'css/custom.css',
]
def setup(app):
    app.add_stylesheet('theme_overrides.css')

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
# html_static_path = ['_static']
