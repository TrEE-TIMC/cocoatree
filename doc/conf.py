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
import sys
# sys.path.insert(0, os.path.abspath('.'))
import sphinx_gallery



# -- Project information -----------------------------------------------------

project = 'cocoatree'
copyright = '2024, compbio@TrEE'
author = 'compbio@TrEE'

# The full version, including alpha/beta/rc tags
release = '0.0'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx_gallery.gen_gallery",
    "sphinx.ext.autodoc",
    "sphinx.ext.mathjax",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.doctest",
    "sphinx.ext.intersphinx",
    "sphinx.ext.imgconverter",
]

mathjax_path = 'https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML'

autosummary_generate = True
autodoc_default_options = {"members": True, "inherited-members": True}

# options for math equations


# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'nameko'
html_theme_path = ["themes"]

html_theme_options = {
    # 'logo': 'logo.png',
    #'github_user': 'tree-timc',
    #'github_repo': 'cocoatree',
}


# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = []

# exclude_patterns = ['modules/generated']

# Sphinx-gallery options
from sphinx_gallery.sorting import ExplicitOrder, FileNameSortKey
examples_dirs = ['../examples']

intersphinx_mapping = {
    'python': ('https://docs.python.org/{.major}'.format(sys.version_info), None),
    'numpy': ('https://docs.scipy.org/doc/numpy/', None),
    'scipy': ('https://docs.scipy.org/doc/scipy/reference', None),
    'matplotlib': ('https://matplotlib.org/', None),
}   

sphinx_gallery_conf = {
    'image_scrapers': ("matplotlib", "cocoatree._scraper.png_scraper"),
    'backreferences_dir': "generated/backreferences",
    'doc_module': ("cocoatree", ),
    'abort_on_example_error': True,
    "reference_url": {"cocoatree": None},
    'binder': {
        # Required keys
        'org': 'tree-timc',
        'repo': 'cocoatree',
        'branch': 'gh-pages',
        'binderhub_url': 'https://mybinder.org',  # Any URL of a binderhub deployment. Must be full URL (e.g. https://mybinder.org).
        'dependencies': '../requirements.txt',
        'use_jupyter_lab': True
        }
}

