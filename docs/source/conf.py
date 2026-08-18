# Configuration file for the Sphinx documentation builder.
import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

# -- Project information

project = 'DOLPHIN'
copyright = 'Ding Lab at McGill University'
author = 'Kailu Song'

# -- General configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    # 'recommonmark',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx.ext.mathjax',
    'sphinx_autodoc_typehints',
    'sphinx_copybutton',
    'nbsphinx',
    'myst_parser'
    # 'myst_nb'
]

# API pages only need signatures and docstrings; importing the full scientific
# stack would make documentation builds require the model runtime and GPU stack.
autodoc_mock_imports = [
    'intervaltree',
    'matplotlib',
    'pybedtools',
    'pyro',
    'pysam',
    'scanpy',
    'sklearn',
    'statsmodels',
    'torch',
    'torch_geometric',
    'tqdm',
]

nbsphinx_execute = 'never'
highlight_language = "python"
pygments_style = "friendly"

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']
suppress_warnings = ['misc.highlighting_failure']

templates_path = ['_templates']

# -- Options for HTML output

html_theme = 'sphinx_rtd_theme'

# -- Options for EPUB output
epub_show_urls = 'footnote'
