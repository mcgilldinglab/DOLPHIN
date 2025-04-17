# Configuration file for the Sphinx documentation builder.

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

nbsphinx_execute = 'never'
highlight_language = "python"
pygments_style = "friendly"

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']

templates_path = ['_templates']

# -- Options for HTML output

html_theme = 'sphinx_rtd_theme'

# -- Options for EPUB output
epub_show_urls = 'footnote'
