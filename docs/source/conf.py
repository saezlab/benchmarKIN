# Configuration file for the Sphinx documentation builder.

# -- Project information

project = 'benchmarKIN'
copyright = '2024, Müller-Dott S., Jaehnig E.J., Pham Munchic K., Yaron-Barir T.M., Dugourd A., Garrido-Rodriguez M., Johnson J.L., Lussana A., Petsalaki E., Savage S.R., Jiang W., Lei J.T., Krug K., Cantley L.C., Mani D.R., Zhang B., Saez-Rodriguez J.'
author = 'Sophia Müller-Dott'

release = '0.1'
version = '0.1.0'

# -- General configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    'myst_parser'
]

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
