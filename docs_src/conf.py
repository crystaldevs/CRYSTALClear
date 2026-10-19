# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html
#
# Build command: sphinx-build -b html ./ ../docs/
#
# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#



import os
import sys

sys.path.insert(0, os.path.abspath('.'))
sys.path.insert(0, os.path.abspath('../..'))
sys.path.insert(0, os.path.abspath('../'))
sys.path.insert(0, os.path.abspath('../CRYSTALClear'))


# -- Project information -----------------------------------------------------

from CRYSTALClear import __author__, __version__
import datetime

year = datetime.datetime.now().date().strftime("%Y")
project = 'CRYSTALClear'
copyright = year + ', ' + __author__
author = __author__

# The full version, including alpha/beta/rc tags
release = __version__


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinx.ext.autodoc', 'sphinx.ext.autosummary',
              'sphinx.ext.napoleon', 'sphinx.ext.githubpages','nbsphinx',
              # As on CRYSTALLine's site: tabs for the install routes, and a
              # copy button on every code block.
              'sphinx_design', 'sphinx_copybutton']
# Set `Returns` section to behave like the `Args` section
# For Google Doc format
napoleon_custom_sections = [('Returns', 'params_style')]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

# The suffix of source filenames.
source_suffix = '.rst'

# The master toctree document.
master_doc = 'index'

# If true, the current module name will be prepended to all description
# unit titles (such as .. function::).
add_module_names = False

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# Sort functions as they are in source code.
autodoc_member_order = 'bysource'

# \AA for Angstrom symbol in math block
rst_prolog = u"""\

.. only:: html

  :math:`\\newcommand{\\AA}{\\text{Å}}`

"""

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.

html_theme = 'sphinx_book_theme'

# Other configuration setups for html.
# https://sphinx-rtd-theme.readthedocs.io/en/latest/configuring.html

html_theme_options = {
    "repository_url": "https://github.com/crystaldevs/CRYSTALClear",
    "repository_branch": "main",
    "path_to_docs": "docs_src",
    "use_repository_button": True,
    "use_edit_page_button": True,
    "use_issues_button": True,
    # Not listed again: the home page is the introduction, which the User Guide
    # already opens with.
    "home_page_in_toc": False,
    "show_navbar_depth": 1,
    "max_navbar_depth": 3,
    # The landing page and the introduction are short, and carry a gallery
    # and a diagram: those get the width a contents column would take.
    "secondary_sidebar_items": {"**": ["page-toc"], "index": [], "introduction": []},
    # The logo's lattice is grey, which a dark page swallows: a lifted copy.
    "logo": {
        "image_light": "_static/CRYSTAL_logo.png",
        "image_dark": "_static/CRYSTAL_logo-dark.png",
    },
}

html_context = {
    # Light by default, as on CRYSTALLine's site. The theme reads this from the
    # context, not from html_theme_options, where it is rejected.
    'default_mode': 'light',
    'display_github': True,
    'github_user': 'crystaldevs',
    'github_repo': 'CRYSTALClear',
    'github_version': 'main',
    'conf_py_path': '/docs_src/',
}

# The colours shared with CRYSTALLine's site.
html_css_files = ['custom.css']

# favicon
html_favicon = '_static/favicon.ico'

# logo
html_logo = "_static/CRYSTAL_logo.png"

# title
html_title = "CRYSTALClear"


# Add any paths that contain custom themes here, relative to this directory.
html_theme_path = ['.']

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# The URL to doc site
html_baseurl = 'https://crystaldevs.github.io/CRYSTALClear/'

# -- Options for LaTeX intepreter---------------------------------------------
latex_elements = {
    'preamble': r'''
    \usepackage{amsmath}
    \usepackage{amsfonts}
    \usepackage{amssymb}
    ''',
    'inputenc': r'\usepackage[utf8]{inputenc}',
    'utf8extra': '',
    'latex_engine': 'pdflatex',
}

# -- To put jupyter notebooks into the documentation site ---------------
nbsphinx_allow_errors = True
