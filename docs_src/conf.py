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


# -- Release notes, taken from the GitHub releases ---------------------------
# As on CRYSTALLine's site: the notes are written once, on the release, and
# there is no second copy here for someone to remember to update. The page is
# built from the API each time the site is, so every push to main (every merged
# pull request) brings it up to date. A build with no network still succeeds:
# the page then says where the notes are instead of carrying them.
_RELEASES_API = "https://api.github.com/repos/crystaldevs/CRYSTALClear/releases?per_page=100"
_RELEASES_PAGE = "https://github.com/crystaldevs/CRYSTALClear/releases"
_RELEASES_FILE = "release-notes.rst"


def _fetch_releases(timeout=15.0):
    """The published releases, newest first, or ``None`` if they can't be read."""
    import json
    import urllib.request

    request = urllib.request.Request(
        _RELEASES_API,
        headers={"Accept": "application/vnd.github+json",
                 "User-Agent": "CRYSTALClear-docs"},
    )
    # CI passes its token, which has a much higher rate limit; locally there is
    # none, and sixty requests an hour is plenty.
    token = os.environ.get("GITHUB_TOKEN")
    if token:
        request.add_header("Authorization", "Bearer " + token)
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            return json.load(response)
    except Exception:  # offline, rate-limited, API changed: same answer
        return None


def _inline(text):
    """GitHub's inline Markdown as reStructuredText: links, code, italics."""
    import re

    text = re.sub(r"\[([^\]]+)\]\((https?://[^)]+)\)", r"`\1 <\2>`__", text)
    text = re.sub(r"`([^`]+)`", r"``\1``", text)
    return re.sub(r"(?<![\w*])_(\S.*?\S)_(?![\w])", r"*\1*", text)


def _release_body(markdown):
    """One release's notes as reStructuredText lines, the text as written."""
    import re

    lines = []
    for raw in markdown.replace("\r\n", "\n").split("\n"):
        line = raw.rstrip()
        heading = re.match(r"#{1,6}\s+(.*)", line)
        bullet = re.match(r"[*-]\s+(.*)", line)
        if heading:
            # A rubric, not a section: the page's sections are the releases.
            lines += ["", ".. rubric:: " + heading.group(1), ""]
        elif bullet:
            if lines and lines[-1] and not lines[-1].startswith("* "):
                lines.append("")
            lines.append("* " + _inline(bullet.group(1)))
        else:
            if line and lines and lines[-1].startswith("* "):
                lines.append("")
            lines.append(_inline(line))
    tidy = []
    for line in lines:
        if line or (tidy and tidy[-1]):
            tidy.append(line)
    while tidy and not tidy[-1]:
        tidy.pop()
    return tidy


def _release_notes_page(releases):
    """The page, from what the API gave us, or the fallback when it gave nothing."""
    lines = ["Release notes", "=============", ""]
    if not releases:
        return "\n".join(lines + [
            "The notes for each release are published on GitHub:",
            "",
            _RELEASES_PAGE,
            "",
            "(They are normally reproduced here, but the releases could not be "
            "read when this page was built.)",
            "",
        ])
    lines += ["Every published release, newest first, as written on "
              "`GitHub <" + _RELEASES_PAGE + ">`__.", ""]
    for release in releases:
        if release.get("draft"):
            continue
        tag = release.get("tag_name") or "untagged"
        name = (release.get("name") or "").strip()
        # The tags are dates; the name is where the version people install is.
        title = "{} ({})".format(name, tag) if name and name != tag else tag
        published = (release.get("published_at") or "")[:10]
        if published:
            title += " — " + published
        if release.get("prerelease"):
            title += " (pre-release)"
        body = _release_body(release.get("body") or "")
        lines += [title, "-" * len(title), ""]
        lines += (body or ["*No notes were published with this release.*"]) + [""]
    return "\n".join(lines)


def _write_release_notes(app):
    page = _release_notes_page(_fetch_releases())
    with open(os.path.join(app.srcdir, _RELEASES_FILE), "w", encoding="utf-8") as handle:
        handle.write(page.rstrip() + "\n")


def setup(app):
    app.connect("builder-inited", _write_release_notes)
