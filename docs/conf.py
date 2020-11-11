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
import os
import sys
sys.path.insert(0, os.path.abspath('..'))


# -- Project information -----------------------------------------------------

from corfu import __version__

project = 'corfu'
copyright = '2020, Nicolas Tessore'
author = 'Nicolas Tessore'

# The full version, including alpha/beta/rc tags
release = __version__


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.viewcode',
    'sphinx.ext.napoleon',
    'sphinx.ext.intersphinx',
    'sphinxcontrib.bibtex',
    'matplotlib.sphinxext.plot_directive',
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# By default, highlight as Python 3.
highlight_language = 'python3'

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'trac'

# If true, figures, tables and code-blocks are automatically numbered if they
# have a caption. The numref role is enabled. Obeyed so far only by HTML and
# LaTeX builders. Default is False.
numfig = True


# -- Options for math -------------------------------------------------------

# Set this option to True if you want all displayed math to be numbered. The
# default is False.
math_number_all = True

# If True, displayed math equations are numbered across pages when numfig is
# enabled. The numfig_secnum_depth setting is respected. The eq, not numref,
# role must be used to reference equation numbers. Default is True.
math_numfig = True


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'traditional'

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
#
html_theme_options = {
    'sidebarwidth': 320,
    'body_min_width': 0,
    'body_max_width': 800,
}

html_logo = '_static/corfu-logo-small.svg'

html_sidebars = {
    'index': [],
    '**': ['globaltoc.html', 'links.html', 'searchbox.html'],
}

html_use_index = False
html_copy_source = False
html_show_sourcelink = False

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']


# -- Intersphinx -------------------------------------------------------------

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'scipy': ('https://docs.scipy.org/doc/scipy/reference/', None),
}


#------------------------------------------------------------------------------
# Matplotlib plot_directive options
#------------------------------------------------------------------------------

plot_include_source = False
plot_formats = [('png', 96), 'pdf']
plot_html_show_formats = False
plot_html_show_source_link = True

plot_pre_code = '''
'''

plot_template = '''
{{ source_code }}

{% for img in images %}
.. figure:: {{ build_dir }}/{{ img.basename }}.*
   {% for option in options -%}
   {{ option }}
   {% endfor %}

   {{ caption }}
{% endfor %}
'''

plot_font_size = 13*72/96.0  # 13 px

plot_rcparams = {
    'font.size': plot_font_size,
    'axes.titlesize': plot_font_size,
    'axes.labelsize': plot_font_size,
    'xtick.labelsize': plot_font_size,
    'ytick.labelsize': plot_font_size,
    'legend.fontsize': plot_font_size,
    'legend.frameon': False,
    'figure.figsize': (3.2, 3.2),
    'figure.subplot.bottom': 0.2,
    'figure.subplot.left': 0.2,
    'figure.subplot.right': 0.9,
    'figure.subplot.top': 0.85,
    'figure.subplot.wspace': 0.4,
    'text.usetex': False,
}
