|Logo|

.. begin-header

*********
``corfu``
*********

|GitHub| |PyPI| |Docs|

A Python library for projected angular correlation functions.

Maintainers:
    Nicolas Tessore, Lucia F. de la Bella

.. end-header


Getting Started
===============

Please refer to the `documentation`_, where you can find

- the `theory`_ behind the |corfu| package,
- `results`_ to demonstrate that the package works as intended,
- an explanation of `how to use <usage_>`_ the package, and
- the `function reference <reference_>`_.


Installation
------------

The easiest way to install the |corfu| module is via pip::

    $ pip install corfu

You can install the |corfu| module from a local working directory, such as a
repository clone::

    $ cd corfu
    $ ls
    corfu.py    setup.py    ...
    $ pip install .

Alternative, the |corfu.py|_ file is self-contained and can be copied to where
it is required::

    $ cd myproject
    $ cp ../corfu/corfu.py .
    $ python
    >>> import corfu  # works


Usage
-----

To project 3d (matter) correlation functions or power spectra along lines of
sight and obtain angular correlation functions or angular power spectra:

- Use the |corfu.ptoxi|_ function to obtain the 3d (matter) correlation
  function from the 3d (matter) power spectrum.
- Use the |corfu.uneqt|_ function to obtain the projected angular correlation
  function from the 3d (matter) correlation function.
- Use the |corfu.wtocl|_ function to obtain the angular power spectrum from the
  angular correlation function.

For further guidance, see the `documentation`_, and in particular the `usage`_
page.  For an explanation of how the computation works, see the `theory`_ page.


.. text substitutions

.. |corfu| replace:: ``corfu``


.. documentation links

.. _documentation: https://corfu.readthedocs.io/en/latest/
.. _theory: https://corfu.readthedocs.io/en/latest/theory.html
.. _results: https://corfu.readthedocs.io/en/latest/results.html
.. _usage: https://corfu.readthedocs.io/en/latest/usage.html
.. _reference: https://corfu.readthedocs.io/en/latest/reference.html


.. reference links

.. |corfu.ptoxi| replace:: ``corfu.ptoxi()``
.. _corfu.ptoxi: https://corfu.readthedocs.io/en/latest/reference.html#corfu.ptoxi

.. |corfu.uneqt| replace:: ``corfu.uneqt()``
.. _corfu.uneqt: https://corfu.readthedocs.io/en/latest/reference.html#corfu.uneqt

.. |corfu.wtocl| replace:: ``corfu.wtocl()``
.. _corfu.wtocl: https://corfu.readthedocs.io/en/latest/reference.html#corfu.wtocl


.. file links

.. |corfu.py| replace:: ``corfu.py``
.. _corfu.py: corfu.py


.. layout

.. |Logo| image:: docs/_static/corfu-logo.svg
   :alt: Logo
   :width: 200


.. begin-badges

.. |GitHub| image:: https://img.shields.io/badge/github-ntessore%2Fcorfu-lightgrey
   :target: https://github.com/ntessore/corfu
   :alt: GitHub Repository

.. |PyPI| image:: https://img.shields.io/pypi/v/corfu.svg
   :target: https://pypi.org/project/corfu
   :alt: PyPI Status

.. |Docs| image:: https://readthedocs.org/projects/corfu/badge/?version=latest
   :target: https://corfu.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status

.. end-badges
