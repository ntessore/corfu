Usage
=====

The ``corfu`` package contains functionality to project three-dimensional
correlation functions and power spectra along lines of sight to obtain the
angular correlation function and angular power spectrum on the sphere.  The
general process is as follows.

- If you have a 3d power spectrum, such as the matter power spectrum obtained
  e.g. from CAMB, CLASS, or Eisenstein-Hu, use the :func:`corfu.ptoxi` function
  to obtain the 3d (matter) correlation function.
- If the objective is the computation of an angular power spectrum, use the
  :func:`corfu.theta` function to obtain a set of evaluation points.
- If you have a 3d (matter) correlation function and filter functions along the
  line of sight, use the :func:`corfu.xitow` function to obtain the projected
  angular correlation function.
- If you have the angular correlation function, use the :func:`corfu.wtocl`
  function to obtain the angular power spectrum.
