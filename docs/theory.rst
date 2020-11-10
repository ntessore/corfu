Theory
======

Authors:
    Nicolas Tessore, Lucia F. de la Bella


Projections in cosmic lensing
-----------------------------

Our aim is to compute the angular correlation function :math:`w(\theta)` of the
cosmic convergence field,

.. math::
   :label: w-kappa

    w(\theta)
    = \frac{9H_0^4\Omega_m^2}{4c^4}
        \iint_{0}^{\infty} \! \frac{q_1(x_1) \, x_1}{a(x_1)} \,
            \frac{q_2(x_2) \, x_2}{a(x_2)} \,
                \xi(r_{12}; t_1, t_2) \, dx_1 \, dx_2 \;,

where :math:`H_0` and :math:`\Omega_m` are the cosmological parameters,
:math:`x` is comoving distance, :math:`a(x)` is the scale factor corresponding
to :math:`x`, :math:`q_1(x)` and :math:`q_2(x)` are the lensing efficiencies
given by

.. math::
   :label: q

    q_i(x)
    = \int_{x}^{\infty} \! \frac{x' - x}{x'} \, n_i(x') \, dx'

for the distributions :math:`n_1(x)` and :math:`n_2(x)` of observed sources, and
:math:`\xi(r_{12}; t_1, t_2)` is the unequal-time matter correlation function
for the separation

.. math::

    r_{12}
    = \sqrt{x_1^2 + x_2^2 - 2x_1x_2 \cos\theta}

at cosmic times :math:`t_1` and :math:`t_2` corresponding to :math:`x_1` and
:math:`x_2`.  For simplicity, we have assumed a flat universe in writing
:eq:`w-kappa` and :eq:`q`; however, everything that follows applies identically
to nonflat geometries.

Other important two-point functions in weak lensing are galaxy clustering and galaxy-galaxy lensing.
Galaxy clustering quantifies correlations between galaxy number density fields

.. math::
   :label: w_gk
    
     w_{gg}(\theta)
         = \frac{1}{c^2} \iint_{0}^{\infty}  \! b(x_1) n_1(x_1)H(x_1) \, b(x_2) n_2(x_2)H(x_2) \, \xi(r_{12}; t_1, t_2) \, dx_1 \, dx_2 \; ,

where :math:`b(x)` is the bias parameter and :math:`H(x)` the Hubble parameter at a redshift corresponding to a comoving distance :math:`x`.


Galaxy-galaxy lensing quantifies the correlation between the shape of background (or source) and foreground (lens) galaxy number density. In the weak lensing regime, the observed galaxy shape is the sum of an intrinsic (unlensed component) and a shear due to gravitational lensing. For simplicity, we only consider the shear component:

.. math::
   :label: w_gk
    
     w_{\kappa g}(\theta)
        = \frac{3H_0^2\Omega_m}{2c^3} \iint_{0}^{\infty} \! \frac{q_1(x_1) \, x_1}{a(x_1)} \, b(x_2) n_2(x_2)H(x_2) \, \xi(r_{12}; t_1, t_2) dx_1 \, dx_2 \;. 


Generic form:

.. math::
   :label: w

    w(\theta)
    = \iint_{0}^{\infty} \! f_1(x_1) \, f_2(x_2) \,
        \xi(r_{12}; t_1, t_2) \, dx_1 \, dx_2 \;.

Transformation:

.. math::
   :label: ptoxi

    \xi(r; t_1, t_2)
    = \frac{1}{2\pi^2} \int_{0}^{\infty} \! P(k; t_1, t_2) \,
                \frac{\sin kr}{kr} \, k^2 \, dk

Addition theorem for spherical Bessel functions (:cite:`1972hmfw_book_____A`,
10.1.45):

.. math::
   :label: addthm

    \frac{\sin kr_{12}}{kr_{12}}
    = \sum_{l} (2l + 1) \, j_l(kx_1) \, j_l(kx_2) \, P_l(\cos\theta) \;.

Inserting :eq:`ptoxi` and :eq:`addthm` into :eq:`w`, and exchanging the order of
summation and integration, yields the relation

.. math::
   :label: cltow

    w(\theta)
    = \sum_{l} \frac{2l + 1}{4\pi} \, C_l \, P_l(\cos\theta)

between angular correlation function :math:`w(\theta)` and the angular power
spectrum

.. math::
   :label: cl

    C_l
    = \frac{2}{\pi} \iiint_{0}^{\infty} \! f_1(x_1) \, f_2(x_2) \,
        P(k; t_1, t_2) \, j_l(kx_1) \, j_l(kx_2) \, k^2 \,
            dk \, dx_1 \, dx_2 \;.

In practice, it is not cosmic convergence but cosmic shear that is observable.
The two-point statistics are related through their respective angular power
spectra :math:`C_l^{\kappa\kappa}` and :math:`C_l^{\gamma\gamma}`, with

.. math::

    C_l^{\gamma\gamma}
    = \frac{(l-1) \, (l+2)}{l \, (l+1)} \, C_l^{\kappa\kappa} \;.

The results we obtain for cosmic convergence are therefore readily applied to
cosmic shear.


How to compute the projection
-----------------------------

Evaluating the angular power spectrum is
far from straightforward. The reason is that the integral :eq:`cl` contains the
product of two highly oscillatory spherical Bessel functions :math:`j_l(kx_1) \,
j_l(kx_2)`. However, this is not dissimilar to the well-known problem of
integrating a spherical function against a pair of highly oscillatory spherical
harmonics :math:`Y_{lm} \, Y_{l'm'}`.  Encouragingly, the integral over the
sphere is routinely evaluated, by first expanding the function to be integrated
into spherical harmonics, and subsequently using Gaunt's triple integral,

.. math::

    \int_{S^2} \! Y_{l_1m_1}(\hat{n}) \, Y_{l_2m_2}(\hat{n}) \,
        Y_{l_3m_3}(\hat{n}) \, d\hat{n} = Y_{l_1l_2l_3m_1m_2m_3} \;,

where the Gaunt coefficient :math:`Y_{l_1l_2l_3m_1m_2m_3}` is readily computed.
We posit that a similar result would represent the most useful analytical
solution for the angular power spectrum integral :eq:`cl`.  It turns out
that we can derive the general form of such a solution, if it exists, without
performing any actual calculations.

To this end, let us assume for a moment that there exists a set of basis
functions :math:`\tilde{\jmath}_{l'}(k)` such that *i)* the unequal-time
power spectrum can be expanded in this basis,

 .. math::
   :label: Pk-j

    P(k; t_1, t_2)
    = \sum_{l'} a_{l'}(t_1, t_2) \, \tilde{\jmath}_{l'}(k) \;,

where the modes :math:`a_{l'}(t_1, t_2)` of the expansion are necessarily
time-dependent, and *ii)* that the overlap integral between two spherical
Bessel functions and the basis functions reduces to coefficients
:math:`\mathfrak{J}_{ll'}` which can be evaluated,

.. math::
   :label: jjj

    \int_{0}^{\infty} \! j_l(kx_1) \, j_l(kx_2) \,
                                        \tilde{\jmath}_{l'}(k) \, k^2 \, dk
    = \frac{\pi}{2} \, \mathfrak{J}_{ll'}(x_1, x_2) \;.

Since :math:`x_1` and :math:`x_2` appear on the left-hand side of :eq:`jjj` as
independent variables, the coefficients on the right-hand side must necessarily
be functions :math:`\mathfrak{J}_{ll'}(x_1, x_2)`.  Let us finally assume the
best-case scenario in which the matrix of functions :math:`\mathfrak{J}_{ll'}`
is diagonal.  Putting together our hypothetical analytical solution :eq:`Pk-j`
and :eq:`jjj`, we find that the angular power spectrum reduces at most to an
integral

.. math::
   :label: cl-ana

    C_l
    = \iint_{0}^{\infty} \! f_1(x_1) \, f_2(x_2) \, a_l(t_1, t_2) \,
                                \mathfrak{J}_{ll}(x_1, x_2) \, dx_1 \, dx_2 \;,

which is precisely the same form as :eq:`w`, the unequal-time angular
correlation function!

This short exercise shows that even if there was a convenient analytical
solution, similar to the case of the spherical harmonics, the remaining integral
:eq:`cl-ana` would still be at least as difficult to compute as the unequal-time
angular correlation function.


Exact projection in real space
------------------------------

The preceding section has shown that the most convenient approach to the exact
unequal-time projection of angular correlations is via the real-space integral,

.. math::
   :label: w2

    w(\theta)
    = \iint_{0}^{\infty} \! f_1(x_1) \, f_2(x_2) \,
                                    \xi(r_{12}; t_1, t_2) \, dx_1 \, dx_2 \;.

Evaluating this double integral should, at least in principle, be
straightforward. The filter functions :math:`f_1` and :math:`f_2` are determined
by observations, and thus known over a fixed grid of points :math:`x_1` and
:math:`x_2`. The filter grid hence provides a natural resolution for numerical
integration and, if the remaining factor :math:`\xi(r_{12}; t_1, t_2)` can be
computed, the value of :eq:`w2` is found by any suitable quadrature scheme, e.g.
the trapezoidal rule if the filters are finely enough resolved.  In practice,
this is exactly how we perform the integration, with one additional complication
that we describe shortly.

But first, we compute the matter correlation function :math:`\xi(r)` from the
matter power spectrum :math:`P(k)`, which usually is the fundamental input.
This is done by expressing the integral relationship :eq:`ptoxi` between the two
functions in terms of the Bessel function :math:`J_{1/2}`,

.. math::
   :label: ptoxi_exact

    \xi(r; t_1, t_2)
    = \frac{1}{(2\pi)^{3/2}} \int_{0}^{\infty} \! P(k; t_1, t_2) \,
                                \frac{J_{1/2}(kr)}{\sqrt{kr}} \, k^2 \, dk \;.

Integrals of this form can be efficiently evaluated over a logarithmic range of
:math:`r` values with the FFTLog algorithm :cite:`2000MNRAS_312__257H`.

Having obtained the function :math:`\xi(r; t_1, t_2)`, we need to integrate it
against the filter functions :math:`f_1(x_1)` and :math:`f_2(x_2)`.  These are
usually obtained, either directly or indirectly, from observations, and thus
given on a fixed grid of :math:`x_1` and :math:`x_2` values.

At first sight, this seems an easy proposition: The correlations :math:`\xi(r;
t_1, t_2)` change only slowly with :math:`t_1` and :math:`t_2`, and are readily
interpolated to the values :math:`t(x_1)` and :math:`t(x_2)` of the filter grid.

To understand how numerical issues arise when our various functions are defined
on grids of either :math:`x_1`, :math:`x_2`, or :math:`r`, we construct a change
of variables from the filter grid :math:`x_1, x_2` to a polar coordinate system
where :math:`r = \sqrt{x_1^2 + x_2^2 - 2x_1x_2 \cos\theta}` is the radial
coordinate.  To find the angular coordinate, we only have to write :math:`r^2`
as a sum of squares; a symmetric choice is

.. math::

    r^2
    = \biggl\{(x_1 - x_2) \, \sqrt{\frac{1 + \cos\theta}{2}}\biggr\}^2
    + \biggl\{(x_1 + x_2) \, \sqrt{\frac{1 - \cos\theta}{2}}\biggr\}^2 \;.

We therefore introduce the angle :math:`\alpha` as :math:`\tan(\alpha) =
\sqrt{\frac{1 - \cos\theta}{1 + \cos\theta}} \frac{x_1 + x_2}{x_1 - x_2}`.
Conversely, the distances :math:`x_1, x_2` for given :math:`r, \alpha` are

.. math::

    x_{1,2}
    = \sqrt{\frac{1 + \cos(\theta)}{2}} \, \frac{r \sin\alpha}{\sin\theta}
    \pm \sqrt{\frac{1 - \cos\theta}{2}} \, \frac{r \cos\alpha}{\sin\theta} \;.


.. _fig_exact-grid:
.. figure:: figures/exact-grid.*
   :alt: integration grids

   The different grids for the exact integration.


Limber's Approximation
----------------------

Write the angular correlation function :eq:`w` as the integral over the mean
radial distance :math:`x = (x_1 + x_2)/2` and the radial separation :math:`R =
x_1 - x_2`:

.. math::
   :label: w-limber-variables

    w(\theta)
    = \int_{0}^{\infty} \! \int_{-2x}^{2x} \! f_1(x+R/2) \, f_2(x-R/2) \, \xi(r_{12}; t_1, t_2) \, dR \, dx \;,

where the distance between the points in terms of :math:`x` and :math:`r` is now

.. math::

    r_{12}
    = \sqrt{2 x^2 \, (1-\cos\theta) + R^2 \, (1 + \cos\theta)/2} \;.

Limber :cite:`1953ApJ___117__134L,1954ApJ___119__655L` introduced an
approximation for the integral :eq:`w-limber-variables` using the assumption
*i)* that the filters and correlation function change slowly and can be
approximated by their midpoint values,

.. math::

    f_1(x+R/2) \, f_2(x-R/2) \, \xi(r_{12}; t_1, t_2)
    \approx f_1(x) \, f_2(x) \, \xi(r_{12}; t) \;,

where :math:`t = t(x)`; *ii)* that the angle :math:`\theta` between the
points is small, :math:`\theta \ll 1`, so that the distance :math:`r_{12}` can
be approximated as

.. math::

    r_{12}
    \approx \sqrt{x^2\theta^2 + R^2} \;;

and *iii)* that the integral over :math:`R` can be extended over the entire
real line.  Limber's approximation for the correlation function
:eq:`w-limber-variables` is thus

.. math::
   :label: w-limber

    w_{\rm L}(\theta)
    = \int_{0}^{\infty} \! f_1(x) \, f_2(x) \, \xi_{\rm L}(x\theta; t)
                                                        \, x\theta \, dx \;,

where :math:`\xi_{\rm L}` is Limber's matter correlation function, defined as

.. math::
   :label: xi_limber

    \xi_{\rm L}(r; t) =
    \frac{1}{r} \, \int_{-\infty}^{\infty} \! \xi(\sqrt{r^2 + R^2}; t) \, dR \;.

If :math:`\xi_{\rm L}` is known, the angular correlation function :eq:`w-limber`
is a single integral over a slowly changing combined filter function
:math:`f_{12}(x) = f_1(x) \, f_2(x)`.  If the filter is determined by
observations, the resolution of the integral fixed, and it can be evaluated by
any suitable method.

It remains to compute :math:`\xi_{\rm L}` from the integral :eq:`xi_limber`.
Using the relation :eq:`ptoxi` to express the matter correlation function
:math:`\xi` as an integral over the matter power spectrum and exchanging the
order of the integrals yields a representation of the Bessel function
:math:`J_0`,

.. math::

    \frac{1}{\pi} \int_{-\infty}^{\infty} \!
                        \frac{\sin(\sqrt{x^2 + y^2})}{\sqrt{x^2 + y^2}} \, dy
    = J_0(x) \;.

We thus obtain an expression for Limber's matter correlation function in terms
of the matter power spectrum,

.. math::
   :label: ptoxi_limber

    \xi_{\rm L}(r; t)
    = \frac{1}{2\pi} \int_{0}^{\infty} \! P(k; t) \,
                                            \frac{J_0(kr)}{kr} \, k^2 \, dk \;.

This is of similar form as the unequal-time matter correlation function
:eq:`ptoxi_exact`, and we can express both in the generic form

.. math::
   :label: ptoxi_generic

    \xi_{\mu}(r; \ldots)
    = \frac{1}{(2\pi)^{1+\mu}} \int_{0}^{\infty} \! P(k; \ldots) \, 
                                \frac{J_\mu(kr)}{(kr)^{1-\mu}} \, k^2 \, dk \;,

setting :math:`\mu = 0` in the Limber case, and :math:`\mu = 1/2` in the exact
case.  In practice, this allows us to use a single generic implementation of the
FFTLog algorithm to compute either the unequal-time matter correlation function
:eq:`ptoxi_exact` or Limber's matter correlation function :eq:`ptoxi_limber`.


Angular Power Spectrum
----------------------

The angular correlation function :math:`w(\theta)` of a scalar field is related
to its angular power spectrum :math:`C_l` by the sum :eq:`cltow`. The inverse
relation is the integral :cite:`2019arXiv190409973T`

.. math::
   :label: w_to_cl

   C_l = 2\pi \int_{0}^{\pi} \! w(\theta) \,
                                    P_l(\cos\theta) \sin(\theta) \, d\theta \;,

which is yet another difficult oscillatory integral to compute. Clearly, an
alternative approach is needed.

Given a set of angles :math:`\theta_1, \theta_2, \ldots`, the computed angular
correlation function forms the vector :math:`w = (w_k)` with components
:math:`w_k = w(\theta_k)`.  Let :math:`M = (m_{kl})` be the matrix with entries
:math:`m_{kl} = (2l + 1)/(4\pi) \, P_l(\cos\theta_k)` up to some maximum number
:math:`l_{\max}`.  The truncated sum :eq:`cltow` can hence be written

.. math::
   :label: cl_to_w

    w_k
    = \sum_{l=0}^{l_{\max}} m_{kl} \, C_l

or, in matrix form, :math:`w = Mc`, if :math:`c = (C_l)` is the vector of
angular power spectrum entries.

Hence, to obtain :math:`C_l` for :math:`l \le l_{\max}` from :math:`w(\theta)`,
compute sufficiently many values :math:`w_k`, and use a least squares solution
of the matrix equation.


References
----------

.. bibliography:: theory.bib
   :style: plain
