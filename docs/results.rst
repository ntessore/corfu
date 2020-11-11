Results
=======

Correlation function
--------------------

We show that computing the 3d correlation function :math:`\xi(r)` from the 3d
power spectrum :math:`P(k)` works.  This is implemented in the
:func:`corfu.ptoxi` function, using the FFTLog algorithm.

First, we import a number of required packages.  We will need the special
functions :data:`~scipy.special.hyp2f1` and :data:`~scipy.special.gamma`.

.. literalinclude:: results/ptoxi.py
   :lines: 1-3

For testing, we will define a 3d power spectrum :math:`P(k) = k^a \, e^{-bk}`
where :math:`a` and :math:`b` are arbitrary positive parameters.  We set up a
grid of wavenumbers :math:`k` over which we evaluate the power spectrum.

.. literalinclude:: results/ptoxi.py
   :lines: 5

We fix the values of the parameters :math:`a` and :math:`b`.

.. literalinclude:: results/ptoxi.py
   :lines: 7

Finally, we are able to compute the input 3d power spectrum.

.. literalinclude:: results/ptoxi.py
   :lines: 9

The basic setup is complete.  We now use ``corfu`` to convert the 3d power
spectrum to a 3d correlation function.

.. literalinclude:: results/ptoxi.py
   :lines: 11

The computation is done with the :func:`corfu.ptoxi` function.  Here, we compute
:math:`\xi(r)` for the exact (non-Limber) projection.  The FFTLog bias parameter
:math:`q = 1.0` is chosen to improve the result while not causing numerical
issues.

.. literalinclude:: results/ptoxi.py
   :lines: 13

The true 3d correlation function is the integral :eq:`ptoxi-exact`, which we can
evaluate analytically for our sample power spectrum,

.. math::

    \begin{aligned}
    \xi(r)
    & = \frac{1}{(2\pi)^{3/2}} \int_{0}^{\infty} \! k^a \, e^{-bk} \,
                                \frac{J_{1/2}(kr)}{\sqrt{kr}} \, k^2 \, dk \\
    &= \frac{(b^2+r^2)^{-\frac{a+2}{2}} \,
                \Gamma(a+2) \sin\bigl((a+2) \arctan(r/b)\bigr)}{2\pi^2 r} \;.
    \end{aligned}

Using the :math:`r` values that the call to :func:`corfu.ptoxi` returned, the
truth values are readily computed using SciPy's :data:`~scipy.special.gamma`
function:

.. literalinclude:: results/ptoxi.py
   :lines: 15

We repeat the steps for Limber's correlation function :math:`\xi_{\rm L}(r)`,
which is returned by :func:`corfu.ptoxi` when setting ``limber=True``:

.. literalinclude:: results/ptoxi.py
   :lines: 17

The analytical results for Limber's correlation function :eq:`ptoxi-limber` can
similarly be obtained,

.. math::

    \begin{aligned}
    \xi_{\rm L}(r)
    & = \frac{1}{2\pi} \int_{0}^{\infty} \! k^a \, e^{-bk} \,
                                \frac{J_0(kr)}{kr} \, k^2 \, dk \\
    &= \frac{b^{-a-2} \, \Gamma(a+2) \, _2F_1\bigl(\frac{a+2}{2}, \frac{a+3}{2}; 1; -\frac{r^2}{b^2}\bigr)}{2\pi r} \;,
    \end{aligned}

where :math:`_2F_1(a, b; c; z)` is the hypergeometric function, available
through SciPy's :data:`scipy.special.hyp2f1`:

.. literalinclude:: results/ptoxi.py
   :lines: 19

The results of the comparison are shown in :numref:`ptoxi-results`. It is clear
that :func:`corfu.ptoxi` produces the desired output, up to aliasing due to the
circular nature of the underlying FFTLog transform.

.. _ptoxi-results:

.. plot:: results/ptoxi.py

   Computing the 3d correlation function :math:`\xi(r)` and Limber's correlation function :math:`\xi_{\rm L}(r)` from the 3d power spectrum :math:`P(k)` using :func:`corfu.ptoxi`.
