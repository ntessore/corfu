'''corfu: projected angular correlation functions

author: Nicolas Tessore <n.tessore@ucl.ac.uk>
license: MIT

'''

__version__     = '0.2.dev'

__all__ = [
    'ptoxi',
    'theta',
    'xitow',
    'xitow_limber',
    'wtocl',
    'cltow',
]

import numpy as np
from scipy.special import loggamma, poch
from scipy.interpolate import splrep, splev, RectBivariateSpline

# constants
PI = 3.1415926535897932384626433832795028841971693993751E+00
PI_HALF = 1.5707963267948966192313216916397514420985846996876E+00
TWO_PI = 6.2831853071795864769252867665590057683943387987502E+00
FOUR_PI = 1.2566370614359172953850573533118011536788677597500E+01
LN_2 = 6.9314718055994530941723212145817656807550013436026E-01
LN_10 = 2.3025850929940456840179914546843642076011014886288E+00


def ptoxi(k, p, q=0.0, d=0.0, limber=False):
    '''compute 3d correlation function from power spectrum

    Parameters
    ----------
    k : array_like (N,)
        Wavenumbers.
    p : array_like (..., N)
        Power spectrum. Can be multidimensional. Last axis must agree with the
        wavenumber axis.
    q : float, optional
        Exponent of power law bias for fast Hankel transform.
    d : float, optional
        Logarithmic shift of output sequence.
    limber : bool, optional
        Compute Limber correlation function for equal time approximation.
        Default is `False`.

    '''

    assert np.ndim(k) == 1, 'k must be 1d array'
    assert np.shape(p)[-1] == len(k), 'last axis of p must agree with size of k'

    # Limber correlation function or exact
    if limber:
        mu = 0
    else:
        mu = 0.5

    # set up log space k
    n = len(k)
    lnk1 = np.log(k[0])
    lnkn = np.log(k[-1])
    lnkc = (lnk1 + lnkn)/2
    dlnk = (lnkn - lnk1)/(n-1)
    jc = (n-1)/2
    j = np.arange(n)

    # make sure given k is linear in log space
    if not np.allclose(k, np.exp(lnkc + (j-jc)*dlnk)):
        raise ValueError('k array not linear in log space')

    # tweak d to fulfil low-ringing condition
    xp = (mu+1+q)/2
    xm = (mu+1-q)/2
    y = PI_HALF/dlnk
    zp = loggamma(xp + 1j*y)
    zm = loggamma(xm + 1j*y)
    u = (LN_2 - d)/dlnk + (zp.imag + zm.imag)/PI
    d = d + (u - np.round(u))*dlnk

    # compute Hankel transform coefficients
    y = np.linspace(0, np.pi*(n//2)/(n*dlnk), n//2+1)
    u = np.empty(n//2+1, dtype=complex)
    v = np.empty(n//2+1, dtype=complex)
    u.imag[:] = y
    u.real[:] = xm
    loggamma(u, out=v)
    u.real[:] = xp
    loggamma(u, out=u)
    y *= 2*(LN_2 - d)
    u.real -= v.real
    u.real += LN_2*q
    u.imag += v.imag
    u.imag += y
    np.exp(u, out=u)

    # fix last coefficient to be real
    u.imag[-1] = 0

    # deal with special cases
    if not np.isfinite(u[0]):
        # write u_0 = 2^q Gamma(xp)/Gamma(xm) = 2^q poch(xm, xp-xm)
        # poch() handles special cases for negative integers correctly
        u[0] = 2**q * poch(xm, xp-xm)
        # the transform may still be singular
        if np.isinf(u[0]):
            raise ValueError(f'singular transform for q = {q}')

    # ensure that kr is good for n even
    if not n&1 and not np.isclose(u[-1].imag, 0):
        raise ValueError('unable to construct low-ringing transform, '
                         'try odd number of points or different q')

    # input array for transform
    xi = np.copy(p)  # allocates memory

    # factor of (k/k_c)^{mu+1-q} for input array
    xi *= np.exp((mu+1-q)*(j-jc)*dlnk)

    # Hankel transform via real FFT
    xi = np.fft.rfft(xi, axis=-1)
    xi *= u
    xi = np.fft.irfft(xi, n, axis=-1)
    xi[..., :] = xi[..., ::-1]

    # factor of (r/r_c)^{mu+1-q} (k_c r_c)^{mu+1-q} for output array
    xi *= np.exp((mu+1-q)*((j-jc)*dlnk + d))

    # set up r in log space
    r = np.exp(d)/k[::-1]

    # prefactor for correlation function
    xi /= TWO_PI**(1+mu)
    xi /= r**3

    # done, return separations and correlations
    return r, xi


def theta(n):
    r'''compute points and weights for the angular correlation function

    Returns :math:`n` angles :math:`\theta_1, \ldots, \theta_n` and associated
    weights at which the angular correlation function should be evaluated for
    an optimal estimate of the angular power spectrum using :func:`corfu.wtocl`.

    Parameters
    ----------
    n : int
        Number of points to return.

    Returns
    -------
    theta : array_like (n,)
        Angles in radians.
    weights : array_like (n,)
        Weights for optimal angular power spectrum estimation.

    Warnings
    --------
    This function currently uses :func:`numpy.polynomial.legendre.leggauss`. In
    the future, a better implementation for higher-order modes will be required.

    Notes
    -----
    The points and weights correspond to a :math:`n`-point Gauss-Legendre
    quadrature rule, for which the angular power spectrum estimate
    :math:`\hat{C}_l` only contains errors from modes larger than :math:`2n-l`.

    '''

    x, w = np.polynomial.legendre.leggauss(n)

    return np.arccos(x[::-1]), w[::-1]


def xitow(theta, f1, f2, xi, progress=False):
    '''unequal time projection'''

    assert np.ndim(theta) == 1, 'theta must be 1d array'
    assert len(f1) == 2, 'f1 must be tuple of radii and weights'
    assert np.ndim(f1[0]) == 1, 'f1[0] must be 1d array'
    assert np.ndim(f1[1]) == 1, 'f1[1] must be 1d array'
    assert np.shape(f1[0]) == np.shape(f1[1]), 'shapes of f1[0] and f1[1] must match'
    assert len(f2) == 2, 'f2 must be tuple of radii and weights'
    assert np.ndim(f2[0]) == 1, 'f2[0] must be 1d array'
    assert np.ndim(f2[1]) == 1, 'f2[1] must be 1d array'
    assert np.shape(f2[0]) == np.shape(f2[1]), 'shapes of f2[0] and f2[1] must match'
    assert len(xi) == 4, 'xi must be tuple of radii, radii, separations, and correlations'
    assert np.ndim(xi[0]) == 1, 'xi[0] must be 1d array'
    assert np.ndim(xi[1]) == 1, 'xi[1] must be 1d array'
    assert np.ndim(xi[2]) == 1, 'xi[2] must be 1d array'
    assert np.ndim(xi[3]) == 3, 'xi[3] must be 3d array'
    assert np.shape(xi[3]) == (len(xi[0]), len(xi[1]), len(xi[2])), 'shape of xi[3] must match xi[0], xi[1], xi[2]'

    # use tqdm to report on progress
    if progress:
        if progress == 'gui':
            from tqdm import tqdm_gui as prog
        elif progress == 'notebook':
            from tqdm.notebook import tqdm as prog
        else:
            from tqdm import tqdm as prog
    else:
        def prog(x, total): return x

    # expand inputs
    x1_f1, f1 = f1
    x2_f2, f2 = f2
    x1_xi, x2_xi, r_xi, xi = xi

    # index array for x1_xi
    n1_xi = np.arange(len(x1_xi))

    # log(separation) values for log-linear interpolation of xi
    lnr_xi = np.log(r_xi)

    # minimum and maximum of log(separation)
    min_lnr_xi = np.min(lnr_xi)
    max_lnr_xi = np.max(lnr_xi)

    # select support = radii where the filters are nonzero
    supp_f1 = np.nonzero(f1)
    supp_f2 = np.nonzero(f2)

    # precompute trig values
    sin_theta = np.sin(theta)
    cos_theta = np.cos(theta)

    # differences of xi along x1 axis
    dxi = np.diff(xi, axis=0)

    # correlation function at fixed x1
    xi1 = np.empty(np.shape(xi)[1:])

    # number of x2 solutions for each (theta, r, x1)
    n_x2 = np.empty(len(theta), dtype=int)
    c_x2 = np.zeros(len(theta)+1, dtype=int)

    # partial results for integration along line of sight x1
    w1 = np.zeros((len(theta), len(f1)))

    # treat each x1 value independently
    for i, x1 in prog(zip(supp_f1[0], x1_f1[supp_f1]), total=len(supp_f1[0])):

        # linearly interpolate xi along x1 axis
        n, u = divmod(np.interp(x1, x1_xi, n1_xi), 1)
        xi1[:, :] = dxi[int(n)]
        xi1 *= u
        xi1 += xi[int(n)]

        # construct interpolator for xi(x2, r | x1)
        xi1_f = RectBivariateSpline(x2_xi, lnr_xi, xi1, kx=1, ky=1)

        # precompute some numbers
        x1_sin_theta = x1*sin_theta
        x1_cos_theta = x1*cos_theta

        # all x2 solutions for this x1 and all theta, r
        x2 = []

        # go through theta, collect solutions for x2 given theta, r, x1
        # combine with the values x2_f2 where the filter is given
        for j, (x1_st, x1_ct) in enumerate(zip(x1_sin_theta, x1_cos_theta)):
            dx = r_xi[r_xi >= x1_st]**2
            dx -= x1_st**2
            np.sqrt(dx, out=dx)
            x2_1 = dx  # same memory
            x2_2 = -dx[dx <= x1_ct]
            x2_1 += x1_ct
            x2_2 += x1_ct
            x2.extend([x2_f2[supp_f2], x2_1, x2_2])
            n_x2[j] = len(x2[-3]) + len(x2[-2]) + len(x2[-1])

        # stack all x2 value subarrays
        x2 = np.concatenate(x2)

        # cumulative counts of x2 values in subarrays, for indexing
        np.cumsum(n_x2, out=c_x2[1:])

        # sort x2 subarrays
        for k, l in zip(c_x2, c_x2[1:]):
            x2[k:l].sort(kind='mergesort')

        # get function values at all distict x2 values
        f2_x2 = np.interp(x2, x2_f2, f2, left=0, right=0)

        # find support of f2_x2
        supp_f2_x2 = np.nonzero(f2_x2)

        # x2 values in support
        x2_s = x2[supp_f2_x2]

        # expand arrays for subarrays
        x1_sin_theta = np.repeat(x1_sin_theta, n_x2)[supp_f2_x2]
        x1_cos_theta = np.repeat(x1_cos_theta, n_x2)[supp_f2_x2]

        # separation array
        lnr = x1_sin_theta**2 + (x1_cos_theta - x2_s)**2
        np.log(lnr, out=lnr)
        lnr *= 0.5

        # bounds check
        assert np.min(lnr) >= min_lnr_xi, 'r < min(r_xi)'
        assert np.max(lnr) <= max_lnr_xi, 'r > max(r_xi)'

        # interpolate correlation values over support
        w2_s = xi1_f.ev(x2_s, lnr)

        # multiply by filter over x2 where nonzero
        w2_s *= f2_x2[supp_f2_x2]

        # full correlation function along line of sight x2
        w2 = np.zeros(len(x2))
        w2[supp_f2_x2] = w2_s

        # integrate along x2 line of sight for this (theta, x1)
        for j, (k, l) in enumerate(zip(c_x2, c_x2[1:])):
            w1[j, i] = np.trapz(w2[k:l], x2[k:l])

    # multiply by filter over x1
    w1 *= f1

    # integrate along x1 line of sight
    w = np.trapz(w1, x1_f1)

    # done, return projection
    return w


def xitow_limber(theta, f12, xi):
    '''equal time projection'''

    assert np.ndim(theta) == 1, 'theta must be 1d array'
    assert len(f12) == 2, 'f12 must be tuple of radii and weights'
    assert np.ndim(f12[0]) == 1, 'f12[0] must be 1d array'
    assert np.ndim(f12[1]) == 1, 'f12[1] must be 1d array'
    assert np.shape(f12[0]) == np.shape(f12[1]), 'shapes of f12[0] and f12[1] must match'
    assert len(xi) == 3, 'xi must be tuple of radii, separations, and correlations'
    assert np.ndim(xi[0]) == 1, 'xi[0] must be 1d array'
    assert np.ndim(xi[1]) == 1, 'xi[1] must be 1d array'
    assert np.ndim(xi[2]) == 2, 'xi[2] must be 2d array'
    assert np.shape(xi[2]) == (len(xi[0]), len(xi[1])), 'shape of xi[2] must match xi[0], xi[1]'

    # expand inputs
    xf, f12 = f12
    xxi, rxi, xi = xi

    # index array for xxi
    nxi = np.arange(len(xxi))

    # log(separation) values for log-linear interpolation of xi
    lnrxi = np.log(rxi)

    # select support = radii where the filter is nonzero
    supp = np.nonzero(f12)

    # partial results for integration along line of sight
    w12 = np.zeros((len(theta), len(xf)))

    # correlation functions at fixed x
    xi12 = np.empty(np.shape(xi)[1:])

    # loop over x in support
    for i, x in zip(supp[0], xf[supp]):

        # linearly interpolate xi along x axis
        n, u = divmod(np.interp(x, xxi, nxi), 1)
        np.multiply(xi[int(n)], 1 - u, out=xi12)
        xi12 += xi[int(n)+1]*u

        # separation array
        lnr = x*theta
        np.log(lnr, out=lnr)

        # bounds check
        minr = np.exp(np.min(lnr))
        maxr = np.exp(np.max(lnr))
        assert minr >= np.min(rxi), f'minimum separation {minr} not in r for xi'
        assert maxr <= np.max(rxi), f'maximum separation {maxr} not in r for xi'

        # set correlations for this x
        w12[:, i] = np.interp(lnr, lnrxi, xi12)

    # multiply by filter over x and theta*x
    w12 *= f12
    w12 *= xf
    w12 *= theta[:, None]

    # integrate along line of sight
    w = np.trapz(w12, xf, axis=1)

    # done, return projection
    return w


def wtocl(w, theta, weights, lmax=None):
    assert np.ndim(w) >= 1, 'w must be at least 1d array'
    assert np.ndim(theta) == 1, 'theta must be 1d array'
    assert len(theta) == np.shape(w)[-1], 'shapes of theta and w must agree'
    assert np.ndim(weights) == 1, 'weights must be 1d array'
    assert len(weights) == len(theta), 'length of weights and theta must agree'
    if lmax is not None:
        assert np.isscalar(lmax) and lmax > 0, 'lmax must be a positive number'

    # default lmax is n, else convert to integer
    if lmax is None:
        lmax = len(theta)
    else:
        lmax = int(lmax)

    # ell values
    ell = np.arange(lmax+1)

    # compute the Legendre Vandermonde matrix
    p = np.polynomial.legendre.legvander(np.cos(theta), lmax)

    # compute Cls
    c = TWO_PI*np.dot(weights*w, p)

    # done
    return c


def cltow(cl, theta):
    assert np.ndim(cl) == 1, 'cl must be 1d array'
    assert np.ndim(theta) == 1, 'theta must be 1d array'

    # get evaluation points for Legendre series
    x = np.cos(theta)

    # ell numbers
    ell = np.arange(len(cl))

    # series coefficients from Cls
    c = (2*ell+1)/FOUR_PI*cl

    # evaluate Legendre polynomial
    w = np.polynomial.legendre.legval(x, c)

    # done, return correlation function
    return w
