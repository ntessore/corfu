import numpy as np
import matplotlib.pyplot as plt
from scipy.special import hyp2f1, gamma

k = np.logspace(-2, 2, 1000)

a, b = 2.0, 0.5

p = k**a*np.exp(-b*k)

import corfu

r, xi_ptoxi = corfu.ptoxi(k, p, q=1.0, limber=False)

xi_truth = (r**2+b**2)**(-(a+2)/2)*gamma(a+2)*np.sin((a+2)*np.arctan(r/b))/(2*np.pi**2*r)

r_L, xi_L_ptoxi = corfu.ptoxi(k, p, q=1.0, limber=True)

xi_L_truth = b**(-a-2)*gamma(a+2)*hyp2f1((a+2)/2, (a+3)/2, 1, -(r_L/b)**2)/(2*np.pi*r_L)

fig, ax = plt.subplots(1, 3, figsize=(8.0, 3.2), sharey=True)

ax[0].loglog(k, p, '-k')
ax[0].set_xlabel(r'$k$')
ax[0].set_title(r'$P(k)$')

ax[1].loglog(r, +xi_ptoxi, '-k', label='ptoxi')
ax[1].loglog(r, -xi_ptoxi, '--k')
ax[1].loglog(r, +xi_truth, '-r', label='truth', zorder=-1)
ax[1].loglog(r, -xi_truth, '--r', zorder=-1)
ax[1].set_xlabel(r'$r$')
ax[1].set_title(r'$\xi(r)$')
ax[1].legend()

ax[2].loglog(r_L, +xi_L_ptoxi, '-k', label='ptoxi')
ax[2].loglog(r_L, -xi_L_ptoxi, '--k')
ax[2].loglog(r_L, +xi_L_truth, '-r', label='truth', zorder=-1)
ax[2].loglog(r_L, -xi_L_truth, '--r', zorder=-1)
ax[2].set_xlabel(r'$r$')
ax[2].set_title(r'$\xi_{\rm L}(r)$')
ax[2].legend()

plt.tight_layout()

plt.show()
