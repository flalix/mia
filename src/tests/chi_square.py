#-*- coding: utf-8 -*-
'''
Created on Aug 31, 2015

@author: Flavio Lichtenstein
@local: Unifesp DIS - Bioinformatica
'''

import numpy as np
from scipy.stats import chi2
import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 1)
# Calculate a few first moments:

df = 2
mean, var, skew, kurt = chi2.stats(df, moments='mvsk')
# Display the probability density function (pdf):

x = np.linspace(chi2.ppf(0.01, df), chi2.ppf(0.99, df), 100)
ax.plot(x, chi2.pdf(x, df), 'r-', lw=5, alpha=0.6, label='chi2 pdf')

# Alternatively, the distribution object can be called (as a function) to fix the shape, location and scale parameters. This returns a “frozen” RV object holding the given parameters fixed.
# Freeze the distribution and display the frozen pdf:

rv = chi2(df)
ax.plot(x, rv.pdf(x), 'k-', lw=2, label='frozen pdf')

vals = chi2.ppf([0.001, 0.5, 0.999], df)
print np.allclose([0.001, 0.5, 0.999], chi2.cdf(vals, df))

# Generate random numbers:

r = chi2.rvs(df, size=10000)

ax.hist(r, normed=True, histtype='stepfilled', alpha=0.2)
ax.legend(loc='best', frameon=False)
plt.show()

