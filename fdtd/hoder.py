#
# authors:        L. Pezzini
# e-mail :        luca.pezzini@edu.unito.it
# date:           11.12.2020
# MIT license
#

#
#    Learn how to define high-order central finite-difference operators
#    Investigate the behaviour of the operators with increasing length
#

import math
import numpy as np
import matplotlib.pyplot as plt

# Define function to calculate Taylor operators
def central_difference_coefficients(nop, n):
    """
    Calculate the central finite difference stencil for an arbitrary number
    of points and an arbitrary order derivative.

    :param nop: The number of points for the stencil. Must be
        an odd number.
    :param n: The derivative order. Must be a positive number.
    """
    m = np.zeros((nop, nop))
    for i in range(nop):
        for j in range(nop):
            dx = j - nop // 2
            m[i, j] = dx ** i

    s = np.zeros(nop)
    s[n] = math.factorial(n)

    # The following statement return oper = inv(m) s
    oper = np.linalg.solve(m, s)
    # Calculate operator
    return oper

# Calculate and plot Taylor operator

# Give length of operator (odd)
nop = 25
# Give order of derivative (0 - interpolation, 1 - first derivative, 2 - second derivative)
n = 1

# Get operator from routine 'central_difference_coefficients'
oper = central_difference_coefficients(nop, n)

# Plot operator
x = np.linspace(-(nop - 1) / 2, (nop - 1) / 2, nop)

# Simple plot with operator
plt.figure(figsize=(10, 4))
plt.plot(x, oper,lw=2,color='blue')
plt.plot(x, oper,lw=2,marker='o',color='blue')
plt.plot(0, 0,lw=2,marker='o',color='red')
#plt.plot (x, nder5-ader, label="Difference", lw=2, ls=":")
plt.title("Taylor Operator with nop =  %i " % nop )
plt.xlabel('x')
plt.ylabel('Operator')
plt.grid()
plt.show()
