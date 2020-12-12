#
# authors:        L. Pezzini
# e-mail :        luca.pezzini@edu.unito.it
# date:           12.12.2020
# MIT license
#

# Implementation of the 1D acoustic wave equation

import numpy as np
import matplotlib
# Show Plot in The Notebook
matplotlib.use("nbagg")
import matplotlib.pyplot as plt

# Sub-plot Configuration
# ----------------------
from matplotlib import gridspec

# Ignore Warning Messages
# -----------------------
import warnings
warnings.filterwarnings("ignore")

# Parameter Configuration
# -----------------------

nx   = 10000        # number of grid points in x-direction
xmax = 10000        # physical domain (m)
dx   = xmax/(nx-1)  # grid point distance in x-direction
c0   = 334.         # wave speed in medium (m/s)
isrc = int(nx/2)    # source location in grid in x-direction
ir   = isrc + 100          # receiver location in grid in x-direction
nt   = 600          # maximum number of time steps
dt   = 0.002        # time step
op   = 3            # Length of FD operator for 2nd space derivative (3 or 5)
print(op, '- point operator')

# Source time function parameters
f0   = 15. # dominant frequency of the source (Hz)
t0   = 4. / f0 # source time shift

# CPL Stability Criterion
# --------------------------
eps = c0 * dt / dx #epsilon value (corrected May 3, 2020)

print('Stability criterion =', eps)

# Source time function (Gaussian)
# -------------------------------
src  = np.zeros(nt + 1)
time = np.linspace(0 * dt, nt * dt, nt)
# 1st derivative of a Gaussian (corrected May 3, 2020)
src  = -8. * (time - t0) * f0 * (np.exp(-1.0 * (4*f0) ** 2 * (time - t0) ** 2))

# Plot source time function

# Plot position configuration
# ---------------------------
plt.ion()
fig1 = plt.figure(figsize=(10, 6))
gs1  = gridspec.GridSpec(1, 2, width_ratios=[1, 1], hspace=0.3, wspace=0.3)

# Plot source time function
# -------------------------
ax1  = plt.subplot(gs1[0])
ax1.plot(time, src) # plot source time function
ax1.set_title('Source Time Function')
ax1.set_xlim(time[0], time[-1])
ax1.set_xlabel('Time (s)')
ax1.set_ylabel('Amplitude')

# Plot source spectrum
# --------------------
ax2  = plt.subplot(gs1[1])
spec = np.fft.fft(src) # source time function in frequency domain
freq = np.fft.fftfreq(spec.size, d = dt ) # time domain to frequency domain (corrected May 3, 2020)
ax2.plot(np.abs(freq), np.abs(spec)) # plot frequency and amplitude
ax2.set_xlim(0, 250) # only display frequency from 0 to 250 Hz
ax2.set_title('Source Spectrum')
ax2.set_xlabel('Frequency (Hz)')
ax2.set_ylabel('Amplitude')

ax2.yaxis.tick_right()
ax2.yaxis.set_label_position("right")

plt.show()
