#  © 2025 Valentina Cesare
#  Licensed under the MIT License.
#  See LICENSE file for full license text.

# Author: Valentina Cesare.
# E-Mail address: valentina.cesare@inaf.it
# Creation date: July 31th 2025.

# 1D PLOTS OF NEWTONIAN AND QUMOND NUMERICAL POTENTIALS ALONG R, z AND PHI DIMENSIONS OF THE 3D COMPUTATIONAL GRID.

import numpy as np
import matplotlib.pyplot as plt

# Lettura dei file
# Coordinates
RHalf = np.loadtxt('R_Half.txt')
z = np.loadtxt('z.txt')
phi = np.loadtxt('phi.txt')

phiNCompR = np.loadtxt('phiNComp_MNPlusSpiral_R_zeq0_phieq0_Smooth_Cutoff.txt')
phiQUMONDCompR = np.loadtxt('phiQUMONDComp_MNPlusSpiral_R_zeq0_phieq0_Smooth_Cutoff.txt')

phiNCompz = np.loadtxt('phiNComp_MNPlusSpiral_z_Rsimeq0_phieq0_Smooth_Cutoff.txt')
phiQUMONDCompz = np.loadtxt('phiQUMONDComp_MNPlusSpiral_z_Rsimeq0_phieq0_Smooth_Cutoff.txt')

phiNCompphi = np.loadtxt('phiNComp_MNPlusSpiral_phi_Rsimeq0_zeq0_Smooth_Cutoff.txt')
phiQUMONDCompphi = np.loadtxt('phiQUMONDComp_MNPlusSpiral_phi_Rsimeq0_zeq0_Smooth_Cutoff.txt')

# Creazione dei plot
plt.figure(figsize=(10, 8))

# R vs phiNComp and phiQUMONDComp together
plt.subplot(3, 1, 1)
plt.plot(RHalf, phiNCompR, 'b.-', label=r'$\phi_{N,Comp}$')
plt.plot(RHalf, phiQUMONDCompR, 'r.-', label=r'$\phi_{QUMOND,Comp}$')
plt.title(r'$R$ vs $\phi_{N}$ and $\phi_{QUMOND}$')
plt.xlabel('$R$ $(kpc)$')
plt.ylabel(r'$\phi(R,z=0,\varphi=0) (km/s)^2$')
plt.legend()
plt.grid(True)

# z vs phiNComp and phiQUMONDComp together
plt.subplot(3, 1, 2)
plt.plot(z, phiNCompz, 'b.-', label=r'$\phi_{N,Comp}$')
plt.plot(z, phiQUMONDCompz, 'r.-', label=r'$\phi_{QUMOND,Comp}$')
plt.title(r'$z$ vs $\phi_{N}$ and $\phi_{QUMOND}$')
plt.xlabel('$z$ $(kpc)$')
plt.ylabel(r'$\phi(R \sim 0,z,\varphi=0) (km/s)^2$')
plt.legend()
plt.grid(True)

# \varphi vs phiNComp and phiQUMONDComp together
plt.subplot(3, 1, 3)
plt.plot(phi, phiNCompphi, 'b.-', label=r'$\phi_{N,Comp}$')
plt.plot(phi, phiQUMONDCompphi, 'r.-', label=r'$\phi_{QUMOND,Comp}$')
plt.title(r'$\varphi$ vs $\phi_{N}$ and $\phi_{QUMOND}$')
plt.xlabel(r'$\varphi \, (rad)$')
plt.ylabel(r'$\phi(R \sim 0,z=0,\varphi) \, (km/s)^2$')
plt.legend()
plt.grid(True)

# Layout e salvataggio in PDF
plt.tight_layout()
plt.savefig('Plots_1D_phiNComp_vs_phiQUMONDComp.pdf')
plt.show()

