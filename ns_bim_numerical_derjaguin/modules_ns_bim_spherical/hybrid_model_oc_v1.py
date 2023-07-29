import pandas as pd
import numpy as np

import scipy
from scipy import optimize, integrate
from scipy import interpolate

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib import container
import matplotlib.ticker as ticker
from matplotlib.ticker import NullFormatter
from matplotlib.ticker import LinearLocator, FormatStrFormatter

import time
import math
import random
import multiprocessing as mp
import os
import subprocess
from shutil import copyfile


################################################################################
class constants:
    """Store physical constants for function access."""
    def __init__(self):
        self.e0   = 1.602e-19        # [C]
        self.Na   = 6.022e23         # [mol-1]
        self.kT   = 4.11e-21         # [J]
        self.eps0 = 8.854e-12        # [C V-1 m-1]
        self.h    = 6.626e-34        # [J s]
        self.F    = 96485            # [C mol-1]
        self.kg_per_Da = 1.6605e-27  # [kg Da-1]
        return

class bulk_solution:
    """Store bulk solution data."""
    def __init__(self, ion_str, eps=80.1, m_z=1):
        self.ion_str = ion_str  # [M]
        self.eps     = eps      # [-]
        self.m_z     = m_z      # (Symmetrical) electorlyte valence
        self.get_kappa()
        self.get_ion_dens()
        return

    def get_kappa(self):
        import numpy
        """Convert ionic strength [M] into the inverse Debye length [m-1]."""
        self.kappa = np.sqrt(2*constants().e0**2*self.ion_str*constants().Na*1.0e3/(constants().kT*self.eps*constants().eps0))
        return

    def get_ion_dens(self):
        """Compute number density of ions [m-3] from ionic strength [M]"""
        self.ion_dens = self.ion_str * constants().Na * 1.0e3
        return

    def update_ion_str(self, ion_str):
        self.ion_str = ion_str  # [M]
        self.get_kappa()
        self.get_ion_dens()
        return

################################################################################
# Misc. function(s)

def inverse_Grahame(sigma, bulk_solution):
    """Compute surface potential from charge density for an isolated plane."""
    eps   = bulk_solution.eps
    eps0  = constants().eps0
    kT    = constants().kT
    e0    = constants().e0
    Na    = constants().Na

    c_bulk = bulk_solution.ion_str * Na * 1.0e3
    psi    = 2*kT/e0*np.arcsinh(sigma/np.sqrt(8*eps*eps0*c_bulk*kT))
    return psi

def inverse_capacitor(sigma, bulk_solution):
    eps   = bulk_solution.eps
    eps0  = constants().eps0
    kT    = constants().kT
    e0    = constants().e0
    Na    = constants().Na

    c_bulk = bulk_solution.ion_str * Na * 1.0e3
    psi    = 2*sigma*kT/(e0*np.sqrt(8*eps*eps0*c_bulk*kT))
    return psi

def print_time(time):
    print("{:.2f}".format(time))
    return



################################################################################
# Plotting

def plot_exp_data(axs, file, label):
    """Plot experimental Keq [m] vs. IS [M] data. Use after plotting prediction."""
    df_exp = pd.read_csv(file)

    # if 'Keq_std' in df_exp.columns:
    #     axs.errorbar(df_exp['IS(M)'], df_exp['Keq'], yerr=df_exp['Keq_std'], fmt='o',
    #                  color=plt.gca().lines[-1].get_color(), markersize=10, label=label, capsize=6)
    # else:
    #     axs.scatter(df_exp['IS(M)'], df_exp['Keq'], label=label, marker='o', s=100,
    #                color=plt.gca().lines[-1].get_color())

    axs.scatter(df_exp['IS(M)'], df_exp['Keq'], label=label, marker='o', s=90)

    handles, labels = axs.get_legend_handles_labels()
    new_handles = []
    for h in handles:
        #only need to edit the errorbar legend entries
        if isinstance(h, container.ErrorbarContainer):
            new_handles.append(h[0])
        else:
            new_handles.append(h)
#     axs.legend(new_handles, labels, loc='lower left', bbox_to_anchor=(-0.05, -0.03), frameon=False, handletextpad=0.1)
    axs.legend(new_handles, labels, loc='upper right', frameon=False, handletextpad=0.1)
    #
    # ax.set_xlabel('Ionic strength [M]')
    # ax.set_ylabel(r'$K_{eq}$ [m]')
    return

def plot_model_Keq(ax, file):
    df = pd.read_csv(file, index_col=0)
    ax.plot(df['ion_str_M'], df['Keq_m'], 's-')
    return

def plot_exp_data_v2(axs, df_exp, label):
    """Plot experimental Keq [m] vs. IS [M] data. Use after plotting prediction."""
    axs.scatter(df_exp['IS(M)'], df_exp['Keq'], label=label, marker='o', s=90)
    handles, labels = axs.get_legend_handles_labels()
    new_handles = []
    for h in handles:
        # I only need to edit the errorbar legend entries.
        if isinstance(h, container.ErrorbarContainer):
            new_handles.append(h[0])
        else:
            new_handles.append(h)
    # axs.legend(new_handles, labels, loc='lower left', bbox_to_anchor=(-0.05, -0.03), frameon=False, handletextpad=0.1)
    axs.legend(new_handles, labels, loc='best', frameon=False, handletextpad=0.1)
    return





################################################################################
