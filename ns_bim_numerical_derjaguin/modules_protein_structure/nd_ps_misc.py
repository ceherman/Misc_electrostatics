import numpy as np
import pandas as pd
from matplotlib import container

################################################################################
# Classes

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
        """Convert ionic strength [M] into the inverse Debye length [m-1]."""
        self.kappa = np.sqrt(2*constants().e0**2*self.ion_str*constants().Na\
        *1.0e3/(constants().kT*self.eps*constants().eps0))
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

class protein_class:
    """Store protein data."""
    def __init__(self, name):
        self.name = name

################################################################################
# Misc. function(s)

def print_time(time):
    print("{:.8f}".format(time))
    return

def read_str_list(input_file_name):
    input_source_file = open(input_file_name, "r")
    current_str_list  = input_source_file.readlines()
    input_source_file.close()
    return current_str_list

def write_str_list(input_file_name, str_list):
    new_file     = open(input_file_name, 'w')
    new_file.writelines(str_list)
    new_file.close()
    return

def get_ion_str_from_protein_name(protein):
    ind_slice = protein.rfind('_') + 1
    return int(protein[ind_slice:])/1000

################################################################################

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
    axs.legend(new_handles, labels, loc='upper right', frameon=False, handletextpad=0.1)
    return


################################################################################


































################################################################################











































################################################################################
