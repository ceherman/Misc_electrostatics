import pandas as pd
import numpy as np

import scipy
from scipy import optimize, integrate
from scipy import interpolate

import random
import multiprocessing as mp
import os
import subprocess

from modules_ns_bim_spherical.geometry import *


# Functions for random/MC assignment of charge location

def get_rand_theta_phi():
    theta = 2.0*np.pi*random.random()
    phi   = np.arccos(1.0 - 2.0*random.random())
    return (theta, phi)

def get_x_y_z(radius):
    (theta, phi) = get_rand_theta_phi()
    x = radius*np.sin(phi)*np.cos(theta)
    y = radius*np.sin(phi)*np.sin(theta)
    z = radius*np.cos(phi)
    return (x, y, z)

def place_charges(df_charges, radius_q, center=[0, 0, 0]):
    for i, cont in df_charges.iterrows():
        (x, y, z) = get_x_y_z(radius_q)
        df_charges.at[i, 'x_new_nm'] = x + center[0]
        df_charges.at[i, 'y_new_nm'] = y + center[1]
        df_charges.at[i, 'z_new_nm'] = z + center[2]
    return df_charges

# Functions for using MC results in numerical Derjaguin model

def get_ion_str_from_mc_file_name(name, out_type='str'):
    ind_1 = name.rfind('_') + 1
    ind_2 = name.rfind('.')
    if out_type == 'str':
        return name[ind_1:ind_2]
    elif out_type == 'float':
        return int(name[ind_1:ind_2])/1000
    else:
        print('Invalid output type.')
        return

def get_mc_result_files(mc_res_dir):
    res_files = os.listdir(mc_res_dir)
    res_files.sort()
    pot_files = [i for i in res_files if i[:12] == 'mc_potential']
    chr_files = [i for i in res_files if i[:9]  == 'mc_charge']
    return (pot_files, chr_files)

def find_element_attributes(df_node_links, node_attributes):
    elem_attribs = []
    for i, cont in df_node_links.iterrows():
        c1 = node_attributes[df_node_links.at[i, 'node_1_id'] - 1]
        c2 = node_attributes[df_node_links.at[i, 'node_2_id'] - 1]
        c3 = node_attributes[df_node_links.at[i, 'node_3_id'] - 1]
        elem_attribs.append(np.array([c1, c2, c3]).mean())
    return elem_attribs

def set_element_attributes(args):
    (ion_str_string, i, df_node_links, df_pot, df_chr) = args

    node_pot = df_pot[i]
    elem_pot = find_element_attributes(df_node_links, node_pot)
    node_chr = df_chr[i]
    elem_chr = find_element_attributes(df_node_links, node_chr)

    df_elements                = df_node_links.copy()
    df_elements['potential_V'] = elem_pot
    df_elements['charge_C_m2'] = elem_chr
    return (ion_str_string, i, df_elements)


def mc_get_elem_dict(prot_dir):
    nsbim_program_dir      = prot_dir + 'DH_NS_BIM/'
    mc_res_dir             = prot_dir + 'mc_results/'
    (pot_files, chr_files) = get_mc_result_files(mc_res_dir)

    node_links_file   = nsbim_program_dir + '/rslt_node_linkages.csv'
    node_coords_file  = nsbim_program_dir + '/rslt_node_coordinates.csv'
    df_node_links     = get_node_links(node_links_file)
    df_node_coords    = get_node_coords(node_coords_file)
    get_elements_node_coords_and_centroids(df_node_links, df_node_coords)
    get_elements_outward_unit_normal(df_node_links)

    df            = pd.read_csv(mc_res_dir + pot_files[0])
    iterations    = df.columns
    elements_dict = {}

    args_list = []
    for (pot_file, chr_file) in zip(pot_files, chr_files):
        df_pot         = pd.read_csv(mc_res_dir + pot_file)
        df_chr         = pd.read_csv(mc_res_dir + chr_file)
        ion_str_string = get_ion_str_from_mc_file_name(pot_file)
        for i in iterations:
            args_list.append((ion_str_string, i, df_node_links, df_pot, df_chr))

    pool      = mp.Pool(processes = mp.cpu_count())
    results   = pool.map(set_element_attributes, args_list)
    pool.close()
    pool.join()

    elements_dict = {}
    for res in results:
        (ion_str_string, i, df_elements) = res
        key                              = ion_str_string + '_' + i
        elements_dict[key]               = df_elements
    return elements_dict

def get_ion_str_from_key(key):
    ind = key.rfind('_')
    return int(key[:ind])/1000

def get_iteration_num_from_key(key):
    ind = key.rfind('_')
    return int(key[ind+1:])

################################################################################
# MC for monomer spheres





























################################################################################
