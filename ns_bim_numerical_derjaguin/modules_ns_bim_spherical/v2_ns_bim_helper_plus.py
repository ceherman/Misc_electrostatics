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

from modules_ns_bim_spherical.hybrid_model_oc_v1 import *
from modules_ns_bim_spherical.numerical_Derjaguin_v1 import *

################################################################################
def my_mk_dir(path):
    try:
        os.makedirs(path)
    except OSError:
        print ("Creation of the directory %s failed" % path)
    else:
        print ("Successfully created the directory %s" % path)
    return

def my_copy(file, destination):
    return copyfile(file, destination)

def get_charge_coords(fort_file, radius_shell):
    df_charges              = pd.read_fwf(fort_file, skiprows=1,
                              names=['x_A', 'y_A', 'z_A', 'q', 'residue', 'res_num',
                                     'chain', 'unk'], infer_nrows=1000)
    df_charges['radii_old'] = np.sqrt(df_charges['x_A']**2 + df_charges['y_A']**2 + \
                                     df_charges['z_A']**2)
    df_charges['x_new_nm']  = df_charges['x_A'] * radius_shell/df_charges['radii_old']
    df_charges['y_new_nm']  = df_charges['y_A'] * radius_shell/df_charges['radii_old']
    df_charges['z_new_nm']  = df_charges['z_A'] * radius_shell/df_charges['radii_old']
    return df_charges

def get_charge_coords_his_pos_charge(fort_file, radius_shell):
    df_charges              = pd.read_fwf(fort_file, skiprows=1,
                              names=['x_A', 'y_A', 'z_A', 'q', 'residue', 'res_num',
                                     'chain', 'unk'], infer_nrows=1000)
    df_charges['radii_old'] = np.sqrt(df_charges['x_A']**2 + df_charges['y_A']**2 + \
                                     df_charges['z_A']**2)
    df_charges['x_new_nm']  = df_charges['x_A'] * radius_shell/df_charges['radii_old']
    df_charges['y_new_nm']  = df_charges['y_A'] * radius_shell/df_charges['radii_old']
    df_charges['z_new_nm']  = df_charges['z_A'] * radius_shell/df_charges['radii_old']

    for row, contents in df_charges.iterrows():
        if contents['residue'] == 'HIS':
            df_charges.at[row, 'q'] = 0.5
    return df_charges

def get_charge_coords_pqr(file):
    df_charges = pd.read_csv(file, skiprows=1, names=['residue', 'res_num',
                             'x_new_nm', 'y_new_nm', 'z_new_nm', 'q', 'num'])
    return df_charges

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

def update_input_source(df_charges, input_file_name):
    new_src_list = []
    for row, contents in df_charges.iterrows():
        charge_id = str(row + 1)
        x = str(contents['x_new_nm'])
        y = str(contents['y_new_nm'])
        z = str(contents['z_new_nm'])
        q = str(contents['q'])
        new_str = charge_id + ',' + x + ',' + y + ',' + z + ',' + q + '\n'
        new_src_list.append(new_str)

    current_str_list = read_str_list(input_file_name)
    new_str_list = current_str_list[:4] + new_src_list + [current_str_list[-1]]
    write_str_list(input_file_name, new_str_list)
    return

def update_input_geom(radius_new, input_file_name):
    current_str_list     = read_str_list(input_file_name)
    r_line               = current_str_list[-2]
    current_str_list[-2] = str(radius_new) + ',' + r_line[(r_line.find(',')+1):]
    write_str_list(input_file_name, current_str_list)
    return

def update_input_phys_dh(df_charges, kappa_inv_nm, input_file_name):
    current_str_list     = read_str_list(input_file_name)
    current_str_list[1]  = str(kappa_inv_nm) + 'd0, ' + current_str_list[1][-7:]
    current_str_list[-2] = str(len(df_charges.index)) + '\n'
    write_str_list(input_file_name, current_str_list)
    return

def update_ns_bim_inputs(bim_dir, df_charges, radius_new, ion_str):
    sol = bulk_solution(ion_str)
    kappa_inv_nm = sol.kappa/1.0e9

    update_input_source(df_charges, bim_dir + 'Input_Source.dat')
    update_input_geom(radius_new, bim_dir + 'Input_Geom.dat')
    update_input_phys_dh(df_charges, kappa_inv_nm, bim_dir + 'Input_Phys_DH.dat')
    return

def run_ns_bim(bim_dir):
    process = subprocess.run(['./run_DH'], cwd=bim_dir,
                         stdout=subprocess.PIPE, universal_newlines=True)
    # print(process.stdout)
    return

def copy_ns_bim_data(bim_dir, destination):
    data_files = ['rslt_node_coordinates.csv', 'rslt_node_linkages.csv',
                  'rslt_protein_potential_variable_radius.csv',
                  'rslt_protein_charge_variable_radius.csv']
    for i in data_files:
        my_copy(bim_dir + i, destination + i)
    return

def meta_run_ns_bim(prot_dir, radius_new, ion_str, destination):
    # Updated to use ion_str
    my_mk_dir(destination)
    offset       = 0.15 # NOTE THIS DOESN'T WORK; HARD-CODED IN NS-BIM # 0.15 # [nm]
    radius_shell = radius_new - offset

    fort_file  = prot_dir + 'fort.3'
    bim_dir    = prot_dir + 'DH_NS_BIM/'

    df_charges = get_charge_coords(fort_file, radius_shell)
    update_ns_bim_inputs(bim_dir, df_charges, radius_new, ion_str)
    run_ns_bim(bim_dir)
    copy_ns_bim_data(bim_dir, destination)
    return

def meta_run_ns_bim_his_pos_charge(prot_dir, radius_new, ion_str, destination):
    # Updated to use ion_str
    my_mk_dir(destination)
    offset       = 0.15 # NOTE THIS DOESN'T WORK; HARD-CODED IN NS-BIM # 0.15 # [nm]
    radius_shell = radius_new - offset

    fort_file  = prot_dir + 'fort.3'
    bim_dir    = prot_dir + 'DH_NS_BIM/'

    df_charges = get_charge_coords_his_pos_charge(fort_file, radius_shell)
    update_ns_bim_inputs(bim_dir, df_charges, radius_new, ion_str)
    run_ns_bim(bim_dir)
    copy_ns_bim_data(bim_dir, destination)
    return

def meta_run_ns_bim_pqr_charge(prot_dir, radius_new, ion_str, destination):
    # Updated to use ion_str
    my_mk_dir(destination)

    file  = prot_dir + 'pqr_charges.csv'
    bim_dir    = prot_dir + 'DH_NS_BIM/'

    df_charges = get_charge_coords_pqr(file)
    update_ns_bim_inputs(bim_dir, df_charges, radius_new, ion_str)
    run_ns_bim(bim_dir)
    copy_ns_bim_data(bim_dir, destination)
    return

def meta_run_ns_bim_pqr_partial_charge(prot_dir, radius_new, ion_str, destination):
    # Updated to use ion_str
    my_mk_dir(destination)

    file  = prot_dir + 'pqr_partial_charges.csv'
    bim_dir    = prot_dir + 'DH_NS_BIM/'

    df_charges = get_charge_coords_pqr(file)
    update_ns_bim_inputs(bim_dir, df_charges, radius_new, ion_str)
    run_ns_bim(bim_dir)
    copy_ns_bim_data(bim_dir, destination)
    return

def meta_run_ns_bim_pqr_charge_surface(prot_dir, radius_new, ion_str, destination):
    # Updated to use ion_str
    my_mk_dir(destination)

    file  = prot_dir + 'pqr_charges_surface.csv'
    bim_dir    = prot_dir + 'DH_NS_BIM/'

    df_charges = get_charge_coords_pqr(file)
    update_ns_bim_inputs(bim_dir, df_charges, radius_new, ion_str)
    run_ns_bim(bim_dir)
    copy_ns_bim_data(bim_dir, destination)
    return

##################################################################
# For use with monomer spheres
def update_input_phys_dh_v2(kappa_inv_nm, input_file_name):
    current_str_list     = read_str_list(input_file_name)
    current_str_list[1]  = str(kappa_inv_nm) + 'd0, ' + current_str_list[1][-7:]
    write_str_list(input_file_name, current_str_list)
    return

def copy_ns_bim_data_v2(bim_dir, destination):
    data_files = ['rslt_nodes_coords_and_pots.csv', 'rslt_node_linkages.csv']
    for i in data_files:
        my_copy(bim_dir + i, destination + i)
    return

def meta_run_ns_bim_v2(prot_dir, ion_str, destination):
    my_mk_dir(destination)
    bim_dir = prot_dir + 'charge_and_potential_base/'

    sol = bulk_solution(ion_str)
    kappa_inv_nm = sol.kappa/1.0e9
    update_input_phys_dh_v2(kappa_inv_nm, bim_dir + 'Input_Phys_DH.dat')

    run_ns_bim(bim_dir)
    copy_ns_bim_data_v2(bim_dir, destination)
    return































################################################################################
# Boneyard

# def copy_ns_bim_data_2(bim_dir, destination, offset):
#     data_files = ['rslt_node_coordinates.csv', 'rslt_node_linkages.csv',
#                   'rslt_protein_potential_variable_radius.csv']
#     for i in data_files:
#         my_copy(bim_dir + i, destination + i)
#     return
#
# def meta_run_ns_bim_2(prot_dir, radius_new, destination, offsets):
#     fort_file  = prot_dir + 'fort.3'
#     bim_dir    = prot_dir + 'DH_NS_BIM/'
#
#     for offset in offsets: # [nm]
#         radius_shell = radius_new - offset
#         df_charges   = get_charge_coords(fort_file, radius_shell)
#         update_ns_bim_inputs(bim_dir, df_charges, radius_new)
#         run_ns_bim(bim_dir)
#         copy_ns_bim_data_2(bim_dir, destination, offset)
#
#     return
