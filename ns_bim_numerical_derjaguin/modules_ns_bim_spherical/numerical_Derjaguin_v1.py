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

from modules_ns_bim_spherical.geometry import *
from modules_ns_bim_spherical.hybrid_model_oc_v1 import *
from modules_ns_bim_spherical.oc_plots_v1 import *
from modules_ns_bim_spherical.orthogonal_collocation_v1 import *
from modules_ns_bim_spherical.v2_ns_bim_helper_plus import *

################################################################################

def map_pot_lin_to_nonlin(lin_pot, bulk_solution):
    sig_lin = constants().eps0*bulk_solution.eps*bulk_solution.kappa*lin_pot
    denom = np.sqrt(8.0*constants().eps0*bulk_solution.eps*bulk_solution.ion_dens*constants().kT)
    return 2.0*constants().kT/constants().e0*np.arcsinh(sig_lin/denom)

def test_a():
    return(constants().e0)

def get_interp_fun(energy_map_file):
    df = pd.read_csv(energy_map_file, index_col=0)
    # Reverse the row order for interpolation
    df = df.reindex(index=df.index[::-1])

    # Get separation (x, in kap_h) and potential (y, in V) data
    x = []
    y = []
    for col in df.columns:
        y.append(float(col))
    for row in df.index:
        x.append(row)

    # Obtain response surface for interpolation
    Z = np.zeros((len(x), len(y)))
    for row in range(len(x)):
        for col in range(len(y)):
            Z[row, col] = df.values[row, col]
    return interpolate.RectBivariateSpline(x, y, Z)

def use_interp_fun(fun, kap_h, pot):
    if kap_h > 9.0 or kap_h < 0.01:  # or pot < -0.345 or pot > 0.3
                                     # (if avoiding potential extrapolation)
        return 0.0
    else:
        return fun(kap_h, pot)[0][0]

def get_integrand(energy):
    # return np.exp(-1.0*energy/constants().kT - 1.0) # Previous, with mistake_1
    return np.exp(-1.0*energy/constants().kT) - 1.0

def get_integrand_list(args):
    (row_plane, elements_df, zs, energy_fun, bulk_solution) = args
    orient_df_1   = get_orientation_df(elements_df, row_plane)

    energy        = np.zeros(len(zs))
    integrand_res = np.zeros(len(zs))

    for i, z in enumerate(zs):
        orient_df_1['z_m']    = (orient_df_1['0_dist_nm'].values)*1.0e-9 + z

        # orient_df_1['f_plus'] = orient_df_1.apply(lambda row:
        #                         use_interp_fun(energy_fun,
        #                                        row['z_m']*bulk_solution.kappa,
        #                                        row['potential_V']), axis=1)
        orient_df_1['f_plus'] = energy_fun.ev(orient_df_1['z_m'].values*bulk_solution.kappa,
                                              orient_df_1['potential_V'].values)

        orient_df_1['f_J']    = unscale_energy(orient_df_1['f_plus'].values,
                                               (orient_df_1['area_nm2'].values)*1.0e-18,
                                               bulk_solution)
        energy[i]        = np.sum(orient_df_1['f_J'])
        integrand_res[i] = get_integrand(energy[i])
    return (row_plane, integrand_res, energy)

def get_integrand_list_cc(args):
    (row_plane, elements_df, zs, energy_fun, bulk_solution) = args
    orient_df_1   = get_orientation_df(elements_df, row_plane)

    energy        = np.zeros(len(zs))
    integrand_res = np.zeros(len(zs))

    for i, z in enumerate(zs):
        orient_df_1['z_m']    = (orient_df_1['0_dist_nm'].values)*1.0e-9 + z
        orient_df_1['f_plus'] = energy_fun.ev(orient_df_1['z_m'].values*bulk_solution.kappa,
                                              orient_df_1['charge_C_m2'].values)

        orient_df_1['f_J']    = unscale_energy(orient_df_1['f_plus'].values,
                                               (orient_df_1['area_nm2'].values)*1.0e-18,
                                               bulk_solution)
        energy[i]        = np.sum(orient_df_1['f_J'])
        integrand_res[i] = get_integrand(energy[i])
    return (row_plane, integrand_res, energy)


def get_integrand_interp_funs(zs, averages):
    index_neg = 0
    while averages[index_neg] > 0.0:
        index_neg += 1

    z_t  = zs[index_neg-1]
    zs_1 = zs[0:index_neg]
    zs_2 = zs[(index_neg-1):]

    a_1_log = np.log(averages[0:index_neg])
    a_2     = averages[(index_neg-1):]

    fun_1_log = interpolate.interp1d(zs_1, a_1_log, kind='cubic', bounds_error=True)
    fun_2     = interpolate.interp1d(zs_2, a_2, kind='cubic', bounds_error=True)
    return(z_t, fun_1_log, fun_2)

def use_integrand_interp(z_value, z_t, fun_1_log, fun_2, last_z, first_z):
    if z_value < first_z:
        return 1.0e20
    elif z_value > last_z:
        return 0.0
    else:
        if z_value < z_t:
            a_log = fun_1_log(z_value)
            a     = np.exp(a_log)
        else:
            a     = fun_2(z_value)
        return a

def get_integrand_interp_funs_wrapper(protein, elements_df, ion_str, resin):
    sol           = bulk_solution(ion_str)
    integrand_dir = 'integrand_results/'
    energy_dir    = 'energy_results/'
    file = resin + '_' + protein + '_' + str(int(sol.ion_str*1.0e3)) + '.csv'

    integrand_df = pd.read_csv(integrand_dir + file, index_col=0)
    zs           = np.zeros(len(integrand_df.columns))
    averages     = np.zeros(len(integrand_df.columns))

    for k in enumerate(integrand_df.items()):
        (index, (z_str, integrands)) = k
        zs[index]       = float(z_str)
        averages[index] = np.average(integrands)

    z_t, fun_1_log, fun_2 = get_integrand_interp_funs(zs, averages)
    return z_t, fun_1_log, fun_2, zs[-1], zs[0]

def eval_Keq_integral(z0, last_z, first_z, z_t, fun_1_log, fun_2):
    Keq_integral_res = integrate.quad(use_integrand_interp, z0, last_z,
                                      args=(z_t, fun_1_log, fun_2, last_z, first_z))
    # assert(Keq_integral_res[1]/Keq_integral_res[0] < 1.0e-3)

    # temp = last_z
    # while (Keq_integral_res[1]/Keq_integral_res[0] > 1.0e-2):
    #     last_z*=0.95
    #     Keq_integral_res = integrate.quad(use_integrand_interp, z0, last_z,
    #                                       args=(z_t, fun_1_log, fun_2, last_z, first_z))
    #     if last_z/temp < 0.05:
    #         break
    return Keq_integral_res[0]

def get_model_Keq(ion_str_list, protein, elements_df, resin, z0):
    Keq_results = np.zeros(len(ion_str_list))
    for j_index, ion_str in enumerate(ion_str_list):
        z_t, fun_1_log, fun_2, last_z, first_z = \
        get_integrand_interp_funs_wrapper(protein, elements_df, ion_str, resin)
        Keq_results[j_index] = eval_Keq_integral(z0, last_z, first_z, z_t, fun_1_log, fun_2)
    return Keq_results

def get_model_Keq_log_spline(ion_str_list, Keq_results):
    ion_str_log = np.log(ion_str_list)
    Keq_log     = np.log(Keq_results)
    return interpolate.interp1d(ion_str_log, Keq_log, kind='cubic', bounds_error=True)

def use_model_Keq_log_spline(fun, ion_str):
    x = np.log(ion_str)
    y = fun(x)
    return np.exp(y)

def get_exp_data(protein, resin):
    if resin == 'sep':
        file = 'experimental_data/' + protein[:3] + '_' + resin + 'h_pH' + \
               protein[-1] + '.csv'
    else:
        file = 'experimental_data/' + protein[:3] + '_' + resin + '_pH' + \
               protein[-1] + '.csv'
    return pd.read_csv(file)

def get_Keq_residual(z0, model_ion_str_list, protein, elements_df, resin):
    model_Keq = get_model_Keq(model_ion_str_list, protein, elements_df, resin, z0)
    spline    = get_model_Keq_log_spline(model_ion_str_list, model_Keq)
    df_exp    = get_exp_data(protein, resin)

    res_vec = np.zeros(len(df_exp['IS(M)']))
    df_exp['model_Keq'] = use_model_Keq_log_spline(spline, df_exp['IS(M)'])
    df_exp['Keq_res']   = np.log(df_exp['Keq']/df_exp['model_Keq'])
    return np.array(df_exp['Keq_res'])

def get_optimal_z0(args):
    (model_ion_str_list, protein, elements_df, resin) = args
    guess = 0.25e-9
    fit = scipy.optimize.least_squares(get_Keq_residual, guess, method='lm',
                                       args=(model_ion_str_list, protein,
                                       elements_df, resin))
    return [protein, resin, fit.x[0]]


################################################################################


def get_integrand_interp_funs_wrapper_v2(file, ion_str, integrand_dir,
                                         integrand_df=pd.DataFrame()):
    if integrand_df.empty:
        integrand_df = pd.read_csv(integrand_dir + file, index_col=0)
    else:
        pass

    zs       = np.zeros(len(integrand_df.columns))
    averages = np.zeros(len(integrand_df.columns))

    for k in enumerate(integrand_df.items()):
        (index, (z_str, integrands)) = k
        zs[index]       = float(z_str)
        averages[index] = np.average(integrands)

    z_t, fun_1_log, fun_2 = get_integrand_interp_funs(zs, averages)
    return z_t, fun_1_log, fun_2, zs[-1], zs[0]

def get_model_Keq_v2(ion_str_list, integrand_df_file, z0, integrand_dir,
                     integrand_df=pd.DataFrame()):
    Keq_results = np.zeros(len(ion_str_list))
    for j_index, ion_str in enumerate(ion_str_list):
        z_t, fun_1_log, fun_2, last_z, first_z = \
        get_integrand_interp_funs_wrapper_v2(integrand_df_file, ion_str, integrand_dir, integrand_df)
        Keq_results[j_index] = eval_Keq_integral(z0, last_z, first_z, z_t, fun_1_log, fun_2)
    return Keq_results

# def get_Keq_residual_v2(z0, model_ion_str_list, protein, resin, integrand_dir):
#     model_Keq = get_model_Keq_v2(ion_str_list, integrand_df_file, z0, integrand_dir
#     spline    = get_model_Keq_log_spline(model_ion_str_list, model_Keq)
#     df_exp    = get_exp_data(protein, resin)
#
#     res_vec = np.zeros(len(df_exp['IS(M)']))
#     df_exp['model_Keq'] = use_model_Keq_log_spline(spline, df_exp['IS(M)'])
#     df_exp['Keq_res']   = np.log(df_exp['Keq']/df_exp['model_Keq'])
#     return np.array(df_exp['Keq_res'])

def get_optimal_z0_v2(args):
    (model_ion_str_list, protein, resin, integrand_dir) = args
    guess = 0.4e-9
    fit = scipy.optimize.least_squares(get_Keq_residual_v2, guess, method='lm',
                                       args=(model_ion_str_list, protein,
                                       resin, integrand_dir))
    return [protein, resin, fit.x[0]]

def get_Keq_residual_vz(z0, model_ion_str_list, protein, resin, integrand_dir,
                        protein_exp):
    model_Keq = get_model_Keq_v2(ion_str_list, integrand_df_file, z0, integrand_dir)
    spline    = get_model_Keq_log_spline(model_ion_str_list, model_Keq)
    df_exp    = get_exp_data(protein_exp, resin)

    res_vec = np.zeros(len(df_exp['IS(M)']))
    df_exp['model_Keq'] = use_model_Keq_log_spline(spline, df_exp['IS(M)'])
    df_exp['Keq_res']   = np.log(df_exp['Keq']/df_exp['model_Keq'])
    return np.array(df_exp['Keq_res'])

def get_optimal_z0_vz(args):
    (model_ion_str_list, protein, resin, integrand_dir, protein_exp) = args
    guess = 0.3e-9
    fit = scipy.optimize.least_squares(get_Keq_residual_vz, guess, method='lm',
                                       args=(model_ion_str_list, protein,
                                       resin, integrand_dir, protein_exp))
    return [protein, resin, fit.x[0]]

################################################################################
# Compute integral contributions from each orientation

def get_integral_contribs(integrand_df_file, z0):
    integrand_df  = pd.read_csv(integrand_df_file, index_col=0)
    zs            = [float(i) for i in integrand_df.columns]
    first_z       = zs[0]
    last_z        = zs[-1]
    zs            = np.array(zs)
    integral_vals = []

    for i, cont in integrand_df.iterrows():
        (z_t, fun_1_log, fun_2) = get_integrand_interp_funs(zs, np.array(cont))
        integral_vals.append(eval_Keq_integral(z0, last_z, first_z, z_t, fun_1_log, fun_2))
    return integral_vals

def get_and_save_integral_contribs(args):
    (resin, protein, elements_df, integrand_df_file, z0, new_file) = args
    vals = get_integral_contribs(integrand_df_file, z0)
    df_new = pd.DataFrame()
    df_new['potential_V'] = elements_df['potential_V']
    df_new['integral'] = vals
    df_new.to_csv(new_file, index=False)
    return

















################################################################################
