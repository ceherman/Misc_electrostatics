import numpy as np
import pandas as pd
from scipy import interpolate, integrate

from modules_protein_structure.nd_ps_misc import *
################################################################################
# Get protein surface data

def get_apbs_data(dir_apbs, protein_name):
    df = pd.read_csv(dir_apbs + protein_name + '/surf_pot.csv',
                    names=['x_A', 'y_A', 'z_A', 'phi_V'])
    df['phi_V'] *= constants().kT/constants().e0
    return df

def get_tabi_data(dir_tabi, protein_name):
    df = pd.read_csv(dir_tabi + protein_name + '/rslt_surface_potential.csv',
                    names=['node_index', 'x_A', 'y_A', 'z_A', 'norm_x', 'norm_y',
                            'norm_z', 'phi_V', 'norm_phi'])
    df['phi_V'] *= 0.0433699
    return df

def get_link(dir_link, protein_name):
    df = pd.read_csv(dir_link + protein_name + '/rslt_node_links.csv', names=['node_1',
                    'node_2', 'node_3'])
    for col in df.columns:
        df[col] -= 1
    return df

def get_elements(df_link, df_data, df_tabi=None):
    max_pot = 0.5

    for ind, contents in df_link.iterrows():
        df_link.at[ind, 'x1_A'] = df_data.at[contents['node_1'], 'x_A']
        df_link.at[ind, 'y1_A'] = df_data.at[contents['node_1'], 'y_A']
        df_link.at[ind, 'z1_A'] = df_data.at[contents['node_1'], 'z_A']
        df_link.at[ind, 'x2_A'] = df_data.at[contents['node_2'], 'x_A']
        df_link.at[ind, 'y2_A'] = df_data.at[contents['node_2'], 'y_A']
        df_link.at[ind, 'z2_A'] = df_data.at[contents['node_2'], 'z_A']
        df_link.at[ind, 'x3_A'] = df_data.at[contents['node_3'], 'x_A']
        df_link.at[ind, 'y3_A'] = df_data.at[contents['node_3'], 'y_A']
        df_link.at[ind, 'z3_A'] = df_data.at[contents['node_3'], 'z_A']

        phi_1 = df_data.at[contents['node_1'], 'phi_V']
        phi_2 = df_data.at[contents['node_2'], 'phi_V']
        phi_3 = df_data.at[contents['node_3'], 'phi_V']

        df_link.at[ind, 'phi_V'] = np.array([phi_1, phi_2, phi_3]).mean()
        if df_link.at[ind, 'phi_V'] > max_pot:
            df_link.at[ind, 'phi_V'] = max_pot

        df_link.at[ind, 'x_cent_A'] = np.mean(np.array([df_link.at[ind, 'x1_A'],
                                                        df_link.at[ind, 'x2_A'],
                                                        df_link.at[ind, 'x3_A']]))
        df_link.at[ind, 'y_cent_A'] = np.mean(np.array([df_link.at[ind, 'y1_A'],
                                                        df_link.at[ind, 'y2_A'],
                                                        df_link.at[ind, 'y3_A']]))
        df_link.at[ind, 'z_cent_A'] = np.mean(np.array([df_link.at[ind, 'z1_A'],
                                                        df_link.at[ind, 'z2_A'],
                                                        df_link.at[ind, 'z3_A']]))
        if df_tabi is not None:
            n_x1 = df_tabi.at[contents['node_1'], 'norm_x']
            n_y1 = df_tabi.at[contents['node_1'], 'norm_y']
            n_z1 = df_tabi.at[contents['node_1'], 'norm_z']
            n_x2 = df_tabi.at[contents['node_2'], 'norm_x']
            n_y2 = df_tabi.at[contents['node_2'], 'norm_y']
            n_z2 = df_tabi.at[contents['node_2'], 'norm_z']
            n_x3 = df_tabi.at[contents['node_3'], 'norm_x']
            n_y3 = df_tabi.at[contents['node_3'], 'norm_y']
            n_z3 = df_tabi.at[contents['node_3'], 'norm_z']

            df_link.at[ind, 'norm_cent_x_A'] = np.array([n_x1, n_x2, n_x3]).mean()
            df_link.at[ind, 'norm_cent_y_A'] = np.array([n_y1, n_y2, n_y3]).mean()
            df_link.at[ind, 'norm_cent_z_A'] = np.array([n_z1, n_z2, n_z3]).mean()

            norm_length = np.sqrt(df_link.at[ind, 'norm_cent_x_A']**2 +\
                                  df_link.at[ind, 'norm_cent_y_A']**2 +\
                                  df_link.at[ind, 'norm_cent_z_A']**2)

            df_link.at[ind, 'norm_cent_x_A'] /= norm_length
            df_link.at[ind, 'norm_cent_y_A'] /= norm_length
            df_link.at[ind, 'norm_cent_z_A'] /= norm_length

        elif df_tabi is None:
            n_x1 = df_data.at[contents['node_1'], 'norm_x']
            n_y1 = df_data.at[contents['node_1'], 'norm_y']
            n_z1 = df_data.at[contents['node_1'], 'norm_z']
            n_x2 = df_data.at[contents['node_2'], 'norm_x']
            n_y2 = df_data.at[contents['node_2'], 'norm_y']
            n_z2 = df_data.at[contents['node_2'], 'norm_z']
            n_x3 = df_data.at[contents['node_3'], 'norm_x']
            n_y3 = df_data.at[contents['node_3'], 'norm_y']
            n_z3 = df_data.at[contents['node_3'], 'norm_z']

            df_link.at[ind, 'norm_cent_x_A'] = np.array([n_x1, n_x2, n_x3]).mean()
            df_link.at[ind, 'norm_cent_y_A'] = np.array([n_y1, n_y2, n_y3]).mean()
            df_link.at[ind, 'norm_cent_z_A'] = np.array([n_z1, n_z2, n_z3]).mean()

            norm_length = np.sqrt(df_link.at[ind, 'norm_cent_x_A']**2 +\
                                  df_link.at[ind, 'norm_cent_y_A']**2 +\
                                  df_link.at[ind, 'norm_cent_z_A']**2)

            df_link.at[ind, 'norm_cent_x_A'] /= norm_length
            df_link.at[ind, 'norm_cent_y_A'] /= norm_length
            df_link.at[ind, 'norm_cent_z_A'] /= norm_length

    return df_link

def translate_prot_centroid(df_prot):
    n_elements = len(df_prot['x_cent_A'])

    prot_x_cent = sum(df_prot['x_cent_A'])/n_elements
    prot_y_cent = sum(df_prot['y_cent_A'])/n_elements
    prot_z_cent = sum(df_prot['z_cent_A'])/n_elements

    df_prot['x_cent_A'] -= prot_x_cent
    df_prot['y_cent_A'] -= prot_y_cent
    df_prot['z_cent_A'] -= prot_z_cent
    return

def initialize_energy_columns(df_prot):
    df_prot['z_m']    = np.nan
    df_prot['f_plus'] = np.nan
    df_prot['f_J']    = np.nan
    return

def get_df_elem_apbs(dir_link, dir_apbs, protein_name, dir_tabi):
    ind_slice = protein_name.rfind('_')
    name_wo_is = protein_name[:ind_slice]

    df_link = get_link(dir_link, name_wo_is)
    df_data = get_apbs_data(dir_apbs, protein_name)
    df_tabi = get_tabi_data(dir_tabi, name_wo_is)
    df_prot = get_elements(df_link, df_data, df_tabi)
    translate_prot_centroid(df_prot)
    initialize_energy_columns(df_prot)
    return df_prot

def get_df_elem_tabi(dir_link, dir_tabi, protein_name):
    df_link = get_link(dir_link, protein_name)
    df_data = get_tabi_data(dir_tabi, protein_name)
    df_prot = get_elements(df_link, df_data)
    translate_prot_centroid(df_prot)
    initialize_energy_columns(df_prot)
    return df_prot


################################################################################
# Helper functions for performing 3D geometry calculations

def get_unit_vector(theta, phi):
    """Find unit vector pointed in the theta, phi direction."""
    a = np.sin(theta)*np.cos(phi)
    b = np.sin(theta)*np.sin(phi)
    c = np.cos(theta)
    return np.array([a, b, c])

def get_q_plane_for_sep(n, protein, min_sep_A):
    get_distance(n, [0, 0, 0], protein)
    translate_length = -1.0*min(protein.df['distance_A'])
    return (translate_length + min_sep_A)*n

def get_directions_in_plane(n):
    # Watch out for possible issue here
    v0 = np.array([0.0, 0.0, 1.0])
    if np.array_equal(n, v0) or np.array_equal(n, -1.0*v0):
        v0 = np.array([0.0, 1.0, 0.0])

    v1 = -1.0*np.cross(n, v0)/np.linalg.norm(np.cross(n, v0))
    v2 = np.cross(n, v1)/np.linalg.norm(np.cross(n, v1))
    return (v1, v2)

def get_distance_helper(n, q_plane, x, y, z):
    """Find the distance between a point (x, y, z) and a plane (defined by
    the unit normal n and the point q_plane)."""
    a  = n[0]
    b  = n[1]
    c  = n[2]
    x0 = q_plane[0]
    y0 = q_plane[1]
    z0 = q_plane[2]
    d  = a*x0 + b*y0 + c*z0
    distance = a*x + b*y + c*z - d
    return -1.0*distance

def get_projection_helper(n, q_plane, x, y, z):
    """Find the boundary element's projection into the resin plane."""
    a = n[0]
    b = n[1]
    c = n[2]
    d = q_plane[0]
    e = q_plane[1]
    f = q_plane[2]

    # Projected points px, py, pz
    t = (a*(d-x) + b*(e-y) + c*(f-z))/(a**2 + b**2 + c**2)
    (px, py, pz) = ((x + t*a), (y + t*b), (z + t*c))
    return (px, py, pz)

def get_triangle_area_helper(px1, py1, pz1, px2, py2, pz2, px3, py3, pz3):
    # Triangle side vectors
    (sx1, sy1, sz1) = (px2 - px1, py2 - py1, pz2 - pz1)
    (sx2, sy2, sz2) = (px3 - px1, py3 - py1, pz3 - pz1)

    # Cross product of the triangle sides
    cx = sy1*sz2 - sz1*sy2
    cy = sz1*sx2 - sx1*sz2
    cz = sx1*sy2 - sy1*sx2
    norm = np.sqrt(cx**2 + cy**2 + cz**2)
    return 0.5*norm

def get_vector_decomp_helper(a, x, y, z):
    b = np.array([x, y, z])
    x = np.linalg.solve(a, b)
    return x[0]


################################################################################
# 3D geometry calculations on full protein structure

def get_distance(n, q_plane, protein):
    protein.df['distance_A'] = get_distance_helper(n, q_plane,
    protein.df['x_cent_A'].values, protein.df['y_cent_A'].values,
    protein.df['z_cent_A'].values)

    protein.df['z_m'] = protein.df['distance_A']*1.0e-10
    return

def get_resin_projection(n, q_plane, protein):
    (protein.df['px1_A'], protein.df['py1_A'], protein.df['pz1_A']) =\
    get_projection_helper(n, q_plane, protein.df['x1_A'].values,
                          protein.df['y1_A'].values, protein.df['z1_A'].values)

    (protein.df['px2_A'], protein.df['py2_A'], protein.df['pz2_A']) =\
    get_projection_helper(n, q_plane, protein.df['x2_A'].values,
                          protein.df['y2_A'].values, protein.df['z2_A'].values)

    (protein.df['px3_A'], protein.df['py3_A'], protein.df['pz3_A']) =\
    get_projection_helper(n, q_plane, protein.df['x3_A'].values,
                          protein.df['y3_A'].values, protein.df['z3_A'].values)
    return

def get_projected_areas(protein):
    protein.df['area_A2'] = get_triangle_area_helper(protein.df['px1_A'].values,
    protein.df['py1_A'].values, protein.df['pz1_A'].values, protein.df['px2_A'].values,
    protein.df['py2_A'].values, protein.df['pz2_A'].values, protein.df['px3_A'].values,
    protein.df['py3_A'].values, protein.df['pz3_A'].values)
    return

def get_normal_vector_decomp_coeff(n, protein):
    (v1, v2) = get_directions_in_plane(n)
    a = np.column_stack((n, v1, v2))
    protein.df['decomp_n_coeff'] = get_vector_decomp_helper(a, protein.df['norm_cent_x_A'],
    protein.df['norm_cent_y_A'], protein.df['norm_cent_z_A'])
    return

def get_elements_facing_resin(protein):
    return protein.df[protein.df['decomp_n_coeff'] > 0.0]

def get_facing_projection(n, q_plane, protein):
    get_distance(n, q_plane, protein)
    get_resin_projection(n, q_plane, protein)
    get_projected_areas(protein)
    get_normal_vector_decomp_coeff(n, protein)

    # Remove elements not facing the resin
    facing_df = get_elements_facing_resin(protein)
    return facing_df.copy()

################################################################################
# Integrand functions

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

def get_integrand(energy):
    return np.exp(-1.0*energy/constants().kT) - 1.0

def unscale_energy(f_plus, area, bulk_solution):
    return f_plus*area*bulk_solution.ion_dens*constants().kT/bulk_solution.kappa

def get_integrand_list(args):
    (protein, index, theta, phi, zs, energy_fun, bulk_solution) = args

    n         = get_unit_vector(theta, phi)
    q_contact = get_q_plane_for_sep(n, protein, 0.0)
    facing_df = get_facing_projection(n, q_contact, protein)
    start_df  = facing_df.copy()

    integrand_res = np.zeros(len(zs))
    energy        = np.zeros(len(zs))

    for i, z in enumerate(zs):
        facing_df['z_m']    = start_df['z_m'] + z
        facing_df['f_plus'] = energy_fun.ev(facing_df['z_m'].values*bulk_solution.kappa,
                                            facing_df['phi_V'].values)
        facing_df['f_J']    = unscale_energy(facing_df['f_plus'].values,
                                             (facing_df['area_A2'].values)*1.0e-20,
                                             bulk_solution)
        energy[i]        = np.sum(facing_df['f_J'])
        integrand_res[i] = get_integrand(energy[i])
    return (index, integrand_res, energy)








################################################################################
# Integral functions to estimate Keq from the tabulated integrand results

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

def get_integrand_interp_funs_wrapper(protein, ion_str, resin):
    sol           = bulk_solution(ion_str)
    integrand_dir = 'integrand_results/'
    energy_dir    = 'energy_results/'

    # HACK - CHANGE FILE HERE TO PLOT RESULTS FROM DIFFERENT SYSTEMSS
    file = resin + '_' + protein + '_' + str(int(sol.ion_str*1.0e3)) + '_apbs_non.csv'

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

def get_model_Keq(ion_str_list, protein, resin, z0):
    Keq_results = np.zeros(len(ion_str_list))
    for j_index, ion_str in enumerate(ion_str_list):
        z_t, fun_1_log, fun_2, last_z, first_z = \
        get_integrand_interp_funs_wrapper(protein, ion_str, resin)
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

def get_Keq_residual(z0, model_ion_str_list, protein, resin):
    model_Keq = get_model_Keq(model_ion_str_list, protein, resin, z0)
    spline    = get_model_Keq_log_spline(model_ion_str_list, model_Keq)
    df_exp    = get_exp_data(protein, resin)

    res_vec = np.zeros(len(df_exp['IS(M)']))
    df_exp['model_Keq'] = use_model_Keq_log_spline(spline, df_exp['IS(M)'])
    df_exp['Keq_res']   = np.log(df_exp['Keq']/df_exp['model_Keq'])
    return np.array(df_exp['Keq_res'])

def get_optimal_z0(args):
    (model_ion_str_list, protein, resin) = args
    guess = 0.25e-9
    fit = scipy.optimize.least_squares(get_Keq_residual, guess, method='lm',
                                       args=(model_ion_str_list, protein,
                                        resin))
    return [protein, resin, fit.x[0]]



def get_integrand_interp_funs_wrapper_v2(file, ion_str, integrand_dir):
    integrand_df = pd.read_csv(integrand_dir + file, index_col=0)
    zs           = np.zeros(len(integrand_df.columns))
    averages     = np.zeros(len(integrand_df.columns))

    for k in enumerate(integrand_df.items()):
        (index, (z_str, integrands)) = k
        zs[index]       = float(z_str)
        averages[index] = np.average(integrands)

    z_t, fun_1_log, fun_2 = get_integrand_interp_funs(zs, averages)
    return z_t, fun_1_log, fun_2, zs[-1], zs[0]

def get_model_Keq_v2(ion_str_list, integrand_df_file, z0, integrand_dir):
    Keq_results = np.zeros(len(ion_str_list))
    for j_index, ion_str in enumerate(ion_str_list):
        z_t, fun_1_log, fun_2, last_z, first_z = \
        get_integrand_interp_funs_wrapper_v2(integrand_df_file, ion_str, integrand_dir)
        Keq_results[j_index] = eval_Keq_integral(z0, last_z, first_z, z_t, fun_1_log, fun_2)
    return Keq_results


################################################################################
# Boneyard


# df_link.at[ind, 'norm_x1_A'] = df_data.at[contents['node_1'], 'norm_x']
# df_link.at[ind, 'norm_y1_A'] = df_data.at[contents['node_1'], 'norm_y']
# df_link.at[ind, 'norm_z1_A'] = df_data.at[contents['node_1'], 'norm_z']
# df_link.at[ind, 'norm_x2_A'] = df_data.at[contents['node_2'], 'norm_x']
# df_link.at[ind, 'norm_y2_A'] = df_data.at[contents['node_2'], 'norm_y']
# df_link.at[ind, 'norm_z2_A'] = df_data.at[contents['node_2'], 'norm_z']
# df_link.at[ind, 'norm_x3_A'] = df_data.at[contents['node_3'], 'norm_x']
# df_link.at[ind, 'norm_y3_A'] = df_data.at[contents['node_3'], 'norm_y']
# df_link.at[ind, 'norm_z3_A'] = df_data.at[contents['node_3'], 'norm_z']
