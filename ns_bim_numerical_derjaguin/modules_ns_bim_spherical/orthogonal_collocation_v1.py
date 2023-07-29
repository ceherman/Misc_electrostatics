import numpy as np

from modules_ns_bim_spherical.geometry import *
from modules_ns_bim_spherical.hybrid_model_oc_v1 import *
from modules_ns_bim_spherical.oc_plots_v1 import *
from modules_ns_bim_spherical.numerical_Derjaguin_v1 import *

################################################################################
# Foundational functions

def get_b1(c1, kap_h, t):
    return c1*np.exp(-1.0*kap_h*t)

def get_b1_d1(c1, kap_h, t):
    return -1.0*kap_h*get_b1(c1, kap_h, t)

def get_b1_d2(c1, kap_h, t):
    return (kap_h**2)*get_b1(c1, kap_h, t)

def get_b2(c2, kap_h, t):
    return c2*np.exp(-1.0*kap_h * (1-t))

def get_b2_d1(c2, kap_h, t):
    return kap_h*get_b2(c2, kap_h, t)

def get_b2_d2(c2, kap_h, t):
    return (kap_h**2)*get_b2(c2, kap_h, t)

def get_b3(kap_h, t):
    return np.exp(-1.0*kap_h*t) + np.exp(-1.0*kap_h * (1-t))

def get_b3_d1(kap_h, t):
    return -1.0*kap_h*np.exp(-1.0*kap_h*t) + kap_h*np.exp(-1.0*kap_h * (1-t))

def get_b3_d2(kap_h, t):
    return (kap_h**2)*get_b3(kap_h, t)

def get_nu(coeff_vec, t):
    total = 0.0
    for i in range(len(coeff_vec) - 1):
        total += coeff_vec[i]*np.cos(np.pi*i*t)
    return total

def get_nu_d1(coeff_vec, t):
    total = 0.0
    for i in range(len(coeff_vec) - 1):
        total += coeff_vec[i]*i*np.sin(np.pi*i*t)
    return -1.0*np.pi*total

def get_nu_d2(coeff_vec, t):
    total = 0.0
    for i in range(len(coeff_vec) - 1):
        total += coeff_vec[i]*(i**2)*np.cos(np.pi*i*t)
    return -1.0*(np.pi**2)*total

def get_yp(coeff_vec, t, kap_h, c2):
    c1 = coeff_vec[-1]
    b1 = get_b1(c1, kap_h, t)
    b2 = get_b2(c2, kap_h, t)
    b3 = get_b3(kap_h, t)
    nu = get_nu(coeff_vec, t)
    return ((1.0-t)*b1 + t*b2 + t*(1.0-t)*b3*nu)

def get_yp_d1(coeff_vec, t, kap_h, c2):
    c1 = coeff_vec[-1]
    b1 = get_b1(c1, kap_h, t)
    b2 = get_b2(c2, kap_h, t)
    b3 = get_b3(kap_h, t)
    nu = get_nu(coeff_vec, t)

    b1_d1 = get_b1_d1(c1, kap_h, t)
    b2_d1 = get_b2_d1(c2, kap_h, t)
    b3_d1 = get_b3_d1(kap_h, t)
    nu_d1 = get_nu_d1(coeff_vec, t)

    value = -1.0*b1 + (1.0-t)*b1_d1 + b2 + t*b2_d1 + (1.0-2.0*t)*b3*nu \
            + t*(1.0-t)*b3_d1*nu + t*(1.0-t)*b3*nu_d1
    return value

def get_yp_d2(coeff_vec, t, kap_h, c2):
    c1 = coeff_vec[-1]
    b3 = get_b3(kap_h, t)
    nu = get_nu(coeff_vec, t)

    b1_d1 = get_b1_d1(c1, kap_h, t)
    b2_d1 = get_b2_d1(c2, kap_h, t)
    b3_d1 = get_b3_d1(kap_h, t)
    nu_d1 = get_nu_d1(coeff_vec, t)

    b1_d2 = get_b1_d2(c1, kap_h, t)
    b2_d2 = get_b2_d2(c2, kap_h, t)
    b3_d2 = get_b3_d2(kap_h, t)
    nu_d2 = get_nu_d2(coeff_vec, t)

    value = -2.0*b1_d1 + (1.0-t)*b1_d2 + 2.0*b2_d1 + t*b2_d2 - 2.0*b3*nu \
            + (2.0-4.0*t)*b3_d1*nu + (2.0-4.0*t)*b3*nu_d1 + t*(1.0-t)*b3_d2*nu \
            + 2.0*t*(1.0-t)*b3_d1*nu_d1 + t*(1.0-t)*b3*nu_d2
    return value

################################################################################
# Interface for parameter scaling and retrieval

def get_t_from_z(z, h):
    return z/h

def get_z_from_t(t, h):
    return h*t

def get_t_from_z_prime(z_prime, kap_h):
    return z_prime/kap_h

def get_psi_from_phi(phi):
    """Compute scaled surface potential."""
    kT = constants().kT
    e0 = constants().e0
    return e0*phi/kT

def get_phi_from_psi(psi):
    """Compute surface potential [V]."""
    kT = constants().kT
    e0 = constants().e0
    return kT*psi/e0

def get_sigma_star_from_sigma(sigma, bulk_solution):
    """Compute scaled surface charge density."""
    kappa = bulk_solution.kappa
    n0    = bulk_solution.ion_dens
    e0    = constants().e0
    return sigma*kappa/(2.0*n0*e0)

def get_sigma_from_sigma_star(sigma_star, bulk_solution):
    """Compute surface charge density [C m-2]."""
    kappa = bulk_solution.kappa
    n0    = bulk_solution.ion_dens
    e0    = constants().e0
    return 2.0*n0*e0*sigma_star/kappa

def get_y_from_psi(psi):
    return np.tanh(psi/4.0)

def get_psi_from_y(y):
    return 4.0*np.arctanh(y)

def unscale_energy(f_plus, area, bulk_solution):
    return f_plus*area*bulk_solution.ion_dens*constants().kT/bulk_solution.kappa

################################################################################
# Class for storing information related to boundary conditions

class cp_bc_info:
    # Information related to mixed boundary conditions (of consant charge
    # density at surface 1 and constant potential at surface 2)
    def __init__(self, sigma_01, phi_02, h, bulk_solution):
        self.sigma_01      = sigma_01
        self.sigma_01_star = get_sigma_star_from_sigma(sigma_01, bulk_solution)
        self.phi_02        = phi_02
        self.psi_02        = get_psi_from_phi(phi_02)
        self.c2            = get_y_from_psi(self.psi_02)
        self.h             = h
        self.bulk_solution = bulk_solution
        self.kap_h         = h*bulk_solution.kappa
        return

    def get_sigma_02_star(self, coeff_vec):
        t     = 1.0
        yp    = get_yp(coeff_vec, t, self.kap_h, self.c2)
        yp_d1 = get_yp_d1(coeff_vec, t, self.kap_h, self.c2)
        self.sigma_02_star = 4.0*yp_d1/(self.kap_h*(1.0-yp**2))
        return

    def get_psi_01(self, coeff_vec):
        t = 0.0
        self.psi_01 = get_psi_at_t(t, coeff_vec, self)
        return

    def update_pc_bc_info(self, coeff_vec):
        self.get_sigma_02_star(coeff_vec)
        self.get_psi_01(coeff_vec)
        return

class cc_bc_info:
    # Information related to mixed boundary conditions (of consant charge
    # density at both surfaces)
    def __init__(self, sigma_01, sigma_02, h, bulk_solution):
        self.sigma_01      = sigma_01
        self.sigma_02      = sigma_02
        self.sigma_01_star = get_sigma_star_from_sigma(sigma_01, bulk_solution)
        self.sigma_02_star = get_sigma_star_from_sigma(sigma_02, bulk_solution)
        self.h             = h
        self.bulk_solution = bulk_solution
        self.kap_h         = h*bulk_solution.kappa
        return

    def get_psi_01(self, coeff_vec):
        t = 0.0
        self.psi_01 = get_psi_at_t_v3(t, coeff_vec, self)
        return

    def get_psi_02(self, coeff_vec):
        t = 1.0
        self.psi_02 = get_psi_at_t_v3(t, coeff_vec, self)
        return

    def update_cc_bc_info(self, coeff_vec):
        self.get_psi_01(coeff_vec)
        self.get_psi_02(coeff_vec)
        return

################################################################################
# Functions for solving the problem

def get_c1_initial_guess(bc_info, bulk_solution):
    phi = inverse_Grahame(bc_info.sigma_01, bulk_solution) # NLPBE
    # phi = inverse_capacitor(bc_info.sigma_01, bulk_solution) # LPBE
    psi = get_psi_from_phi(phi)
    y   = get_y_from_psi(psi)
    return y

def scaled_1d_pb_residual(coeff_vec, t, bc_info):
    yp    = get_yp(coeff_vec, t, bc_info.kap_h, bc_info.c2)
    yp_d1 = get_yp_d1(coeff_vec, t, bc_info.kap_h, bc_info.c2)
    yp_d2 = get_yp_d2(coeff_vec, t, bc_info.kap_h, bc_info.c2)
    if -1.0 < yp < 1.0:
        psi   = get_psi_from_y(yp)
        value = 2.0*yp*yp_d1**2/((1.0-yp**2)**2) + yp_d2/(1.0-yp**2) \
                - (bc_info.kap_h**2)/4.0*np.sinh(psi) # NLPBE
        # value = 2.0*yp*yp_d1**2/((1.0-yp**2)**2) + yp_d2/(1.0-yp**2) \
        #         - (bc_info.kap_h**2)*np.arctanh(yp) # LPBE
        return value
    else:
        return 1.0e20 # To keep the solver from finding solutions where abs(yp) > 1

def scaled_bc_residual(coeff_vec, bc_info):
    t     = 0.0
    yp    = get_yp(coeff_vec, t, bc_info.kap_h, bc_info.c2)
    yp_d1 = get_yp_d1(coeff_vec, t, bc_info.kap_h, bc_info.c2)
    return yp_d1 + bc_info.kap_h*bc_info.sigma_01_star/4.0*(1.0-yp**2)

def get_residual_vector(coeff_vec, bc_info):
    p = len(coeff_vec) - 1
    residual_vec = np.zeros(len(coeff_vec))

    for i in range(p):
        t_i = 1.0/(2.0*(p)) + i/(p)
        residual_vec[i] = scaled_1d_pb_residual(coeff_vec, t_i, bc_info)

    residual_vec[-1] = scaled_bc_residual(coeff_vec, bc_info)
    return residual_vec

def get_coeff_vec(guess, bc_info):
    fit = scipy.optimize.least_squares(get_residual_vector, guess, method='lm',
                                       args=(bc_info,)) # xtol=1.0e-10, ftol=1.0e-10
    return fit



################################################################################
# Calculating interaction energy from the problem solution

def get_psi_at_t(t, coeff_vec, bc_info):
    yp = get_yp(coeff_vec, t, bc_info.kap_h, bc_info.c2)
    return get_psi_from_y(yp)

# Same as getting psi at z
def get_psi_at_z_prime(z_prime, coeff_vec, bc_info):
    t = get_t_from_z_prime(z_prime, bc_info.kap_h)
    return get_psi_at_t(t, coeff_vec, bc_info)

def sinh_integrand(z_prime, coeff_vec, bc_info):
    psi = get_psi_at_z_prime(z_prime, coeff_vec, bc_info)
    return psi*np.sinh(psi)

def cosh_integrand(z_prime, coeff_vec, bc_info):
    psi = get_psi_at_z_prime(z_prime, coeff_vec, bc_info)
    return np.cosh(psi) - 1.0

def get_integral(integrand_fun, coeff_vec, bc_info):
    value = scipy.integrate.quad(integrand_fun, 0.0, bc_info.kap_h,
                                 args=(coeff_vec, bc_info))
    return value

def get_scaled_u_el(bc_info, sinh_integral):
    value = bc_info.sigma_01_star*bc_info.psi_01 + \
          + bc_info.sigma_02_star*bc_info.psi_02 - sinh_integral
    return value

def get_scaled_ts(sinh_integral, cosh_integral):
    return 2.0*(sinh_integral - cosh_integral)

def get_scaled_f_chem(bc_info):
    return -1.0*bc_info.sigma_02_star*bc_info.psi_02

def get_scaled_energy(coeff_vec, bc_info):
    """Calculate the scaled (energy per unit area)"""
    bc_info.update_pc_bc_info(coeff_vec)

    int_rel_tol   = 1.0e-3
    sinh_integral = get_integral(sinh_integrand, coeff_vec, bc_info)
    cosh_integral = get_integral(cosh_integrand, coeff_vec, bc_info)
    assert (sinh_integral[1]/sinh_integral[0] < int_rel_tol)
    assert (cosh_integral[1]/cosh_integral[0] < int_rel_tol)

    u_el   = get_scaled_u_el(bc_info, sinh_integral[0])
    ts     = get_scaled_ts(sinh_integral[0], cosh_integral[0])
    f_chem = get_scaled_f_chem(bc_info)
    return u_el + ts + f_chem


################################################################################
# Combined function to solve the problem and compute interaction energy

def solve_for_energy(bc_info, guess):
    cost_tol = 1.0e-3
    fit      = get_coeff_vec(guess, bc_info)
    energy   = get_scaled_energy(fit.x, bc_info)
    assert(fit.cost < cost_tol)
    return (energy, fit.x)

def get_ref_energy(sigma_01, phi_02, bulk_solution):
    h         = 10.0/sol.kappa  # 10 DL separation
    bc_info   = cp_bc_info(sigma_01, phi_02, h, bulk_solution)
    guess     = np.full(41, 0.0)
    guess[-1] = get_c1_initial_guess(bc_info, sol)
    (energy, coeff_vec) = solve_for_energy(bc_info, guess)
    return energy

def get_energy_for_seps(args):
    (sigma_01, phi_02, bulk_solution, seps) = args

    energies = []
    for h_d in seps:
        h       = h_d/bulk_solution.kappa
        bc_info = cp_bc_info(sigma_01, phi_02, h, bulk_solution)
        if h_d == seps[0]:
            guess     = np.full(41, 0.0)
            guess[-1] = get_c1_initial_guess(bc_info, sol)
        else:
            pass
        (energy, guess) = solve_for_energy(bc_info, guess)
        energies.append(energy)

    ref = get_ref_energy(sigma_01, phi_02, sol)
    for i in range(len(energies)):
        energies[i] -= ref

    results = []
    for i in range(len(seps)):
        results.append((seps[i], phi_02, energies[i]))

    return results


################################################################################
# v2 of the combined functions to solve the problem and compute the interaction
# energy, saving the individual contributions to the potential of mean force
def get_scaled_u_el_v2(bc_info, sinh_integral):
    term_1 = bc_info.sigma_01_star*bc_info.psi_01
    term_2 = bc_info.sigma_02_star*bc_info.psi_02
    term_3 = -1.0*sinh_integral
    return (term_1, term_2, term_3)

def get_scaled_energy_v2(coeff_vec, bc_info):
    """Calculate the scaled (energy per unit area)"""
    bc_info.update_pc_bc_info(coeff_vec)

    int_rel_tol   = 1.0e-3
    sinh_integral = get_integral(sinh_integrand, coeff_vec, bc_info)
    cosh_integral = get_integral(cosh_integrand, coeff_vec, bc_info)
    assert (sinh_integral[1]/sinh_integral[0] < int_rel_tol)
    assert (cosh_integral[1]/cosh_integral[0] < int_rel_tol)

    (t1_elec, t2_elec, t3_elec) = get_scaled_u_el_v2(bc_info, sinh_integral[0])
    t4_entr                      = get_scaled_ts(sinh_integral[0], cosh_integral[0])
    return (t1_elec, t2_elec, t3_elec, t4_entr)

def solve_for_energy_v2(bc_info, guess):
    cost_tol           = 1.0e-3
    fit                = get_coeff_vec(guess, bc_info)
    (t1_elec, t2_elec, t3_elec, t4_entr) = get_scaled_energy_v2(fit.x, bc_info)
    assert(fit.cost < cost_tol)
    return (t1_elec, t2_elec, t3_elec, t4_entr, fit.x)

def get_ref_energy_v2(sigma_01, phi_02, bulk_solution):
    h         = 10.0/sol.kappa  # 10 DL separation
    bc_info   = cp_bc_info(sigma_01, phi_02, h, bulk_solution)
    guess     = np.full(41, 0.0)
    guess[-1] = get_c1_initial_guess(bc_info, sol)
    (t1_elec, t2_elec, t3_elec, t4_entr, coeff_vec) = solve_for_energy_v2(bc_info, guess)
    return (t1_elec, t2_elec, t3_elec, t4_entr)

def get_energy_for_seps_v2(args):
    (sigma_01, phi_02, bulk_solution, seps) = args

    energy_t1 = []
    energy_t2 = []
    energy_t3 = []
    energy_t4 = []

    for h_d in seps:
        h       = h_d/bulk_solution.kappa
        bc_info = cp_bc_info(sigma_01, phi_02, h, bulk_solution)
        if h_d == seps[0]:
            guess     = np.full(41, 0.0)
            guess[-1] = get_c1_initial_guess(bc_info, sol)
        else:
            pass
        (t1, t2, t3, t4, guess) = solve_for_energy_v2(bc_info, guess)

        energy_t1.append(t1)
        energy_t2.append(t2)
        energy_t3.append(t3)
        energy_t4.append(t4)

    (t1_ref, t2_ref, t3_ref, t4_ref) = get_ref_energy_v2(sigma_01, phi_02, sol)
    for i in range(len(energy_t1)):
        energy_t1[i] -= t1_ref
        energy_t2[i] -= t2_ref
        energy_t3[i] -= t3_ref
        energy_t4[i] -= t4_ref

    results = []
    for i in range(len(seps)):
        results.append((seps[i], phi_02, energy_t1[i], energy_t2[i],
                        energy_t3[i], energy_t4[i]))
    return results
















































################################################################################
# Functions for solving the problem with cc BC

def get_yp_cc(coeff_vec, t, kap_h):
    c1 = coeff_vec[-2]
    c2 = coeff_vec[-1]

    b1 = get_b1(c1, kap_h, t)
    b2 = get_b2(c2, kap_h, t)
    b3 = get_b3(kap_h, t)
    nu = get_nu(coeff_vec, t)
    return ((1.0-t)*b1 + t*b2 + t*(1.0-t)*b3*nu)

def get_yp_d1_cc(coeff_vec, t, kap_h):
    c1 = coeff_vec[-2]
    c2 = coeff_vec[-1]

    b1 = get_b1(c1, kap_h, t)
    b2 = get_b2(c2, kap_h, t)
    b3 = get_b3(kap_h, t)
    nu = get_nu(coeff_vec, t)

    b1_d1 = get_b1_d1(c1, kap_h, t)
    b2_d1 = get_b2_d1(c2, kap_h, t)
    b3_d1 = get_b3_d1(kap_h, t)
    nu_d1 = get_nu_d1(coeff_vec, t)

    value = -1.0*b1 + (1.0-t)*b1_d1 + b2 + t*b2_d1 + (1.0-2.0*t)*b3*nu \
            + t*(1.0-t)*b3_d1*nu + t*(1.0-t)*b3*nu_d1
    return value

def get_yp_d2_cc(coeff_vec, t, kap_h):
    c1 = coeff_vec[-2]
    c2 = coeff_vec[-1]

    b3 = get_b3(kap_h, t)
    nu = get_nu(coeff_vec, t)

    b1_d1 = get_b1_d1(c1, kap_h, t)
    b2_d1 = get_b2_d1(c2, kap_h, t)
    b3_d1 = get_b3_d1(kap_h, t)
    nu_d1 = get_nu_d1(coeff_vec, t)

    b1_d2 = get_b1_d2(c1, kap_h, t)
    b2_d2 = get_b2_d2(c2, kap_h, t)
    b3_d2 = get_b3_d2(kap_h, t)
    nu_d2 = get_nu_d2(coeff_vec, t)

    value = -2.0*b1_d1 + (1.0-t)*b1_d2 + 2.0*b2_d1 + t*b2_d2 - 2.0*b3*nu \
            + (2.0-4.0*t)*b3_d1*nu + (2.0-4.0*t)*b3*nu_d1 + t*(1.0-t)*b3_d2*nu \
            + 2.0*t*(1.0-t)*b3_d1*nu_d1 + t*(1.0-t)*b3*nu_d2
    return value

def scaled_1d_pb_residual_cc(coeff_vec, t, bc_info):
    yp    = get_yp_cc(coeff_vec, t, bc_info.kap_h)
    yp_d1 = get_yp_d1_cc(coeff_vec, t, bc_info.kap_h)
    yp_d2 = get_yp_d2_cc(coeff_vec, t, bc_info.kap_h)
    if -1.0 < yp < 1.0:
        psi   = get_psi_from_y(yp)
        value = 2.0*yp*yp_d1**2/((1.0-yp**2)**2) + yp_d2/(1.0-yp**2) \
                - (bc_info.kap_h**2)/4.0*np.sinh(psi)
        return value
    else:
        return 1.0e20 # To keep the solver from finding solutions where abs(yp) > 1

def scaled_bc_residual_cc(coeff_vec, bc_info):
    t     = 0.0
    yp    = get_yp_cc(coeff_vec, t, bc_info.kap_h)
    yp_d1 = get_yp_d1_cc(coeff_vec, t, bc_info.kap_h)
    res_1 = yp_d1 + bc_info.kap_h*bc_info.sigma_01_star/4.0*(1.0-yp**2)

    t     = 1.0
    yp    = get_yp_cc(coeff_vec, t, bc_info.kap_h)
    yp_d1 = get_yp_d1_cc(coeff_vec, t, bc_info.kap_h)
    res_2 = yp_d1 - bc_info.kap_h*bc_info.sigma_02_star/4.0*(1.0-yp**2)
    return res_1, res_2

def get_residual_vector_cc(coeff_vec, bc_info):
    p = len(coeff_vec) - 2
    residual_vec = np.zeros(len(coeff_vec))

    for i in range(p):
        t_i = 1.0/(2.0*(p)) + i/(p)
        residual_vec[i] = scaled_1d_pb_residual_cc(coeff_vec, t_i, bc_info)

    residual_vec[-2], residual_vec[-1] = scaled_bc_residual_cc(coeff_vec, bc_info)
    return residual_vec

def get_coeff_vec_cc(guess, bc_info):
    fit = scipy.optimize.least_squares(get_residual_vector_cc, guess, method='lm',
                                       args=(bc_info,)) # xtol=1.0e-10, ftol=1.0e-10
    return fit
#
#
#
#
#
#
#
#
#
#
#
#
#
# # ################################################################################
# # # Functions for solving the problem with cc BC
# #
# # def get_yp_cc(coeff_vec, t, kap_h):
# #     c1 = coeff_vec[-1]
# #     c2 = coeff_vec[-2]
# #
# #     b1 = get_b1(c1, kap_h, t)
# #     b2 = get_b2(c2, kap_h, t)
# #     b3 = get_b3(kap_h, t)
# #     nu = get_nu(coeff_vec, t)
# #     return ((1.0-t)*b1 + t*b2 + t*(1.0-t)*b3*nu)
# #
# # def get_yp_d1_cc(coeff_vec, t, kap_h):
# #     c1 = coeff_vec[-1]
# #     c2 = coeff_vec[-2]
# #
# #     b1 = get_b1(c1, kap_h, t)
# #     b2 = get_b2(c2, kap_h, t)
# #     b3 = get_b3(kap_h, t)
# #     nu = get_nu(coeff_vec, t)
# #
# #     b1_d1 = get_b1_d1(c1, kap_h, t)
# #     b2_d1 = get_b2_d1(c2, kap_h, t)
# #     b3_d1 = get_b3_d1(kap_h, t)
# #     nu_d1 = get_nu_d1(coeff_vec, t)
# #
# #     value = -1.0*b1 + (1.0-t)*b1_d1 + b2 + t*b2_d1 + (1.0-2.0*t)*b3*nu \
# #             + t*(1.0-t)*b3_d1*nu + t*(1.0-t)*b3*nu_d1
# #     return value
# #
# # def get_yp_d2_cc(coeff_vec, t, kap_h):
# #     c1 = coeff_vec[-1]
# #     c2 = coeff_vec[-2]
# #
# #     b3 = get_b3(kap_h, t)
# #     nu = get_nu(coeff_vec, t)
# #
# #     b1_d1 = get_b1_d1(c1, kap_h, t)
# #     b2_d1 = get_b2_d1(c2, kap_h, t)
# #     b3_d1 = get_b3_d1(kap_h, t)
# #     nu_d1 = get_nu_d1(coeff_vec, t)
# #
# #     b1_d2 = get_b1_d2(c1, kap_h, t)
# #     b2_d2 = get_b2_d2(c2, kap_h, t)
# #     b3_d2 = get_b3_d2(kap_h, t)
# #     nu_d2 = get_nu_d2(coeff_vec, t)
# #
# #     value = -2.0*b1_d1 + (1.0-t)*b1_d2 + 2.0*b2_d1 + t*b2_d2 - 2.0*b3*nu \
# #             + (2.0-4.0*t)*b3_d1*nu + (2.0-4.0*t)*b3*nu_d1 + t*(1.0-t)*b3_d2*nu \
# #             + 2.0*t*(1.0-t)*b3_d1*nu_d1 + t*(1.0-t)*b3*nu_d2
# #     return value
# #
# # def scaled_1d_pb_residual_cc(coeff_vec, t, bc_info):
# #     yp    = get_yp_cc(coeff_vec, t, bc_info.kap_h)
# #     yp_d1 = get_yp_d1_cc(coeff_vec, t, bc_info.kap_h)
# #     yp_d2 = get_yp_d2_cc(coeff_vec, t, bc_info.kap_h)
# #     if -1.0 < yp < 1.0:
# #         psi   = get_psi_from_y(yp)
# #         value = 2.0*yp*yp_d1**2/((1.0-yp**2)**2) + yp_d2/(1.0-yp**2) \
# #                 - (bc_info.kap_h**2)/4.0*np.sinh(psi) # NLPBE
# #         # value = 2.0*yp*yp_d1**2/((1.0-yp**2)**2) + yp_d2/(1.0-yp**2) \
# #         #         - (bc_info.kap_h**2)*np.arctanh(yp) # LPBE
# #         return value
# #     else:
# #         return 1.0e20 # To keep the solver from finding solutions where abs(yp) > 1
# #
# # def scaled_bc_residual_cc(coeff_vec, bc_info):
# #     t     = 0.0
# #     yp    = get_yp_cc(coeff_vec, t, bc_info.kap_h)
# #     yp_d1 = get_yp_d1_cc(coeff_vec, t, bc_info.kap_h)
# #     res_1 = yp_d1 + bc_info.kap_h*bc_info.sigma_01_star/4.0*(1.0-yp**2)
# #
# #     t     = 1.0
# #     yp    = get_yp_cc(coeff_vec, t, bc_info.kap_h)
# #     yp_d1 = get_yp_d1_cc(coeff_vec, t, bc_info.kap_h)
# #     res_2 = yp_d1 - bc_info.kap_h*bc_info.sigma_02_star/4.0*(1.0-yp**2)
# #     return res_1, res_2
# #
# # def get_residual_vector_cc(coeff_vec, bc_info):
# #     p = len(coeff_vec) - 2
# #     residual_vec = np.zeros(len(coeff_vec))
# #
# #     for i in range(p):
# #         t_i = 1.0/(2.0*(p)) + i/(p)
# #         residual_vec[i] = scaled_1d_pb_residual_cc(coeff_vec, t_i, bc_info)
# #
# #     residual_vec[-1], residual_vec[-2] = scaled_bc_residual_cc(coeff_vec, bc_info)
# #     return residual_vec
# #
# # def get_coeff_vec_cc_v3(guess, bc_info):
# #     fit = scipy.optimize.least_squares(get_residual_vector_cc_v3, guess, method='lm',
# #                                        args=(bc_info,)) # xtol=1.0e-10, ftol=1.0e-10
# #     return fit
# #
# #
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
# ################################################################################
# # v3 of the combined functions to solve the problem with constant protein
# # surface charge density (cc BC). Individual contributions to the potential of
# # mean force are saved
#
# def get_c1_c2_initial_guess(bc_info, bulk_solution):
#     phi_1 = inverse_Grahame(bc_info.sigma_01, bulk_solution) # NLPBE
#     psi_1 = get_psi_from_phi(phi_1)
#     y_1   = get_y_from_psi(psi_1)
#
#     phi_2 = inverse_Grahame(bc_info.sigma_02, bulk_solution) # NLPBE
#     psi_2 = get_psi_from_phi(phi_2)
#     y_2   = get_y_from_psi(psi_2)
#     return y_1, y_2
#
# def get_psi_at_t_v3(t, coeff_vec, bc_info):
#     yp = get_yp_cc(coeff_vec, t, bc_info.kap_h)
#     return get_psi_from_y(yp)
#
# def get_psi_at_z_prime_v3(z_prime, coeff_vec, bc_info):
#     t = get_t_from_z_prime(z_prime, bc_info.kap_h)
#     return get_psi_at_t_v3(t, coeff_vec, bc_info)
#
# def sinh_integrand_v3(z_prime, coeff_vec, bc_info):
#     psi = get_psi_at_z_prime_v3(z_prime, coeff_vec, bc_info)
#     return psi*np.sinh(psi)
#
# def cosh_integrand_v3(z_prime, coeff_vec, bc_info):
#     psi = get_psi_at_z_prime_v3(z_prime, coeff_vec, bc_info)
#     return np.cosh(psi) - 1.0
#
# def get_scaled_energy_v3(coeff_vec, bc_info):
#     """Calculate the scaled (energy per unit area)"""
#     bc_info.update_cc_bc_info(coeff_vec)
#
#     int_rel_tol   = 1.0e-3
#     sinh_integral = get_integral(sinh_integrand_v3, coeff_vec, bc_info)
#     cosh_integral = get_integral(cosh_integrand_v3, coeff_vec, bc_info)
#     assert (sinh_integral[1]/sinh_integral[0] < int_rel_tol)
#     assert (cosh_integral[1]/cosh_integral[0] < int_rel_tol)
#
#     (t1_elec, t2_elec, t3_elec) = get_scaled_u_el_v2(bc_info, sinh_integral[0])
#     t4_entr                     = get_scaled_ts(sinh_integral[0], cosh_integral[0])
#     return (t1_elec, t2_elec, t3_elec, t4_entr)
#
# def solve_for_energy_v3(bc_info, guess):
#     cost_tol           = 1.0e-3
#     fit                = get_coeff_vec_cc(guess, bc_info)
#     (t1_elec, t2_elec, t3_elec, t4_entr) = get_scaled_energy_v3(fit.x, bc_info)
#     assert(fit.cost < cost_tol)
#     return (t1_elec, t2_elec, t3_elec, t4_entr, fit.x)
#
# def get_ref_energy_v3(sigma_01, sigma_02, bulk_solution):
#     h         = 10.0/sol.kappa  # 10 DL separation
#     bc_info   = cc_bc_info(sigma_01, sigma_02, h, bulk_solution)
#     guess     = np.full(42, 0.0)
#     guess[-2], guess[-1] = get_c1_c2_initial_guess(bc_info, sol)
#     (t1_elec, t2_elec, t3_elec, t4_entr, coeff_vec) = solve_for_energy_v3(bc_info, guess)
#     return (t1_elec, t2_elec, t3_elec, t4_entr)
#
# def get_energy_for_seps_v3(args):
#     (sigma_01, sigma_02, bulk_solution, seps) = args
#     energy_t1 = []
#     energy_t2 = []
#     energy_t3 = []
#     energy_t4 = []
#
#     for h_d in seps:
#         h       = h_d/bulk_solution.kappa
#         bc_info = cc_bc_info(sigma_01, sigma_02, h, bulk_solution)
#         if h_d == seps[0]:
#             guess                = np.full(42, 0.0)
#             guess[-2], guess[-1] = get_c1_c2_initial_guess(bc_info, sol)
#         else:
#             pass
#         (t1, t2, t3, t4, guess) = solve_for_energy_v3(bc_info, guess)
#
#         energy_t1.append(t1)
#         energy_t2.append(t2)
#         energy_t3.append(t3)
#         energy_t4.append(t4)
#
#     (t1_ref, t2_ref, t3_ref, t4_ref) = get_ref_energy_v3(sigma_01, sigma_02, sol)
#     for i in range(len(energy_t1)):
#         energy_t1[i] -= t1_ref
#         energy_t2[i] -= t2_ref
#         energy_t3[i] -= t3_ref
#         energy_t4[i] -= t4_ref
#
#     results = []
#     for i in range(len(seps)):
#         results.append((seps[i], phi_02, energy_t1[i], energy_t2[i],
#                         energy_t3[i], energy_t4[i]))
#     return results



























# ################################################################################
# # Functions to compute the full PMF and the 'PMF' I understand to be used by
# # NS-BIM, both based on the analytical DH solution
#
# def dh_analyt_get_c2(bc_info):
#     num = bc_info.psi_02 + bc_info.sigma_01_star*np.exp(bc_info.kap_h)
#     den = np.exp(bc_info.kap_h) + np.exp(-1.0*bc_info.kap_h)
#     return num/den
#
# def dh_analyt_get_c1(bc_info):
#     c2 = dh_analyt_get_c2(bc_info)
#     return c2 - bc_info.sigma_01_star
#
# def dh_analyt_get_psi(t, bc_info):
#     c1 = dh_analyt_get_c1(bc_info)
#     c2 = dh_analyt_get_c2(bc_info)
#     return c1*np.exp(bc_info.kap_h*t) + c2*np.exp(-1.0*bc_info.kap_h*t)
#
# def dh_update_pc_bc_info(bc_info):
#     bc_info.psi_01 = dh_analyt_get_psi(0, bc_info)
#
#     c1 = dh_analyt_get_c1(bc_info)
#     c2 = dh_analyt_get_c2(bc_info)
#     bc_info.sigma_02_star = -1.0*c1*np.exp(bc_info.kap_h) - c2*np.exp(-1.0*bc_info.kap_h)
#     return
#
# def dh_analyt_get_psi_at_z_prime(z_prime, bc_info):
#     t = get_t_from_z_prime(z_prime, bc_info.kap_h)
#     return dh_analyt_get_psi(t, bc_info)
#
# def dh_get_integral(integrand_fun, bc_info):
#     return scipy.integrate.quad(integrand_fun, 0.0, bc_info.kap_h, args=(bc_info))
#
# def dh_sinh_integrand(z_prime, bc_info):
#     psi = dh_analyt_get_psi_at_z_prime(z_prime, bc_info)
#     # From Maclaurin series
#     return psi**2
#
# def dh_cosh_integrand(z_prime, bc_info):
#     psi = dh_analyt_get_psi_at_z_prime(z_prime, bc_info)
#     # From Maclaurin series
#     return (psi**2)/2.0
#
# def dh_get_scaled_energy(bc_info):
#     int_rel_tol   = 1.0e-3
#     sinh_integral = dh_get_integral(dh_sinh_integrand, bc_info)
#     cosh_integral = dh_get_integral(dh_cosh_integrand, bc_info)
#     assert (sinh_integral[1]/sinh_integral[0] < int_rel_tol)
#     assert (cosh_integral[1]/cosh_integral[0] < int_rel_tol)
#
#     u_el   = get_scaled_u_el(bc_info, sinh_integral[0])
#     ts     = get_scaled_ts(sinh_integral[0], cosh_integral[0])
#     f_chem = get_scaled_f_chem(bc_info)
#     return u_el + ts + f_chem
#
# def dh_get_bim_and_sl_pmf(ion_str, sigma_01, phi_02, kap_h):
#     sol        = bulk_solution(ion_str)
#     h          = kap_h/sol.kappa
#     bc_info    = cp_bc_info(sigma_01, phi_02, h, sol)
#     dh_update_pc_bc_info(bc_info)
#
#     pmf = dh_get_scaled_energy(bc_info)
#     bim = bc_info.sigma_01_star*bc_info.psi_01 + bc_info.sigma_02_star*bc_info.psi_02
#     return (pmf, bim)
#
#
# ################################################################################
# # Boneyard
#
# # Working correctly, just not needed
# def update_c1_for_coeff_vec(coeff_vec, bc_info):
#     # Used to obtain the initial guess for c1 based on the delta_i in coeff_vec
#     group = bc_info.kap_h*bc_info.sigma_01_star/4.0
#     t     = 0.0
#     b2    = get_b2(bc_info.c2, bc_info.kap_h, t)
#     b3    = get_b3(bc_info.kap_h, t)
#     nu    = get_nu(coeff_vec, t)
#
#     a = group
#     b = 1.0 + bc_info.kap_h
#     c = -1.0*group - b2 - b3*nu
#     option_1 = (-1.0*b + np.sqrt(b**2 - 4.0*a*c))/(2.0*a)
#     option_2 = (-1.0*b + np.sqrt(b**2 - 4.0*a*c))/(2.0*a)
#
#     tol = 1.0e-3
#     coeff_vec[-1] = option_1
#     if scaled_bc_residual(coeff_vec, bc_info) < tol:
#         pass
#     else:
#         coeff_vec[-1] = option_2
#         assert(scaled_bc_residual(coeff_vec, bc_info) < tol)
#     return coeff_vec
#
#
# def get_energy(coeff_vec, bc_info, area, bulk_solution):
#     f_plus = get_scaled_energy(coeff_vec, bc_info)
#     return unscale_energy(f_plus, area, bulk_solution)
#
# def get_psi_at_z(z, coeff_vec, bc_info):
#     t = get_t_from_z(z, bc_info.h)
#     return get_psi_at_t(t, coeff_vec, bc_info)
