import numpy as np
import matplotlib.pyplot as plt

import modules_ns_bim_spherical.orthogonal_collocation_v1 as oc

################################################################################

def align_yaxis(ax1, ax2):
    y_lims = np.array([ax.get_ylim() for ax in [ax1, ax2]])

    # force 0 to appear on both axes, comment if don't need
    y_lims[:, 0] = y_lims[:, 0].clip(None, 0)
    y_lims[:, 1] = y_lims[:, 1].clip(0, None)

    # normalize both axes
    y_mags = (y_lims[:,1] - y_lims[:,0]).reshape(len(y_lims),1)
    y_lims_normalized = y_lims / y_mags

    # find combined range
    y_new_lims_normalized = np.array([np.min(y_lims_normalized), np.max(y_lims_normalized)])

    # denormalize combined range to get new axes
    new_lim1, new_lim2 = y_new_lims_normalized * y_mags
    ax1.set_ylim(new_lim1)
    ax2.set_ylim(new_lim2)
    return


def plot_psi_and_residual(fit_1, bc_info):
    seps = np.linspace(0.0, 1, 2000)
    seps = seps.tolist()
    psi = []
    res = []

    for t in seps:
        psi.append(oc.get_psi_at_t(t, fit_1.x, bc_info))
        res.append(oc.scaled_1d_pb_residual(fit_1.x, t, bc_info))

    fig, ax = plt.subplots()
    font = {'weight' : 'normal',
            'size'   : 20}
    plt.rc('font', **font)
    fig.set_size_inches(6, 5, forward=True)
    ax2 = ax.twinx()
    ax.set_xlabel('t')
    ax.set_ylabel(r'$\psi$')
    ax2.set_ylabel('PB residual')
    ax.set_xlim(0.0, 1)
    ax.hlines(0, 0, 1)

    v_tick_len = 0.15
    p = len(fit_1.x) - 1
    for i in range(p):
        t_i = 1.0/(2.0*(p)) + i/(p)
        ax.axvline(x=t_i, ymin=0, ymax=0.05, color='k')

    lns1 = ax.plot(seps, psi, '-', label=r'$\psi$, {:.0f} $\delta_i$'.format(len(fit_1.x)-1))
    lns2 = ax2.plot(seps, res, '--', color='k', label='Residual')

    lns = lns1 + lns2
    labs = [l.get_label() for l in lns]
    ax.legend(lns, labs, loc='best', frameon=False, handlelength=1.0)

    align_yaxis(ax, ax2)
    return

def plot_energy_profile(hs, energy):
    fig, ax = plt.subplots()
    font    = {'weight' : 'normal', 'size'   : 20}
    plt.rc('font', **font)
    fig.set_size_inches(6, 5, forward=True)
    ax.set_ylabel('Energy [-]')
    ax.set_xlabel(r'Separation [$\kappa$H]')
    ax.set_xlim(0, 10)
    ax.hlines(0, 0, 100)
    # ax.set_ylim(-2, 1)

    ax.plot(hs, energy)
    # fig.savefig('oc_energy_profile.png', bbox_inches='tight', dpi = 300)
    return

def plot_psi_and_residual_cc(fit_1, bc_info):
    seps = np.linspace(0.0, 1, 2000)
    seps = seps.tolist()
    psi = []
    res = []

    for t in seps:
        psi.append(oc.get_psi_at_t_v3(t, fit_1, bc_info))
        res.append(oc.scaled_1d_pb_residual_cc(fit_1, t, bc_info))

    fig, ax = plt.subplots()
    font = {'weight' : 'normal',
            'size'   : 20}
    plt.rc('font', **font)
    fig.set_size_inches(6, 5, forward=True)
    ax2 = ax.twinx()
    ax.set_xlabel('t')
    ax.set_ylabel(r'$\psi$')
    ax2.set_ylabel('PB residual')
    ax.set_xlim(0.0, 1)
    ax.hlines(0, 0, 1)

    v_tick_len = 0.15
    p = len(fit_1) - 1
    for i in range(p):
        t_i = 1.0/(2.0*(p)) + i/(p)
        ax.axvline(x=t_i, ymin=0, ymax=0.05, color='k')

    lns1 = ax.plot(seps, psi, '-', label=r'$\psi$, {:.0f} $\delta_i$'.format(len(fit_1)-1))
    lns2 = ax2.plot(seps, res, '--', color='k', label='Residual')

    lns = lns1 + lns2
    labs = [l.get_label() for l in lns]
    ax.legend(lns, labs, loc='best', frameon=False, handlelength=1.0)

    align_yaxis(ax, ax2)
    return



















################################################################################
