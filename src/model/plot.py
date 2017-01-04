"""Plotting functions related to dalecv2.
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.mlab as mlab
import matplotlib.gridspec as gridspec
import mod_class as mc
import data_class as dc
import scipy
import pickle
import seaborn as sns
from matplotlib.patches import Rectangle

# ------------------------------------------------------------------------------
# Plot observation time series
# ------------------------------------------------------------------------------

def plot_phi(pvals, dC):
    """Plots phi_onset and phi_fall fns controlling leaf on and leaf off. Takes
    a parameter set (pvals) and a dataClass (dC).
    """
    sns.set_context(rc={'lines.linewidth':.8, 'lines.markersize':6})
    fig, ax = plt.subplots(nrows=1, ncols=1)
    m = mc.DalecModel(dC)
    mod_lst = m.mod_list(pvals)
    phi_on = m.oblist('phi_onset', mod_lst)
    phi_off = m.oblist('phi_fall', mod_lst)

    palette = sns.color_palette("colorblind", 11)

    ax.plot(dC.dates, phi_on, color=palette[0], label='phi_onset')
    ax.plot(dC.dates, phi_off, color=palette[2], label='phi_fall')
    plt.legend()
    ax.set_xlabel('Year')
    ax.set_ylabel('Rate of leaf on/leaf off')
    plt.gcf().autofmt_xdate()
    return ax, fig


def plot_gpp_sensitivity_lai():
    sns.set_context(rc={'lines.linewidth':.8, 'lines.markersize':6})
    fig, ax = plt.subplots(nrows=1, ncols=1)
    palette = sns.color_palette("colorblind", 11)
    d = dc.DalecDataTwin(1999, 2000, '')
    m= mc.DalecModel(d)
    m.x = 200
    gpp = []
    cf_list = np.arange(60., 600.)
    for cf in cf_list:
        gpp.append(m.acm(cf, 60., d.edinburgh_mean[10], d.acm))

    ax.plot(cf_list/60., gpp, color=palette[0])
    ax.set_xlabel('LAI')
    ax.set_ylabel('GPP (g C m-2 day-1)')
    return ax, fig


def plot_gpp_sensitivity_ceff():
    sns.set_context(rc={'lines.linewidth':.8, 'lines.markersize':6})
    fig, ax = plt.subplots(nrows=1, ncols=1)
    palette = sns.color_palette("colorblind", 11)
    d = dc.DalecDataTwin(1999, 2000, '')
    m= mc.DalecModel(d)
    m.x = 200
    gpp = []
    ceff_list = np.arange(10., 100.)
    for ceff in ceff_list:
        gpp.append(m.acm(180., 60., ceff, d.acm))

    ax.plot(ceff_list, gpp, color=palette[0])
    ax.set_xlabel('Canopy efficiency parameter')
    ax.set_ylabel('GPP (g C m-2 day-1)')
    return ax, fig


def plot_drive_dat(dat, dC):
    """Plots a specified observation using obs eqn in obs module. Takes an
    observation string, a dataClass (dC) and a start and finish point.
    """
    sns.set_context(rc={'lines.linewidth': 0.8, 'lines.markersize': 6})
    fig, ax = plt.subplots(nrows=1, ncols=1)
    drive_dict = {'t_mean': dC.t_mean, 't_max': dC.t_max, 't_min': dC.t_min,
                  'I': dC.I, 't_day': dC.t_day, 't_night': dC.t_night}
    palette = sns.color_palette("colorblind", 11)

    ax.plot(dC.dates, drive_dict[dat], color=palette[0])
    ax.set_xlabel('Year')
    ax.set_ylabel(dat)
    plt.gcf().autofmt_xdate()
    return ax, fig


def plot_ob_dict(ob, dC):
    """Plots a specified observation using obs eqn in obs module. Takes an
    observation string, a dataClass (dC) and a start and finish point.
    """
    sns.set_context(rc={'lines.linewidth': 0.8, 'lines.markersize': 6})
    fig, ax = plt.subplots(nrows=1, ncols=1)

    palette = sns.color_palette("colorblind", 11)

    ax.plot(dC.dates, dC.ob_dict[ob], 'o', color=palette[0])
    ax.set_xlabel('Year')
    ax.set_ylabel(ob)
    plt.gcf().autofmt_xdate()
    return ax, fig


def plot_ob_dict_east_west(ob, dC_east, dC_west):
    """Plots a specified observation using obs eqn in obs module. Takes an
    observation string, a dataClass (dC) and a start and finish point.
    """
    sns.set_context(rc={'lines.linewidth': 0.8, 'lines.markersize': 6})
    fig, ax = plt.subplots(nrows=1, ncols=1)

    palette = sns.color_palette("colorblind", 11)

    ax.plot(dC_east.dates, dC_east.ob_dict[ob], 'o', color=palette[0], label='East', markeredgecolor='black',
            markeredgewidth=0.5)
    ax.plot(dC_west.dates, dC_west.ob_dict[ob], 'o', color=palette[2], label='West', markeredgecolor='black',
            markeredgewidth=0.5)
    plt.legend()
    ax.set_xlabel('Year')
    ax.set_ylabel(ob)
    plt.gcf().autofmt_xdate()
    return ax, fig


def plot_obs(ob, pvals, dC):
    """Plots a specified observation using obs eqn in obs module. Takes an
    observation string, a dataClass (dC) and a start and finish point.
    """
    sns.set_context(rc={'lines.linewidth': 0.8, 'lines.markersize': 6})
    fig, ax = plt.subplots(nrows=1, ncols=1)
    m = mc.DalecModel(dC)
    mod_lst = m.mod_list(pvals)
    obs_lst = m.oblist(ob, mod_lst)

    palette = sns.color_palette("colorblind", 11)

    ax.plot(dC.dates, obs_lst, color=palette[0])
    ax.set_xlabel('Year')
    ax.set_ylabel(ob)
    plt.gcf().autofmt_xdate()
    return ax, fig


def plot_obs_east_west(ob, xa_east, xa_west, d_e, d_w, y_label='None', xb='None', ob_std_e=0, ob_std_w=0, y_lim='None',
                       axes='None'):
    """Plots a specified observation using obs eqn in obs module. Takes an
    observation string, a dataClass (dC) and a start and finish point.
    """
    sns.set_context(rc={'lines.linewidth': .8, 'lines.markersize': 6})
    if axes == 'None':
        fig, ax = plt.subplots(nrows=1, ncols=1)
        ret_val = ax, fig
    else:
        ax = axes
        ret_val = ax
    me = mc.DalecModel(d_e)
    mw = mc.DalecModel(d_w)
    mod_lst_e = me.mod_list(xa_east)
    obs_lst_e = me.oblist(ob, mod_lst_e)
    mod_lst_w = mw.mod_list(xa_west)
    obs_lst_w = mw.oblist(ob, mod_lst_w)

    palette = sns.color_palette("colorblind", 11)

    ax.plot(d_e.dates, obs_lst_e, color=palette[0], label='Unthinned')
    if ob_std_e is not 0:
        ax.fill_between(d_e.dates, obs_lst_e-ob_std_e, obs_lst_e+ob_std_e, facecolor=palette[0],
                        alpha=0.5, linewidth=0.0)
    ax.plot(d_w.dates, obs_lst_w, color=palette[2], label='Thinned')
    if ob_std_w is not 0:
        ax.fill_between(d_w.dates, obs_lst_w-ob_std_e, obs_lst_w+ob_std_e, facecolor=palette[2],
                        alpha=0.5, linewidth=0.0)
    if xb != 'None':
        mod_lst_xb = mw.mod_list(xb)
        obs_lst_xb = mw.oblist(ob, mod_lst_xb)
        ax.plot(d_w.dates, obs_lst_xb, '--', color=palette[3], label='Prior model')
    if ob in d_e.ob_dict.keys():
        ax.errorbar(d_e.dates, d_e.ob_dict[ob], yerr=d_e.ob_err_dict[ob], fmt='o', color=palette[0],
                    markeredgecolor='black', markeredgewidth=0.5)
        ax.errorbar(d_w.dates, d_w.ob_dict[ob], yerr=d_w.ob_err_dict[ob], fmt='o', color=palette[2],
                    markeredgecolor='black', markeredgewidth=0.5)
    plt.legend()
    ax.set_xlabel('Date')
    if y_label == 'None':
        ax.set_ylabel(ob)
    else:
        ax.set_ylabel(y_label)
    if y_lim != 'None':
        axes = plt.gca()
        axes.set_ylim(y_lim)
    plt.gcf().autofmt_xdate()
    return ret_val


def plot_prior(ob, xa, d, y_label='None', y_lim='None', axes='None'):
    """Plots a specified observation using obs eqn in obs module. Takes an
    observation string, a dataClass (dC) and a start and finish point.
    """
    sns.set_context(rc={'lines.linewidth': .8, 'lines.markersize': 6})
    if axes == 'None':
        fig, ax = plt.subplots(nrows=1, ncols=1)
        ret_val = ax, fig
    else:
        ax = axes
        ret_val = ax
    m = mc.DalecModel(d)
    mod_lst = m.mod_list(xa)
    obs_lst = m.oblist(ob, mod_lst)


    palette = sns.color_palette("colorblind", 11)

    ax.plot(d.dates, obs_lst, color=palette[1])


    if ob in d.ob_dict.keys():
        ax.errorbar(d.dates, d.ob_dict[ob], yerr=d.ob_err_dict[ob], fmt='o', color=palette[2],
                    markeredgecolor='black', markeredgewidth=0.5)
    ax.set_xlabel('Date')
    if y_label == 'None':
        ax.set_ylabel(ob)
    else:
        ax.set_ylabel(y_label)
    if y_lim != 'None':
        axes = plt.gca()
        axes.set_ylim(y_lim)
    plt.gcf().autofmt_xdate()
    return ret_val


def plot_obs_east_west_cum(ob, xa_east, xa_west, d_e, d_w, y_label='None', xb='None', ob_std_e=0, ob_std_w=0,
                           y_lim='None'):
    """Plots a specified observation using obs eqn in obs module. Takes an
    observation string, a dataClass (dC) and a start and finish point.
    """
    sns.set_context(rc={'lines.linewidth': .8, 'lines.markersize': 6})
    fig, ax = plt.subplots(nrows=1, ncols=1)
    me = mc.DalecModel(d_e)
    mw = mc.DalecModel(d_w)
    mod_lst_e = me.mod_list(xa_east)
    obs_lst_e = me.oblist(ob, mod_lst_e)
    mod_lst_w = mw.mod_list(xa_west)
    obs_lst_w = mw.oblist(ob, mod_lst_w)

    palette = sns.color_palette("colorblind", 11)

    ax.plot(d_e.dates, np.cumsum(obs_lst_e), color=palette[0], label='Unthinned')
    if ob_std_e is not 0:
        ax.fill_between(d_e.dates, np.cumsum(obs_lst_e)-ob_std_e, np.cumsum(obs_lst_e)+ob_std_e, facecolor=palette[0],
                        alpha=0.5, linewidth=0.0)
        cum_east = y_label+'_east: ' + str(np.cumsum(obs_lst_e)[-1])+' +/- ' + str(ob_std_e[-1])
        print cum_east
    ax.plot(d_w.dates, np.cumsum(obs_lst_w), color=palette[2], label='Thinned')
    if ob_std_w is not 0:
        ax.fill_between(d_w.dates, np.cumsum(obs_lst_w)-ob_std_w, np.cumsum(obs_lst_w)+ob_std_w, facecolor=palette[2],
                        alpha=0.5, linewidth=0.0)
        cum_west = y_label+'_west: ' + str(np.cumsum(obs_lst_w)[-1])+' +/- ' + str(ob_std_w[-1])
        print cum_west
    if xb != 'None':
        mod_lst_xb = mw.mod_list(xb)
        obs_lst_xb = mw.oblist(ob, mod_lst_xb)
        ax.plot(d_w.dates, obs_lst_xb, '--', color=palette[3], label='Prior model')
    plt.legend()
    ax.set_xlabel('Date')
    if y_label == 'None':
        ax.set_ylabel(ob)
    else:
        ax.set_ylabel(y_label)
    if y_lim != 'None':
        axes = plt.gca()
        axes.set_ylim(y_lim)
    plt.gcf().autofmt_xdate()
    if ob_std_e is not 0:
        return ax, fig, cum_east, cum_west
    else:
        return ax, fig



def plot_obs_cum(ob, xa, d, y_label='None', ob_std=0, y_lim='None'):
    """Plots a specified observation using obs eqn in obs module. Takes an
    observation string, a dataClass (dC) and a start and finish point.
    """
    sns.set_context(rc={'lines.linewidth': .8, 'lines.markersize': 6})
    fig, ax = plt.subplots(nrows=1, ncols=1)
    m = mc.DalecModel(d)
    mod_lst = m.mod_list(xa)
    obs_lst = m.oblist(ob, mod_lst)

    palette = sns.color_palette("colorblind", 11)

    ax.plot(d.dates, np.cumsum(obs_lst), color=palette[0], label='Unthinned')
    if ob_std is not 0:
        ax.fill_between(d.dates, np.cumsum(obs_lst)-ob_std, np.cumsum(obs_lst)+ob_std, facecolor=palette[0],
                        alpha=0.5, linewidth=0.0)
        cum_east = y_label+'_east: ' + str(np.cumsum(obs_lst)[-1])+' +/- ' + str(ob_std[-1])
        print cum_east
    #plt.legend()
    ax.set_xlabel('Date')
    if y_label == 'None':
        ax.set_ylabel(ob)
    else:
        ax.set_ylabel(y_label)
    if y_lim != 'None':
        axes = plt.gca()
        axes.set_ylim(y_lim)
    plt.gcf().autofmt_xdate()
    if ob_std is not 0:
        return ax, fig, cum_east
    else:
        return ax, fig



def plot_obs_east_west_cum_part(xa_east, xa_west, d_e, d_w, ob_std_e, ob_std_w, y_label='None', xb='None',
                                y_lim='None', axes='None'):
    """Plots a specified observation using obs eqn in obs module. Takes an
    observation string, a dataClass (dC) and a start and finish point.
    """
    sns.set_context(rc={'lines.linewidth': .8, 'lines.markersize': 6})
    if axes == 'None':
        fig, ax = plt.subplots(nrows=1, ncols=1)
        ret_val = ax, fig
    else:
        ax = axes
        ret_val = ax
    me = mc.DalecModel(d_e)
    mw = mc.DalecModel(d_w)
    mod_lst_e = me.mod_list(xa_east)
    nee_lst_e = me.oblist('nee', mod_lst_e)
    gpp_lst_e = me.oblist('gpp', mod_lst_e)
    rt_lst_e = me.oblist('rt', mod_lst_e)
    mod_lst_w = mw.mod_list(xa_west)
    nee_lst_w = mw.oblist('nee', mod_lst_w)
    gpp_lst_w = me.oblist('gpp', mod_lst_w)
    rt_lst_w = me.oblist('rt', mod_lst_w)

    palette = sns.color_palette("colorblind", 11)

    ax.plot(d_e.dates, np.cumsum(nee_lst_e), color=palette[0], label='Unthinned')
    ax.fill_between(d_e.dates, np.cumsum(nee_lst_e)-ob_std_e[0], np.cumsum(nee_lst_e)+ob_std_e[0], facecolor=palette[0],
                        alpha=0.5, linewidth=0.0)
    ax.plot(d_e.dates, -np.cumsum(gpp_lst_e), '--', color=palette[0])
    ax.fill_between(d_e.dates, -np.cumsum(gpp_lst_e)-ob_std_e[1], -np.cumsum(gpp_lst_e)+ob_std_e[1], facecolor=palette[0],
                        alpha=0.5, linewidth=0.0)
    ax.plot(d_e.dates, np.cumsum(rt_lst_e), ':', color=palette[0])
    ax.fill_between(d_e.dates, np.cumsum(rt_lst_e)-ob_std_e[2], np.cumsum(rt_lst_e)+ob_std_e[2], facecolor=palette[0],
                        alpha=0.5, linewidth=0.0)

    ax.plot(d_w.dates, np.cumsum(nee_lst_w), color=palette[2], label='Thinned')
    ax.fill_between(d_w.dates, np.cumsum(nee_lst_w)-ob_std_w[0], np.cumsum(nee_lst_w)+ob_std_w[0], facecolor=palette[2],
                        alpha=0.5, linewidth=0.0)
    ax.plot(d_e.dates, -np.cumsum(gpp_lst_w), '--', color=palette[2])
    ax.fill_between(d_e.dates, -np.cumsum(gpp_lst_w)-ob_std_w[1], -np.cumsum(gpp_lst_w)+ob_std_w[1], facecolor=palette[2],
                        alpha=0.5, linewidth=0.0)
    ax.plot(d_e.dates, np.cumsum(rt_lst_w), ':', color=palette[2])
    ax.fill_between(d_e.dates, np.cumsum(rt_lst_w)-ob_std_w[2], np.cumsum(rt_lst_w)+ob_std_w[2], facecolor=palette[2],
                        alpha=0.5, linewidth=0.0)

    plt.legend(loc=2)
    ax.set_xlabel('Date')
    if y_label == 'None':
        ax.set_ylabel('NEE partitioning')
    else:
        ax.set_ylabel(y_label)
    if y_lim != 'None':
        axes = plt.gca()
        axes.set_ylim(y_lim)
    plt.gcf().autofmt_xdate()
    return ret_val


def plot_obs_east_west_part(xa, d, y_lim='None', axes='None'):
    """Plots a specified observation using obs eqn in obs module. Takes an
    observation string, a dataClass (dC) and a start and finish point.
    """
    sns.set_context(rc={'lines.linewidth': .8, 'lines.markersize': 6})
    if axes == 'None':
        fig, ax = plt.subplots(nrows=1, ncols=1)
        ret_val = ax, fig
    else:
        ax = axes
        ret_val = ax
    m = mc.DalecModel(d)
    mod_lst = m.mod_list(xa)
    soilr_lst = m.oblist('soilresp', mod_lst)
    litr_lst = m.oblist('litresp', mod_lst)
    autor_lst = m.oblist('ra', mod_lst)
    gpp_lst = m.oblist('gpp', mod_lst)
    rt_lst = m.oblist('rt', mod_lst)

    palette = sns.color_palette("colorblind", 11)

    ax.stackplot(d.dates, np.cumsum(soilr_lst), np.cumsum(litr_lst), np.cumsum(autor_lst),
                 colors=(palette[0], palette[1], palette[2]),)
    print np.cumsum(soilr_lst)[-1], np.cumsum(litr_lst)[-1], np.cumsum(autor_lst)[-1]
    #ax.text(d.dates[-90], 50, 'Soil respiration', color='w')
    #ax.text(d.dates[-100], 350, 'Litter respiration', color='w')
    #ax.text(d.dates[-120], 940, 'Autotrohpic respiration', color='w')
    #ax.stackplot(d.dates, -np.cumsum(gpp_lst), color=palette[4])
    p1 = Rectangle((0, 0), 1, 1, fc=palette[0])
    p2 = Rectangle((0, 0), 1, 1, fc=palette[1])
    p3 = Rectangle((0, 0), 1, 1, fc=palette[2])
    plt.legend([p3, p2, p1], ['Autotrophic respiration', 'Litter respiration', 'Soil respiration'], loc=2)
    ax.set_xlabel('Date')
    ax.set_ylabel(r'Cumulative respiration partitioning (g C m$^{-2}$)')
    ax.set_ylim([0, 1800])
    plt.gcf().autofmt_xdate()
    return ret_val


def part_plot(xa_east, xa_west, d_e, d_w):
    sns.set_context('poster', font_scale=1., rc={'lines.linewidth': .8, 'lines.markersize': 1.})
    sns.set_style('whitegrid')
    f, ((ax1, ax2)) = plt.subplots(1, 2)
    # Observed values
    plot_obs_east_west_part(xa_east, d_e, axes=ax1)
    ax1.set_title(r'a) Unthinned forest')#, y=1.06)

    plot_obs_east_west_part(xa_west, d_w, axes=ax2)
    ax2.set_title(r'b) Thinned forest')#, y=1.06)

    f.tight_layout()
    #f.subplots_adjust(hspace=.5)
    return f


def plot_pheno_obs(d, axes='None'):
    """Plots a specified observation using obs eqn in obs module. Takes an
    observation string, a dataClass (dC) and a start and finish point.
    """
    sns.set_context(rc={'lines.linewidth': .8, 'lines.markersize': 6})
    pheno = mlab.csv2rec('/Users/ewan/projects/ah_data/pheno_obs_15.csv')
    if axes == 'None':
        fig, ax = plt.subplots(nrows=1, ncols=1)
        ret_val = ax, fig
    else:
        ax = axes
        ret_val = ax
    palette = sns.color_palette("colorblind", 11)

    ax.plot(d.dates, pheno['canopy_roi'][0:-1], 'o', color=palette[1])

    ax.set_xlabel('Date')
    ax.set_ylabel(r'Green fraction')
    #ax.set_ylim([0.34, 0.4])
    plt.gcf().autofmt_xdate()
    return ret_val


def plot_4dvar(ob, dC, xb=None, xa=None, erbars=1, y_label='None', awindl=None, obdict_a=None):
    """Plots a model predicted observation value for two initial states (xb,xa)
    and also the actual observations taken of the physical quantity. Takes a ob
    string, two initial states (xb,xa), a dataClass and a start and finish
    time step.
    """
    sns.set_context(rc={'lines.linewidth':.8, 'lines.markersize':6})
    fig, ax = plt.subplots(nrows=1, ncols=1)
    m = mc.DalecModel(dC)
    palette = sns.color_palette("colorblind", 11)
    if xb != None:
        mod_lst = m.mod_list(xb)
        obs_lst = m.oblist(ob, mod_lst)
        ax.plot(dC.dates, obs_lst, color=palette[0])
    if xa != None:
        mod_lst = m.mod_list(xa)
        obs_lst = m.oblist(ob, mod_lst)
        ax.plot(dC.dates, obs_lst, color=palette[1])

    ob_dict = dC.ob_dict
    ob_err_dict = dC.ob_err_dict
    if ob in ob_dict.keys():
        if erbars == True:
            ax.errorbar(dC.dates, ob_dict[ob], yerr=ob_err_dict[ob],
                        fmt='o', label=ob+'_o', color=palette[2], alpha=0.7)
        else:
            ax.plot(dC.dates, ob_dict[ob], 'o', label=ob+'_o', color=palette[2])
    if obdict_a != None:
        ax.plt.plot(dC.dates, obdict_a[ob], 'o')

    if awindl != None:
        ax.axvline(x=dC.dates[awindl], color='k', ls='dashed')

    ax.set_xlabel('Year')
    if y_label == 'None':
        ax.set_ylabel(ob)
    else:
        ax.set_ylabel(y_label)
    plt.gcf().autofmt_xdate()
    return ax, fig


def plot_scatter(ob, pvals, dC, awindl, bfa='a'):
    """Plots scatter plot of obs vs model predicted values. Takes an initial
    parameter set, a dataClass (must have only desired ob for comparison
    specified in dC), assimilation window length and whether a comparison of
    background 'b', forecast 'f' or analysis 'a' is desired.
    """
    sns.set_context('poster', font_scale=1.5, rc={'lines.linewidth': 1., 'lines.markersize': 6.})
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 10))
    sns.set_style('ticks')
    palette = sns.color_palette("colorblind", 11)
    m = mc.DalecModel(dC)
    mod_lst = m.mod_list(pvals)
    obs_lst = m.oblist(ob, mod_lst)
    y_obs = dC.ob_dict[ob]
    plt_ob_lst = (y_obs/y_obs)*obs_lst
    if bfa == 'b' or bfa == 'a':
        selection = xrange(0, awindl)
    elif bfa == 'f':
        selection = xrange(awindl, len(obs_lst))
    else:
        raise Exception('Please check function input for bfa variable')
    ob_lst = plt_ob_lst[selection][np.isnan(y_obs[selection]) != True]
    y_obs = y_obs[selection][np.isnan(y_obs[selection]) != True]
    if ob in ['lai', 'c_woo']:
        ob_lst = ob_lst.compressed()
        y_obs = y_obs.compressed()

    one_one = np.arange((min(min(y_obs), min(ob_lst))),(max(max(y_obs), max(ob_lst))))
    plt.plot(one_one, one_one, color=palette[0])

    ax.plot(y_obs, ob_lst, 'o', color=palette[1])
    error = np.sqrt(np.sum((y_obs - ob_lst)**2) / len(y_obs))
    yhx = np.mean(y_obs - ob_lst)
    mod_obs_bar = np.mean(ob_lst)
    std_mod_obs = np.nanstd(ob_lst)
    obs_bar = np.mean(y_obs)
    std_obs = np.std(y_obs)
    sse = np.sum([(y_obs[x]-ob_lst[x])**2 for x in xrange(len(y_obs))])
    sst = np.sum([(y_obs[x]-obs_bar)**2 for x in xrange(len(y_obs))])
    r2 = 1 - (sse/sst)
    rms = np.sqrt(np.sum([((ob_lst[x]-mod_obs_bar)-(y_obs[x]-obs_bar))**2 for x in range(len(y_obs))]) / len(y_obs))
    corr_coef = (np.sum([((ob_lst[x]-mod_obs_bar)*(y_obs[x]-obs_bar)) for x in range(len(y_obs))]) / len(y_obs)) / \
                (std_mod_obs*std_obs)

    plt.xlabel(ob.upper()+r' observations (g C m$^{-2}$ day$^{-1}$)')
    plt.ylabel(ob.upper()+' model (g C m$^{-2}$ day$^{-1}$)')
    plt.title('mean(y-hx)=%.2f, rms=%.2f, corr_coef=%.2f' %( yhx, rms, corr_coef))
    print bfa+'_error=%f, mean(y-hx)=%f, rms=%f, corr_coef=%f, r2=%f' %(error, yhx, rms, corr_coef, r2)
    #plt.xlim((-20, 15))
    #plt.ylim((-20, 15))
    return ax, fig,


def plot_a_inc(xb, xa, lab='a_inc'):
    """Plot error between truth and xa/xb shows as a bar chart.
    """
    sns.set_context('poster', font_scale=1.5, rc={'lines.linewidth':1, 'lines.markersize':10})
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(15,5))
    sns.set_style('ticks')

    n = 23
    width = 0.5
    ind = np.arange(n)
    #fig = plt.figure()
    #ax = fig.add_subplot(111)
    rects1 = ax.bar(ind, (xa-xb)/xb, width, color=sns.xkcd_rgb["faded green"],
                    label=lab)
    ax.set_ylabel('Normalised analysis increment')
    #ax.set_title('% error in parameter values for xa and xb')
    ax.set_xticks(ind+width/2)
    keys = [r'$\theta_{min}$', r'$f_{auto}$', r'$f_{fol}$', r'$f_{roo}$', r'$c_{lspan}$', r'$\theta_{woo}$',
            r'$\theta_{roo}$', r'$\theta_{lit}$', r'$\theta_{som}$', r'$\Theta$', r'$c_{eff}$', r'$d_{onset}$',
            r'$f_{lab}$', r'$c_{ronset}$', r'$d_{fall}$', r'$c_{rfall}$', r'$c_{lma}$', r'$C_{lab}$', r'$C_{fol}$',
            r'$C_{roo}$', r'$C_{woo}$', r'$C_{lit}$', r'$C_{som}$']
    ax.set_xticklabels(keys, rotation=90)
    ax.legend()
    return ax, fig


def plot_inc_east_west(xb, xa_east, xa_west):
    """Plot error between truth and xa/xb shows as a bar chart.
    """
    sns.set_context('poster', font_scale=1.5, rc={'lines.linewidth': 1, 'lines.markersize': 10})
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(15, 5))
    sns.set_style('ticks')
    palette = sns.color_palette("colorblind", 11)
    n = 23
    width = 0.35
    ind = np.arange(n)
    #fig = plt.figure()
    #ax = fig.add_subplot(111)
    rects1 = ax.bar(ind, 100*(xa_east-xb)/xb, width, color=sns.xkcd_rgb["faded green"],
                    label='Unthinned')
    rects2 = ax.bar(ind+width, 100*(xa_west-xb)/xb, width, color=sns.xkcd_rgb["pale red"],
                    label='Thinned')
    ax.set_ylabel('Normalised analysis increment (%)')
    #ax.set_title('% error in parameter values for xa and xb')
    ax.set_xticks(ind+width)
    keys = [r'$\theta_{min}$', r'$f_{auto}$', r'$f_{fol}$', r'$f_{roo}$', r'$c_{lspan}$', r'$\theta_{woo}$',
            r'$\theta_{roo}$', r'$\theta_{lit}$', r'$\theta_{som}$', r'$\Theta$', r'$c_{eff}$', r'$d_{onset}$',
            r'$f_{lab}$', r'$c_{ronset}$', r'$d_{fall}$', r'$c_{rfall}$', r'$c_{lma}$', r'$C_{lab}$', r'$C_{fol}$',
            r'$C_{roo}$', r'$C_{woo}$', r'$C_{lit}$', r'$C_{som}$']
    ax.set_xticklabels(keys, rotation=90)
    ax.legend(loc=2)
    return ax, fig


def plot_var_red_east_west(b_mat, a_cov_east, a_cov_west):
    """Plot error between truth and xa/xb shows as a bar chart.
    """
    sns.set_context('poster', font_scale=1.5, rc={'lines.linewidth': 1, 'lines.markersize': 10})
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(15, 5))
    sns.set_style('ticks')
    palette = sns.color_palette("colorblind", 11)
    n = 23
    width = 0.35
    ind = np.arange(n)
    #fig = plt.figure()
    #ax = fig.add_subplot(111)
    xa_east = np.sqrt(np.diag(a_cov_east))
    xa_west = np.sqrt(np.diag(a_cov_west))
    xb = np.sqrt(np.diag(b_mat))
    rects1 = ax.bar(ind, -100*(xa_east-xb)/xb, width, color=sns.xkcd_rgb["faded green"],
                    label='Unthinned')
    rects2 = ax.bar(ind+width, -100*(xa_west-xb)/xb, width, color=sns.xkcd_rgb["pale red"],
                    label='Thinned')
    ax.set_ylabel('Reduction in error (%)')
    #ax.set_title('% error in parameter values for xa and xb')
    ax.set_xticks(ind+width)
    keys = [r'$\theta_{min}$', r'$f_{auto}$', r'$f_{fol}$', r'$f_{roo}$', r'$c_{lspan}$', r'$\theta_{woo}$',
            r'$\theta_{roo}$', r'$\theta_{lit}$', r'$\theta_{som}$', r'$\Theta$', r'$c_{eff}$', r'$d_{onset}$',
            r'$f_{lab}$', r'$c_{ronset}$', r'$d_{fall}$', r'$c_{rfall}$', r'$c_{lma}$', r'$C_{lab}$', r'$C_{fol}$',
            r'$C_{roo}$', r'$C_{woo}$', r'$C_{lit}$', r'$C_{som}$']
    ax.set_xticklabels(keys, rotation=90)
    ax.legend(loc=2)
    return ax, fig


def plot_scatter_twin(ob, pvals, dC, awindl, bfa='a'):
    """Plots scatter plot of obs vs model predicted values. Takes an initial
    parameter set, a dataClass (must have only desired ob for comparison
    specified in dC), assimilation window length and whether a comparison of
    background 'b', forecast 'f' or analysis 'a' is desired.
    """
    sns.set_context('poster', font_scale=1.5, rc={'lines.linewidth': 1., 'lines.markersize': 6.})
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 10))
    #sns.set_style('ticks')
    palette = sns.color_palette("colorblind", 11)
    m = mc.DalecModel(dC)
    mod_lst = m.mod_list(pvals)
    mod_lst_truth = m.mod_list(dC.x_truth)
    obs_lst = m.oblist(ob, mod_lst)
    y_obs = m.oblist(ob, mod_lst_truth)
    plt_ob_lst = (y_obs/y_obs)*obs_lst
    if bfa == 'b' or bfa == 'a':
        selection = xrange(0, awindl)
    elif bfa == 'f':
        selection = xrange(awindl, len(obs_lst))
    else:
        raise Exception('Please check function input for bfa variable')
    ob_lst = plt_ob_lst[selection][np.isnan(y_obs[selection]) != True]
    y_obs = y_obs[selection][np.isnan(y_obs[selection]) != True]

    one_one = np.arange(int(min(min(y_obs), min(ob_lst))), int(max(max(y_obs), max(ob_lst))))
    plt.plot(one_one, one_one, color=palette[0])
    print int(min(min(y_obs), min(ob_lst))), int(max(max(y_obs), max(ob_lst)))

    ax.plot(y_obs, ob_lst, 'o', color=palette[1])
    error = np.sqrt(np.sum((y_obs - ob_lst)**2) / len(y_obs))
    yhx = np.mean(y_obs - ob_lst)
    mod_obs_bar = np.mean(ob_lst)
    std_mod_obs = np.nanstd(ob_lst)
    obs_bar = np.mean(y_obs)
    std_obs = np.std(y_obs)
    rms = np.sqrt(np.sum([((ob_lst[x]-mod_obs_bar)-(y_obs[x]-obs_bar))**2 for x in range(len(y_obs))]) / len(y_obs))
    corr_coef = (np.sum([((ob_lst[x]-mod_obs_bar)*(y_obs[x]-obs_bar)) for x in range(len(y_obs))]) / len(y_obs)) / \
                (std_mod_obs*std_obs)

    plt.xlabel(ob.upper()+r' observations (g C m$^{-2}$ day$^{-1}$)')
    plt.ylabel(ob.upper()+' model (g C m$^{-2}$ day$^{-1}$)')
    plt.title('mean(y-hx)=%.2f, rms=%.2f, corr_coef=%.2f' %( yhx, rms, corr_coef))
    print bfa+'_error=%f, mean(y-hx)=%f, rms=%f, corr_coef=%f' %(error, yhx, rms, corr_coef)
    #plt.xlim((-20, 15))
    #plt.ylim((-20, 15))
    return ax, fig


def plot_4dvar_twin(ob, dC, xb=None, xa=None, erbars=1, awindl=None, obdict_a=None):
    """Plots a model predicted observation value for two initial states (xb,xa)
    and also the actual observations taken of the physical quantity. Takes a ob
    string, two initial states (xb,xa), a dataClass and a start and finish
    time step.
    """
    sns.set_context(rc={'lines.linewidth':.8, 'lines.markersize':6})
    fig, ax = plt.subplots(nrows=1, ncols=1)
    palette = sns.color_palette("colorblind", 11)
    m = mc.DalecModel(dC)
    mod_lst = m.mod_list(dC.x_truth)
    obs_lst = m.oblist(ob, mod_lst)
    ax.plot(dC.dates, obs_lst, color=palette[3])
    if xb != None:
        mod_lst = m.mod_list(xb)
        obs_lst = m.oblist(ob, mod_lst)
        ax.plot(dC.dates, obs_lst, ':', color=palette[0])
    if xa != None:
        mod_lst = m.mod_list(xa)
        obs_lst = m.oblist(ob, mod_lst)
        ax.plot(dC.dates, obs_lst, color=palette[1])

    mod_lst = m.mod_list(dC.x_truth)
    obs_lst = m.oblist(ob, mod_lst)
    ax.plot(dC.dates, obs_lst, '--', color=palette[3])

    # ob_dict = obdict_a
    # ob_err_dict = dC.ob_err_dict
    # if ob in ob_dict.keys():
    #    if erbars == True:
    #        ax.errorbar(dC.dates, ob_dict[ob], yerr=ob_err_dict[ob],
    #                     fmt='o', label=ob+'_o', color=palette[2], alpha=0.7)
    #    else:
    #        ax.plot(dC.dates, ob_dict[ob], 'o', label=ob+'_o', color=palette[2])
    if obdict_a != None:
        ax.plot(dC.dates[0:len(obdict_a[ob])], obdict_a[ob], 'o', color=palette[2])

    if awindl != None:
        ax.axvline(x=dC.dates[awindl], color='k', ls='dashed')

    ax.set_xlabel('Year')
    ax.set_ylabel(ob)
    plt.gcf().autofmt_xdate()

    return ax, fig


def plottwinerr(truth, xb, xa, xb_lab='xb', xa_lab='xa'):
    """Plot error between truth and xa/xb shows as a bar chart.
    """
    sns.set_context('poster', font_scale=1.5, rc={'lines.linewidth': 1, 'lines.markersize': 10})
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(15,5))
    #sns.set_style('ticks')
    n = 23
    width = 0.35
    ind = np.arange(n)
    rects1 = ax.bar(ind, 100*abs(truth-xb)/truth, width, color=sns.xkcd_rgb["faded green"], label=xb_lab+'_err')
    rects2 = ax.bar(ind+width, 100*abs(truth-xa)/truth, width, color=sns.xkcd_rgb["pale red"], label=xa_lab+'_err')
    ax.set_ylabel('% error')
    ax.set_title('% error in parameter values for '+xa_lab+' and '+xb_lab)
    ax.set_xticks(ind+width)
    keys = [r'$\theta_{min}$', r'$f_{auto}$', r'$f_{fol}$', r'$f_{roo}$', r'$c_{lspan}$', r'$\theta_{woo}$',
            r'$\theta_{roo}$', r'$\theta_{lit}$', r'$\theta_{som}$', r'$\Theta$', r'$c_{eff}$', r'$d_{onset}$',
            r'$f_{lab}$', r'$c_{ronset}$', r'$d_{fall}$', r'$c_{rfall}$', r'$c_{lma}$', r'$C_{lab}$', r'$C_{fol}$',
            r'$C_{roo}$', r'$C_{woo}$', r'$C_{lit}$', r'$C_{som}$']
    ax.set_xticklabels(keys, rotation=90)
    ax.legend()
    return ax, fig


def plot_a_inc_all(xb, xadiag, xaedc, xarcor, xaedcrcor):
    """Plot error between truth and xa/xb shows as a bar chart.
    """
    sns.set_context('poster', font_scale=1.5, rc={'lines.linewidth':1, 'lines.markersize':10})
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(15,5))
    #sns.set_style('ticks')
    n = 23
    width = 0.22
    ind = np.arange(n)
    #fig = plt.figure()
    #ax = fig.add_subplot(111)
    rects1 = ax.bar(ind, (xadiag-xb)/xb, width, color=sns.xkcd_rgb["faded green"],
                    label='A')
    rects2 = ax.bar(ind+width, (xaedc-xb)/xb, width, color=sns.xkcd_rgb["pale red"],
                    label='B')
    rects3 = ax.bar(ind+width*2, (xarcor-xb)/xb, width, color=sns.xkcd_rgb["dusty purple"],
                    label='C')
    rects4 = ax.bar(ind+width*3, (xaedcrcor-xb)/xb, width, color=sns.xkcd_rgb["amber"],
                    label='D')
    ax.set_ylabel('Normalised analysis increment')
    #ax.set_title('% error in parameter values for xa and xb')
    ax.set_xticks(ind+width*2)
    keys = [r'$\theta_{min}$', r'$f_{auto}$', r'$f_{fol}$', r'$f_{roo}$', r'$c_{lspan}$', r'$\theta_{woo}$',
            r'$\theta_{roo}$', r'$\theta_{lit}$', r'$\theta_{som}$', r'$\Theta$', r'$c_{eff}$', r'$d_{onset}$',
            r'$f_{lab}$', r'$c_{ronset}$', r'$d_{fall}$', r'$c_{rfall}$', r'$c_{lma}$', r'$C_{lab}$', r'$C_{fol}$',
            r'$C_{roo}$', r'$C_{woo}$', r'$C_{lit}$', r'$C_{som}$']
    ax.set_xticklabels(keys, rotation=90)
    ax.legend()
    return ax, fig


def plot_table(xb, xa_east, xa_west):
    keys = [r'$\theta_{min}$', r'$f_{auto}$', r'$f_{fol}$', r'$f_{roo}$', r'$c_{lspan}$', r'$\theta_{woo}$',
            r'$\theta_{roo}$', r'$\theta_{lit}$', r'$\theta_{som}$', r'$\Theta$', r'$c_{eff}$', r'$d_{onset}$',
            r'$f_{lab}$', r'$c_{ronset}$', r'$d_{fall}$', r'$c_{rfall}$', r'$c_{lma}$', r'$C_{lab}$', r'$C_{fol}$',
            r'$C_{roo}$', r'$C_{woo}$', r'$C_{lit}$', r'$C_{som}$']
    lower_bnd = np.array([1e-5, 0.3, 0.01,
                          0.01, 1.0001, 2.5e-5,
                          1e-4, 1e-4, 1e-7,
                          0.018, 10., 60.,
                          0.01, 10., 220.,
                          10., 10., 10.,
                          1e-4, 10., 100.,
                          10., 100.,])

    upper_bnd = np.array([1e-2, 0.7, 0.5,
                          0.5, 10., 1e-3,
                          1e-2, 1e-2, 1e-3,
                          0.1, 100., 185.,
                          0.5, 100., 332.,
                          170., 400., 1000.,
                          1000., 1000., 1e5,
                          1000., 2e5])
    col_lab = ['xb', 'xa_east', 'xa_west', 'low_bnd', 'high_bnd']
    dat = np.hstack([np.array([xb]).T, np.array([xa_east]).T, np.array([xa_west]).T, np.array([lower_bnd]).T,
                    np.array([upper_bnd]).T])
    fig, ax = plt.subplots(nrows=1, ncols=1)
    ax.axis('tight')
    ax.axis('off')
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    ax.table(cellText=dat, rowLabels=keys, colLabels=col_lab, loc='center',)
    return ax, fig


def plot_infmat(infmat, cmin=-0.3, cmax=0.3):
    """Plots influence matrix.
    """
    sns.set(style="ticks")
    sns.set_context('poster', font_scale=1.2)
    fig, ax = plt.subplots(figsize=(11,9))
    #ax.set_aspect('equal')
    plt.imshow(infmat, interpolation='nearest', cmap='bwr', vmin=cmin, vmax=cmax, aspect='auto')
    plt.colorbar(label='Observation influence')
    ax.set_ylabel('Day of year')
    ax.set_xlabel('Day of year')
    #sns.heatmap(rmat, ax=ax, xticklabels=np.arange(len(rmat)), yticklabels=np.arange(len(rmat)))
    return ax, fig


def plot_gaussian_dist(mu, sigma, bounds, xt=None, axx=None):
    """
    Plots a Gausian
    :param mu: mean
    :param sigma: standard deviation
    :param bounds: paramter range
    :param truth: optional truth value
    :param axx: optional axes
    :return: plot
    """
    points = np.linspace(bounds[0], bounds[1], 10000)

    if axx == None:
        if type(mu) is list:
            for m in len(mu):
                plt.plot(points, mlab.normpdf(points, mu[m], sigma[m]))
        else:
            plt.plot(points, mlab.normpdf(points, mu, sigma))
        plt.axvline(xt, linestyle='--', linewidth=50, ms=10)
    else:
        if type(mu) is list:
            for m in len(mu):
                axx.plot(points, mlab.normpdf(points, mu[m], sigma[m]))
        else:
            axx.plot(points, mlab.normpdf(points, mu, sigma))
        axx.axvline(xt, linestyle='--')
        return axx


def plot_many_guassian(mulst, siglst, bndlst, mulst2=None, siglst2=None, truth=None):
    matplotlib.rcParams.update({'figure.autolayout': True})
    sns.set_context('paper')
    #sns.set_style('ticks')
    # define the figure size and grid layout properties
    figsize = (15, 10)
    cols = 5
    gs = gridspec.GridSpec(len(mulst) // cols + 1, cols)

    # plot each markevery case for linear x and y scales
    fig1 = plt.figure(num=1, figsize=figsize)
    ax = []
    keys = [r'$\theta_{min}$', r'$f_{auto}$', r'$f_{fol}$', r'$f_{roo}$', r'$c_{lspan}$', r'$\theta_{woo}$',
            r'$\theta_{roo}$', r'$\theta_{lit}$', r'$\theta_{som}$', r'$\Theta$', r'$c_{eff}$', r'$d_{onset}$',
            r'$f_{lab}$', r'$c_{ronset}$', r'$d_{fall}$', r'$c_{rfall}$', r'$c_{lma}$', r'$C_{lab}$', r'$C_f$',
            r'$C_r$', r'$C_w$', r'$C_l$', r'$C_s$']
    for i, case in enumerate(keys):
        row = (i // cols)
        col = i % cols
        ax.append(fig1.add_subplot(gs[row, col]))
        ax[-1].set_title(case)
        if truth is not None:
            plot_gaussian_dist(mulst[i], siglst[i], bndlst[i], axx=ax[-1], xt=truth[i])
        else:
            plot_gaussian_dist(mulst[i], siglst[i], bndlst[i], axx=ax[-1])
        ax[-1].set_xlim((bndlst[i][0], bndlst[i][1]))
        if mulst2 is not None:
            plot_gaussian_dist(mulst2[i], siglst2[i], bndlst[i], axx=ax[-1])


def plot_rmat(rmat):
    """Plots a R matrix.
    """
    sns.set(style="whitegrid")
    sns.set_context('poster', font_scale=1.2)
    fig, ax = plt.subplots(figsize=(11,9))
    ax.set_aspect('equal')
    sns.heatmap(rmat, ax=ax, vmax=1., xticklabels=False, yticklabels=False,
                linewidths=.5, cbar=True, cbar_kws={'label': 'Correlation'})
    #sns.heatmap(rmat, ax=ax, xticklabels=np.arange(len(rmat)), yticklabels=np.arange(len(rmat)))
    return ax, fig


def plot_bmat(bmat):
    """Plots a B matrix.
    """
    sns.set(style="whitegrid")
    sns.set_context('poster', font_scale=1.2)
    fig, ax = plt.subplots(figsize=(11,9))
    ax.set_aspect('equal')
    keys = [r'$\theta_{min}$', r'$f_{auto}$', r'$f_{fol}$', r'$f_{roo}$', r'$c_{lspan}$', r'$\theta_{woo}$',
            r'$\theta_{roo}$', r'$\theta_{lit}$', r'$\theta_{som}$', r'$\Theta$', r'$c_{eff}$', r'$d_{onset}$',
            r'$f_{lab}$', r'$c_{ronset}$', r'$d_{fall}$', r'$c_{rfall}$', r'$c_{lma}$', r'$C_{lab}$', r'$C_{fol}$',
            r'$C_{roo}$', r'$C_{woo}$', r'$C_{lit}$', r'$C_{som}$']
    ax.set_yticks(np.arange(23))
    ax.set_yticklabels(keys)
    ax.set_xticks(np.arange(23))
    ax.set_xticklabels(keys, rotation=90)
    cmap = sns.diverging_palette(220, 10, as_cmap=True)
    mask = np.eye(23, dtype=bool)
    sns.heatmap(bmat, xticklabels=keys, yticklabels=keys, ax=ax,
                cmap=cmap, vmax=.6, square=True, linewidths=.5, cbar=True,
                cbar_kws={'label': 'Correlation'})

    #ax.set_label('Correlation')
    return ax, fig


# Paper Plots disturbance

def plot_east_west_paper(ob, xa_east, xa_west, d_e, d_w, e_ens, w_ens, y_label='None', y_lim='None', axes='None'):
    ob_ens_e = ob_plist(d_e, e_ens, ob)
    ob_ens_w = ob_plist(d_w, w_ens, ob)
    ob_std_e = ob_mean_std(ob_ens_e)[1]
    ob_std_w = ob_mean_std(ob_ens_w)[1]
    return plot_obs_east_west(ob, xa_east, xa_west, d_e, d_w, y_label=y_label, ob_std_e=ob_std_e, ob_std_w=ob_std_w,
                              y_lim=y_lim, axes=axes)


def plot_east_west_paper_cum(ob, xa_east, xa_west, d_e, d_w, e_ens, w_ens, y_label='None', y_lim='None'):
    ob_ens_e = ob_plist(d_e, e_ens, ob)
    ob_ens_w = ob_plist(d_w, w_ens, ob)
    ob_std_e = ob_mean_std_cum(ob_ens_e)[1]
    ob_std_w = ob_mean_std_cum(ob_ens_w)[1]
    return plot_obs_east_west_cum(ob, xa_east, xa_west, d_e, d_w, y_label=y_label, ob_std_e=ob_std_e, ob_std_w=ob_std_w,
                              y_lim=y_lim)


def plot_east_west_paper_part(xa_east, xa_west, d_e, d_w, e_ens, w_ens, y_label='NEE partitioning', y_lim='None',
                              axes='None'):
    nee_ens_e = ob_plist(d_e, e_ens, 'nee')
    nee_ens_w = ob_plist(d_w, w_ens, 'nee')
    gpp_ens_e = ob_plist(d_e, e_ens, 'gpp')
    gpp_ens_w = ob_plist(d_w, w_ens, 'gpp')
    rt_ens_e = ob_plist(d_e, e_ens, 'rt')
    rt_ens_w = ob_plist(d_w, w_ens, 'rt')
    ob_std_e = []
    ob_std_w = []
    ob_std_e.append(ob_mean_std_cum(nee_ens_e)[1])
    ob_std_w.append(ob_mean_std_cum(nee_ens_w)[1])
    ob_std_e.append(ob_mean_std_cum(gpp_ens_e)[1])
    ob_std_w.append(ob_mean_std_cum(gpp_ens_w)[1])
    ob_std_e.append(ob_mean_std_cum(rt_ens_e)[1])
    ob_std_w.append(ob_mean_std_cum(rt_ens_w)[1])
    return plot_obs_east_west_cum_part(xa_east, xa_west, d_e, d_w, ob_std_e, ob_std_w, y_label=y_label, y_lim=y_lim,
                                       axes=axes)


def plot_east_west_paper2(ob, xa_east, xa_west, d_e, d_w, y_label='None'):
    e_ens = pickle.load(open('east_ens.p', 'r'))
    w_ens = pickle.load(open('west_ens.p', 'r'))
    ob_ens_e = ob_ensemble(d_e, e_ens, ob)
    ob_ens_w = ob_ensemble(d_w, w_ens, ob)
    ob_std_e = ob_mean_std(ob_ens_e)[1]
    ob_std_w = ob_mean_std(ob_ens_w)[1]
    return plot_obs_east_west(ob, xa_east, xa_west, d_e, d_w, y_label=y_label, ob_std_e=ob_std_e, ob_std_w=ob_std_w)

def test_bnds(dC, pvals):
    """
    Test if a parameter set falls within given bounds.
    :param dC: DALEC2 data class
    :param pvals: parameter set to test
    :return: True or False (Pass or Fail)
    """
    tst = 0
    for bnd in enumerate(dC.bnds_tst):
        if bnd[1][0] <= pvals[bnd[0]] <= bnd[1][1]:
            tst += 1
        else:
            continue
    if tst == 23:
        return True
    else:
        return False


def create_ensemble(dC, covmat, pvals):
    ensemble = []
    while len(ensemble) < 500:
        rand = np.random.multivariate_normal(pvals, covmat, 100)
        for xa in rand:
            if test_bnds(dC, xa) == True and pvals_test_uc(dC, xa) == True:
                ensemble.append(xa)
            else:
                continue
    return ensemble


def plist_ens(dC, ensemble):
    m = mc.DalecModel(dC)
    plist_ens = np.ones((len(ensemble), len(dC.I), 23))
    for pvals in enumerate(ensemble):
        plist = m.mod_list(pvals[1])
        plist_ens[pvals[0]] = plist[0:-1]
    return plist_ens


def ob_plist(dC, plist_ens, ob):
    m = mc.DalecModel(dC)
    ob_ens = np.ones((len(plist_ens), len(dC.I)))
    for pvals in enumerate(plist_ens):
        ob_ens[pvals[0]] = m.oblist(ob, pvals[1])
    return ob_ens


def ob_ensemble(dC, ensemble, ob):
    m = mc.DalecModel(dC)
    nee_ens = np.ones((len(ensemble), len(dC.I)))
    for pvals in enumerate(ensemble):
        plist = m.mod_list(pvals[1])
        nee_ens[pvals[0]] = m.oblist(ob, plist)
    return nee_ens


def nee_cumsum_ensemble(dC, ensemble, ob):
    m = mc.DalecModel(dC)
    nee_ens = np.ones((len(ensemble), len(dC.I)))
    nee_cumsum = [np.cumsum(m.oblist(ob, ensemble[x])) for x in xrange(len(ensemble))]
    for x in xrange(len(ensemble)):
        nee_ens[x] = nee_cumsum[x]
    return nee_ens


def ob_mean_std(ob_ens):
    nee_mean = np.nanmean(ob_ens, axis=0)
    nee_std = np.nanstd(ob_ens, axis=0)
    return nee_mean, nee_std


def ob_mean_std_cum(ob_ens):
    nee_cum = np.zeros_like(ob_ens)
    for x in xrange(len(ob_ens)):
        nee_cum[x] = np.cumsum(ob_ens[x])
    nee_mean = np.nanmean(nee_cum, axis=0)
    nee_std = np.nanstd(nee_cum, axis=0)
    return nee_mean, nee_std


def mean_confidence_interval(data, confidence=0.95):
    a = 1.0*np.array(data)
    n = len(a)
    m, se = np.mean(a, axis=0), scipy.stats.sem(a, axis=0)
    h = se * scipy.stats.t.ppf((1+confidence)/2., n-1,)
    return m, m-h, m+h, h


def pvals_test_uc(dC, pvals):
    """ Test if pvals pass constraints.
    """
    #calculate carbon fluxes
    f_auto = pvals[1]
    f_fol = (1-f_auto)*pvals[2]
    f_lab = (1-f_auto-f_fol)*pvals[12]
    f_roo = (1-f_auto-f_fol-f_lab)*pvals[3]
    f_woo = (1 - f_auto - f_fol - f_lab - f_roo)

    #universal constraint tests
    uc = [10*pvals[16] > pvals[18],
          pvals[8] < pvals[7],
          pvals[0] > pvals[8],
          pvals[5] < 1/(365.25*pvals[4]),
          pvals[6] > pvals[8]*np.exp(pvals[9]*dC.t_mean[0]),
          0.2*f_roo < (f_fol+f_lab) < 5*f_roo,
          pvals[14]-pvals[11]>45]
    if all(uc) == True:
        return True
    else:
        return False


def remove_nan_vals(ob_ens):
    nan_arr = np.unique(np.where(np.isnan(ob_ens) == True)[0])
    nee_arr = np.delete(ob_ens, nan_arr, 0)
    return nee_arr



# Misc functions

def cov2cor(X):
    """ Takes a covariance matrix and returns the correlation matrix
    :param X: Covariance matrix
    :return: Correlation matrix
    """
    D = np.zeros_like(X)
    d = np.sqrt(np.diag(X))
    np.fill_diagonal(D, d)
    DInv = np.linalg.inv(D)
    R = np.dot(np.dot(DInv, X), DInv)
    return R


def cov2cor_abs(X):
    """ Takes a covariance matrix and returns the correlation matrix
    :param X: Covariance matrix
    :return: Correlation matrix
    """
    D = np.zeros_like(X)
    d = np.sqrt(np.diag(abs(X)))
    np.fill_diagonal(D, d)
    DInv = np.linalg.inv(D)
    R = np.dot(np.dot(DInv, abs(X)), DInv)
    return R
