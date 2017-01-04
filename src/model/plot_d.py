"""Plotting functions related to dalecv2.
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.mlab as mlab
import matplotlib.gridspec as gridspec
import dalec_d as mc
import data_class_dalec_d as dc
import seaborn as sns

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


def plot_obs_east_west(ob, xa_east, xa_west, d_e, d_w):
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

    ax.plot(d_e.dates, obs_lst_e, color=palette[0], label='East')
    ax.plot(d_w.dates, obs_lst_w, color=palette[2], label='West')
    if ob in d_e.ob_dict.keys():
        ax.plot(d_e.dates, d_e.ob_dict[ob], 'o', color=palette[0], markeredgecolor='black',
                markeredgewidth=0.5)
        ax.plot(d_w.dates, d_w.ob_dict[ob], 'o', color=palette[2], markeredgecolor='black',
                markeredgewidth=0.5)
    plt.legend()
    ax.set_xlabel('Date')
    ax.set_ylabel(ob)
    plt.gcf().autofmt_xdate()
    return ax, fig


def plot_4dvar(ob, dC, xb=None, xa=None, erbars=1, awindl=None, obdict_a=None):
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
    ax.set_ylabel(ob)
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

    one_one = np.arange(int(min(min(y_obs), min(ob_lst))),int(max(max(y_obs), max(ob_lst))))
    plt.plot(one_one, one_one, color=palette[0])

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
    n = 23
    width = 0.35
    ind = np.arange(n)
    #fig = plt.figure()
    #ax = fig.add_subplot(111)
    rects1 = ax.bar(ind, (xa_east-xb)/xb, width, color=sns.xkcd_rgb["faded green"],
                    label='East')
    rects2 = ax.bar(ind+width, (xa_west-xb)/xb, width, color=sns.xkcd_rgb["pale red"],
                    label='West')
    ax.set_ylabel('Normalised analysis increment')
    #ax.set_title('% error in parameter values for xa and xb')
    ax.set_xticks(ind+width)
    keys = [r'$\theta_{min}$', r'$f_{auto}$', r'$f_{fol}$', r'$f_{roo}$', r'$c_{lspan}$', r'$\theta_{woo}$',
            r'$\theta_{roo}$', r'$\theta_{lit}$', r'$\theta_{som}$', r'$\Theta$', r'$c_{eff}$', r'$d_{onset}$',
            r'$f_{lab}$', r'$c_{ronset}$', r'$d_{fall}$', r'$c_{rfall}$', r'$c_{lma}$', r'$C_{lab}$', r'$C_{fol}$',
            r'$C_{roo}$', r'$C_{woo}$', r'$C_{lit}$', r'$C_{som}$']
    ax.set_xticklabels(keys, rotation=90)
    ax.legend()
    return ax, fig


def plot_scatter_twin(ob, pvals, dC, awindl, bfa='a'):
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
    sns.set_style('ticks')
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
    sns.set_style('ticks')
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
    sns.set_style('ticks')
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



