import data_class as dc
import mod_class as mc
import plot as p
import pickle


# ------------------------------------------------------------------------------
# Twin runs perturbed obs
# ------------------------------------------------------------------------------


def day_night_twin_run(start, end, obs, f_name, obs_loc):
    d = dc.DalecDataTwin(start, end, obs, nc_file='../../alice_holt_data/ah_data_daily_test.nc')
    pik_obs = pickle.load(open(obs_loc, 'r'))
    d.ob_dict = pik_obs
    m = mc.DalecModel(d)
    assim_results, xa = m.find_min_tnc_cvt(d.xb, f_name+'_assim_res')
    d2 = dc.DalecDataTwin(start, 2013, obs)
    # Plot 4dvar time series
    ax, fig = p.plot_4dvar_twin('nee', d2, xa=xa)
    fig.savefig(f_name+'_nee.pdf', bbox_inches='tight')
    ax, fig = p.plot_4dvar_twin('nee_day', d2, xa=xa)
    fig.savefig(f_name+'_need.pdf', bbox_inches='tight')
    ax, fig = p.plot_4dvar_twin('nee_night', d2, xa=xa)
    fig.savefig(f_name+'_neen.pdf', bbox_inches='tight')
    ax, fig = p.plot_4dvar_twin('lai', d2, xa=xa)
    fig.savefig(f_name+'_lai.pdf', bbox_inches='tight')

    # Plot scatter plots of obs
    ax, fig = p.plot_scatter_twin('nee', xa, d2, len(d.I), 'f')
    fig.savefig(f_name+'_nee_scat.pdf', bbox_inches='tight')
    ax, fig = p.plot_scatter_twin('nee_day', xa, d2, len(d.I), 'f')
    fig.savefig(f_name+'_need_scat.pdf', bbox_inches='tight')
    ax, fig = p.plot_scatter_twin('nee_night', xa, d2, len(d.I), 'f')
    fig.savefig(f_name+'_neen_scat.pdf', bbox_inches='tight')
    ax, fig = p.plot_scatter_twin('lai', xa, d2, len(d.I), 'f')
    fig.savefig(f_name+'_lai_scat.pdf', bbox_inches='tight')

    # Plot error in analysis and background
    ax, fig = p.plottwinerr(d.x_truth, d.xb, xa)
    fig.savefig(f_name+'_twin_err.pdf', bbox_inches='tight')
    return 'all experimented'

def exps_obs_err(start, end, f_name):
    day_night_twin_run(start, end, 'nee', f_name+'nee', 'obs_exps/'+str(start)+'_'+str(end)+'_nee_twin.p')
    day_night_twin_run(start, end, 'nee_day', f_name+'need', 'obs_exps/'+str(start)+'_'+str(end)+'_nee_day_twin.p')
    day_night_twin_run(start, end, 'nee_night', f_name+'neen', 'obs_exps/'+str(start)+'_'+str(end)+'_nee_night_twin.p')
    day_night_twin_run(start, end, 'nee_day, nee_night', f_name+'needn',
                       'obs_exps/'+str(start)+'_'+str(end)+'_nee_day_night_twin.p')
    return 'all experimented'


# ------------------------------------------------------------------------------
# Twin runs perfect obs
# ------------------------------------------------------------------------------


def day_night_twin_run_no_obs_error(start, end, obs, f_name, its=1000):
    d = dc.DalecDataTwin(start, end, obs, err_scale=0.0)
    m = mc.DalecModel(d)
    assim_results, xa = m.find_min_tnc_cvt(d.xb, f_name+'_assim_res', maxits=its)
    d2 = dc.DalecDataTwin(start, 2013, obs)
    # Plot 4dvar time series
    ax, fig = p.plot_4dvar_twin('nee', d2, xa=xa)
    fig.savefig(f_name+'_nee.pdf', bbox_inches='tight')
    ax, fig = p.plot_4dvar_twin('nee_day', d2, xa=xa)
    fig.savefig(f_name+'_need.pdf', bbox_inches='tight')
    ax, fig = p.plot_4dvar_twin('nee_night', d2, xa=xa)
    fig.savefig(f_name+'_neen.pdf', bbox_inches='tight')
    ax, fig = p.plot_4dvar_twin('lai', d2, xa=xa)
    fig.savefig(f_name+'_lai.pdf', bbox_inches='tight')

    # Plot scatter plots of obs
    ax, fig = p.plot_scatter_twin('nee', xa, d2, len(d.I), 'f')
    fig.savefig(f_name+'_nee_scat.pdf', bbox_inches='tight')
    ax, fig = p.plot_scatter_twin('nee_day', xa, d2, len(d.I), 'f')
    fig.savefig(f_name+'_need_scat.pdf', bbox_inches='tight')
    ax, fig = p.plot_scatter_twin('nee_night', xa, d2, len(d.I), 'f')
    fig.savefig(f_name+'_neen_scat.pdf', bbox_inches='tight')
    ax, fig = p.plot_scatter_twin('lai', xa, d2, len(d.I), 'f')
    fig.savefig(f_name+'_lai_scat.pdf', bbox_inches='tight')

    # Plot error in analysis and background
    ax, fig = p.plottwinerr(d.x_truth, d.xb, xa)
    fig.savefig(f_name+'_twin_err.pdf', bbox_inches='tight')
    return 'all experimented'


def exps_no_obs_err(start, end, f_name):
    day_night_twin_run_no_obs_error(start, end, 'nee', f_name+'nee', 10000)
    day_night_twin_run_no_obs_error(start, end, 'nee_day', f_name+'need', 10000)
    day_night_twin_run_no_obs_error(start, end, 'nee_night', f_name+'neen', 10000)
    day_night_twin_run_no_obs_error(start, end, 'nee_day, nee_night', f_name+'needn', 10000)
    return 'all experimented'
