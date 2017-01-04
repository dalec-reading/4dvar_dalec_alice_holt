"""Dalecv2 model class takes a data class and then uses functions to run the
dalecv2 model.
"""
import numpy as np
import scipy.optimize as spop
import pickle
import algopy
import emcee
import joblib as jl
import random as rand
import multiprocessing
import collections as col


class DalecModel():

    def __init__(self, dataclass, time_step=0, strtrun=0):
        """dataClass and timestep at which to run the dalecv2 model.
        """
        self.dC = dataclass
        self.x = time_step
        if self.dC.k is None:
            self.lenrun = self.dC.len_run
        else:
            self.lenrun = self.dC.len_run*self.dC.k
        self.xb = self.dC.xb
        self.modcoston = True
        self.modobdict = {'gpp': self.gpp, 'nee': self.nee, 'nee_day': self.nee_day,
                          'nee_night': self.nee_night, 'rt': self.rt,
                          'cf': self.cf, 'clab': self.clab, 'c_roo': self.cr,
                          'c_woo': self.cw, 'cl': self.cl, 'cs': self.cs,
                          'lai': self.lai, 'clma': self.clma,'rh': self.rh, 'ra': self.ra,
                          'd_onset': self.d_onset, 'phi_onset': self.phi_on,
                          'phi_fall': self.phi_off}
        self.startrun = strtrun
        self.endrun = self.lenrun
        self.yoblist, self.yerroblist, self.ytimestep = self.obscost()
        self.rmatrix = self.rmat(self.yerroblist)
        self.obs_time_step = self.no_obs_at_time()
        self.diag_b = np.diag(np.diag(self.dC.B))
        self.b_tilda = np.dot(np.dot(np.linalg.inv(np.sqrt(self.diag_b)), self.dC.B),
                              np.linalg.inv(np.sqrt(self.diag_b)))
        self.nume = 100
        self.names = self.dC.edinburgh_mean.dtype.names


# ------------------------------------------------------------------------------
# Model functions
# ------------------------------------------------------------------------------

    @staticmethod
    def fit_polynomial(ep, mult_fac):
        """ Polynomial used to find phi_f and phi (offset terms used in
        phi_onset and phi_fall), given an evaluation point for the polynomial
        and a multiplication term.
        :param ep: evaluation point
        :param mult_fac: multiplication term
        :return: fitted polynomial value
        """
        cf = [2.359978471e-05, 0.000332730053021, 0.000901865258885,
              -0.005437736864888, -0.020836027517787, 0.126972018064287,
              -0.188459767342504]
        poly_val = cf[0]*ep**6 + cf[1]*ep**5 + cf[2]*ep**4 + cf[3]*ep**3 + cf[4]*ep**2 + \
            cf[5]*ep**1 + cf[6]*ep**0
        phi = poly_val*mult_fac
        return phi

    def temp_term(self, Theta, temperature):
        """ Calculates the temperature exponent factor for carbon pool
        respiration's given a value for Theta parameter.
        :param Theta: temperature dependence exponent factor
        :return: temperature exponent respiration
        """
        temp_term = np.exp(Theta*temperature)
        return temp_term

    def acm(self, cf, clma, ceff, acm):
        """ Aggregated canopy model (ACM) function
        ------------------------------------------
        Takes a foliar carbon (cf) value, leaf mass per area (clma) and canopy
        efficiency (ceff) and returns the estimated value for Gross Primary
        Productivity (gpp) of the forest at that time.
        :param cf: foliar carbon (g C m-2)
        :param clma: leaf mass area (g C m-2_
        :param ceff: canopy efficiency parameter
        :return: GPP value
        """
        t_range = 0.5*(self.dC.t_max[self.x] - self.dC.t_min[self.x])
        L = cf / clma
        q = acm[1] - acm[2]
        gc = (abs(self.dC.phi_d))**acm[8] / \
             (t_range + acm[4]*self.dC.R_tot)
        p = ((ceff*L) / gc)*np.exp(acm[6]*self.dC.t_max[self.x])
        ci = 0.5*(self.dC.ca + q - p + np.sqrt((self.dC.ca + q - p)**2
                  - 4*(self.dC.ca*q - p*acm[1])))
        E0 = (acm[5]*L**2) / (L**2 + acm[7])
        delta = -23.4*np.cos((360.*(self.dC.D[self.x] + 10) / 365.) *
                             (np.pi/180.))*(np.pi/180.)
        s = 24*np.arccos((- np.tan(self.dC.lat)*np.tan(delta))) / np.pi
        if s >= 24.:
            s = 24.
        elif s <= 0.:
            s = 0.
        else:
            s = s
        gpp = (E0*self.dC.I[self.x]*gc*(self.dC.ca - ci))*(acm[0]*s +
                                                           acm[3]) / \
              (E0*self.dC.I[self.x] + gc*(self.dC.ca - ci))
        return gpp

    def phi_onset(self, d_onset, cronset):
        """Leaf onset function (controls labile to foliar carbon transfer)
        takes d_onset value, cronset value and returns a value for phi_onset.
        """
        release_coeff = np.sqrt(2.)*cronset / 2.
        mag_coeff = (np.log(1.+1e-3) - np.log(1e-3)) / 2.
        offset = self.fit_polynomial(1+1e-3, release_coeff)
        phi_onset = (2. / np.sqrt(np.pi))*(mag_coeff / release_coeff) * \
            np.exp(-(np.sin((self.dC.D[self.x] - d_onset + offset) /
                     self.dC.radconv)*(self.dC.radconv / release_coeff))**2)
        return phi_onset

    def phi_fall(self, d_fall, crfall, clspan):
        """Leaf fall function (controls foliar to litter carbon transfer) takes
        d_fall value, crfall value, clspan value and returns a value for
        phi_fall.
        """
        release_coeff = np.sqrt(2.)*crfall / 2.
        mag_coeff = (np.log(clspan) - np.log(clspan - 1.)) / 2.
        offset = self.fit_polynomial(clspan, release_coeff)
        phi_fall = (2. / np.sqrt(np.pi))*(mag_coeff / release_coeff) * \
            np.exp(-(np.sin((self.dC.D[self.x] - d_fall + offset) /
                   self.dC.radconv)*self.dC.radconv / release_coeff)**2)
        return phi_fall

    def dalecv2(self, p):
        """DALECV2 carbon balance model
        -------------------------------
        evolves carbon pools to the next time step, taking the 6 carbon pool
        values and 17 parameters at time t and evolving them to time t+1.
        Outputs both the 6 evolved C pool values and the 17 constant parameter
        values.

        phi_on = phi_onset(d_onset, cronset)
        phi_off = phi_fall(d_fall, crfall, clspan)
        gpp = acm(cf, clma, ceff)
        temp = temp_term(Theta)

        clab2 = (1 - phi_on)*clab + (1-f_auto)*(1-f_fol)*f_lab*gpp
        cf2 = (1 - phi_off)*cf + phi_on*clab + (1-f_auto)*f_fol*gpp
        cr2 = (1 - theta_roo)*cr + (1-f_auto)*(1-f_fol)*(1-f_lab)*f_roo*gpp
        cw2 = (1 - theta_woo)*cw + (1-f_auto)*(1-f_fol)*(1-f_lab)*(1-f_roo)*gpp
        cl2 = (1-(theta_lit+theta_min)*temp)*cl + theta_roo*cr + phi_off*cf
        cs2 = (1 - theta_som*temp)*cs + theta_woo*cw + theta_min*temp*cl
        """
        out = algopy.zeros(len(p), dtype=p['cf'])
        # ACM
        gpp = self.acm(p['cf'], p['clma'], p['ceff'], self.dC.acm)
        # Labile release and leaf fall factors
        phi_on = self.phi_onset(p['d_onset'], p['cronset'])
        phi_off = self.phi_fall(p['d_fall'], p['crfall'], p['clspan'])
        # Temperature factor
        temp = self.temp_term(p['Theta'], self.dC.t_mean[self.x])

        out[-6] = (1-phi_on)*p['clab'] + (1-p['f_auto'])*(1-p['f_fol'])*p['f_lab']*gpp
        out[-5] = (1-phi_off)*p['cf'] + phi_on*p['clab'] + (1-p['f_auto'])*(1-p['f_fol'])*p['f_lab']*gpp
        out[-4] = (1-p['theta_roo'])*p['cr'] + (1-p['f_auto'])*(1-p['f_fol'])*(1-p['f_lab'])*p['f_roo']*gpp
        out[-3] = (1-p['theta_woo'])*p['cw'] + (1-p['f_auto'])*(1-p['f_fol'])*(1-p['f_lab'])*(1-p['f_roo'])*gpp
        out[-2] = (1-(p['theta_lit']+p['theta_min'])*temp)*p['cl'] + p['theta_roo']*p['cr'] + phi_off*p['cf']
        out[-1] = (1-p['theta_som']*temp)*p['cs'] + p['theta_woo']*p['cw'] + p['theta_min']*temp*p['cl']
        for x in xrange(len(p)-6):
            out[x] = p[self.names[x]]
        return out

    def dalecv2_diff(self, p):
        """DALECV2 carbon balance model
        -------------------------------
        evolves carbon pools to the next time step, taking the 6 carbon pool
        values and 17 parameters at time t and evolving them to time t+1.
        Outputs both the 6 evolved C pool values and the 17 constant parameter
        values.

        phi_on = phi_onset(d_onset, cronset)
        phi_off = phi_fall(d_fall, crfall, clspan)
        gpp = acm(cf, clma, ceff)
        temp = temp_term(Theta)

        clab2 = (1 - phi_on)*clab + (1-f_auto)*(1-f_fol)*f_lab*gpp
        cf2 = (1 - phi_off)*cf + phi_on*clab + (1-f_auto)*f_fol*gpp
        cr2 = (1 - theta_roo)*cr + (1-f_auto)*(1-f_fol)*(1-f_lab)*f_roo*gpp
        cw2 = (1 - theta_woo)*cw + (1-f_auto)*(1-f_fol)*(1-f_lab)*(1-f_roo)*gpp
        cl2 = (1-(theta_lit+theta_min)*temp)*cl + theta_roo*cr + phi_off*cf
        cs2 = (1 - theta_som*temp)*cs + theta_woo*cw + theta_min*temp*cl
        """
        out = algopy.zeros(6, dtype=p[-6])

        phi_on = self.phi_onset(p[11], p[13])
        phi_off = self.phi_fall(p[14], p[15], p[4])
        gpp = self.acm(p[18], p[16], p[10], self.dC.acm)
        temp = self.temp_term(p[9], self.dC.t_mean[self.x])

        out[0] = (1 - phi_on)*p[17] + (1-p[1])*(1-p[2])*p[12]*gpp
        out[1] = (1 - phi_off)*p[18] + phi_on*p[17] + (1-p[1])*p[2]*gpp
        out[2] = (1 - p[6])*p[19] + (1-p[1])*(1-p[2])*(1-p[12])*p[3]*gpp
        out[3] = (1 - p[5])*p[20] + (1-p[1])*(1-p[2])*(1-p[12])*(1-p[3])*gpp
        out[4] = (1-(p[7]+p[0])*temp)*p[21] + p[6]*p[19] + phi_off*p[18]
        out[5] = (1 - p[8]*temp)*p[22] + p[5]*p[20] + p[0]*temp*p[21]
        return out

    def jac_dalecv2(self, p):
        """Use algopy reverse mode ad calc jac of dv2.
        """
        mat = np.ones((len(p), len(p)))*-9999.
        mat[0:-6] = np.eye(len(p)-6, len(p))
        p_algo = algopy.UTPM.init_jacobian(p)
        p_algo_dic = self.create_ordered_dic(p_algo)
        mat[-6:] = algopy.UTPM.extract_jacobian(self.dalecv2(p_algo_dic))
        return mat, self.dalecv2(p_algo_dic)

    def jac2_dalecv2(self, p):
        """Use algopy reverse mode ad calc jac of dv2.
        """
        mat = np.ones((len(p), len(p)))*-9999.
        mat[0:-6] = np.eye(len(p)-6, len(p))
        p_algo = algopy.UTPM.init_jacobian(p)
        p_algo_lst = self.create_ordered_lst(p_algo)
        mat[-6:] = algopy.UTPM.extract_jacobian(self.dalecv2_diff(p_algo_lst))
        return mat

    def jac_try(self, p):
        out = algopy.zeros((self.endrun-self.startrun, len(p)), dtype=p[-6])
        self.x = self.startrun
        pp = p
        for t in xrange(self.endrun-self.startrun):
            update = self.dalecv2(self.create_ordered_dic(pp))
            pp = update
            out[t, :] = update
            self.x += 1
        return out

    def jac3_dalecv2(self, p):
        """Use algopy reverse mode ad calc jac of dv2.
        """
        #mat = np.ones((len(p), len(p)))*-9999.
        #mat[0:-6] = np.eye(len(p)-6, len(p))
        p_algo = algopy.UTPM.init_jacobian(p)
        #p_algo_lst = self.create_ordered_lst(p_algo)
        mat = algopy.UTPM.extract_jacobian(self.jac_try(p_algo))
        return mat

    def create_ordered_dic(self, p):
        param_dic = col.OrderedDict()
        idx = 0
        for name in self.dC.param_dict.keys():
            if name in self.names:
                param_dic[name] = p[idx]
                idx += 1
            else:
                param_dic[name] = self.dC.param_dict[name]
        return param_dic

    def create_ordered_lst(self, p):
        param = [-9999]*23
        idx = 0
        for name in enumerate(self.dC.param_dict.keys()):
            if name[1] in self.names:
                param[name[0]] = p[idx]
                idx += 1
            else:
                param[name[0]] = self.dC.param_dict[name[1]]
        return param

    def mod_list(self, p):
        """Creates an array of evolving model values using dalecv2 function.
        Takes a list of initial param values.
        """
        param = np.array(self.create_ordered_lst(p))
        mod_list = np.concatenate((np.array([param]),
                                  np.ones((self.endrun-self.startrun, len(param)))*-9999.))

        self.x = self.startrun
        for t in xrange(self.endrun-self.startrun):
            mod_param = mod_list[t]
            mod_param[-6:] = self.dalecv2(self.create_ordered_dic(mod_param))
            mod_list[(t+1)] = mod_param
            self.x += 1

        self.x -= self.endrun
        return mod_list

    def linmod_list(self, p):
        """Creates an array of linearized models (Mi's) taking a list of
        initial param values and a run length (lenrun).
        """
        param = np.array(self.create_ordered_lst(p))
        mod_list = np.concatenate((np.array([param]),
                                  np.ones((self.endrun-self.startrun, len(param)))*-9999.))
        matlist = np.ones((self.endrun-self.startrun, len(p), len(p)))*-9999.

        self.x = self.startrun
        for t in xrange(self.endrun-self.startrun):
            mod_param = mod_list[t]
            c_pool_update = self.dalecv2_diff(mod_param)
            mod_param[-6:] = c_pool_update
            mod_list[(t+1)] = mod_param
            matlist[t] = self.jac2_dalecv2(p)
            p[-6:] = c_pool_update
            self.x += 1

        self.x -= self.endrun
        return mod_list, matlist

    @staticmethod
    def mfac(matlist, timestep):
        """matrix factorial function, takes a list of matrices and a time step,
        returns the matrix factoral.
        """
        if timestep == -1.:
            return np.eye(23)
        mat = matlist[0]
        for t in xrange(0, timestep):
            mat = np.dot(matlist[t+1], mat)
        return mat

    def evolve_mat(self, mat, matlist):
        evolve_mat = mat
        for m in matlist:
            evolve_mat = np.dot(np.dot(m, evolve_mat), m.T)
        return evolve_mat


# ------------------------------------------------------------------------------
# Observation functions
# ------------------------------------------------------------------------------

    def gpp(self, p):
        """Function calculates gross primary production (gpp).
        """
        gpp = self.acm(p[18], p[16], p[10], self.dC.acm)
        return gpp

    def rt(self, p):
        """Function calculates total ecosystem respiration (rec).
        """
        rec = p[1]*self.acm(p[18], p[16], p[10], self.dC.acm) + \
            (p[7]*p[21] + p[8]*p[22])*self.temp_term(p[9], self.dC.t_mean[self.x])
        return rec

    def nee(self, p):
        """Function calculates Net Ecosystem Exchange (nee).
        """
        nee = -(1. - p[1])*self.acm(p[18], p[16], p[10], self.dC.acm) + \
            (p[7]*p[21] + p[8]*p[22])*self.temp_term(p[9], self.dC.t_mean[self.x])
        return nee

    def nee_day(self, p):
        """Function calculates daytime Net Ecosystem Exchange (nee).
        """
        nee = -(1. - (self.dC.day_len[self.x]/24.)*p[1])*self.acm(p[18], p[16], p[10], self.dC.acm) + \
               (self.dC.day_len[self.x]/24.)*(p[7]*p[21] + p[8]*p[22])*self.temp_term(p[9], self.dC.t_day[self.x])
        return nee

    def nee_night(self, p):
        """Function calculates nighttime Net Ecosystem Exchange (nee).
        """
        nee = (self.dC.night_len[self.x]/24.)*p[1]*self.acm(p[18], p[16], p[10], self.dC.acm) + \
              (self.dC.night_len[self.x]/24.)*(p[7]*p[21] + p[8]*p[22])*self.temp_term(p[9], self.dC.t_night[self.x])
        return nee

    def rh(self, p):
        """Fn calculates rh (soilresp+litrep).
        """
        rh = (p[7]*p[21] + p[8]*p[22])*self.temp_term(p[9], self.dC.t_mean[self.x])
        return rh

    def ra(self, p):
        """Fn calculates ra (autotrophic resp.).
        """
        ra = p[1]*self.acm(p[18], p[16], p[10], self.dC.acm)
        return ra

    def lai(self, p):
        """Fn calculates leaf area index (cf/clma).
        """
        lai = p[18] / p[16]
        return lai

    def clma(self, p):
        """Fn calculates clma (carbon leaf mass per area)
        """
        clma = p[16]
        return clma

    def lf(self, p):
        """Fn calulates litter fall.
        """
        lf = self.phi_fall(p[14], p[15], p[4])*p[18]
        return lf

    def clab(self, p):
        """Fn calulates labile carbon.
        """
        clab = p[17]
        return clab

    def cf(self, p):
        """Fn calulates foliar carbon.
        """
        cf = p[18]
        return cf

    def cr(self, p):
        """Fn calulates root carbon.
        """
        cr = p[19]
        return cr

    def cw(self, p):
        """Fn calulates woody biomass carbon.
        """
        cw = p[20]
        return cw

    def cl(self, p):
        """Fn calulates litter carbon.
        """
        cl = p[21]
        return cl

    def cs(self, p):
        """Fn calulates soil organic matter carbon.
        """
        cs = p[22]
        return cs

    def d_onset(self, p):
        """Fn calculates day of leaf on,
        """
        d_onset = p[11]
        return d_onset

    def phi_on(self, p):
        """Fn calculates day of leaf on,
        """
        phi_on_ob = self.phi_onset(p[11], p[13])
        return phi_on_ob

    def phi_off(self, p):
        """Fn calculates day of leaf on,
        """
        phi_off_ob = self.phi_fall(p[14], p[15], p[4])
        return phi_off_ob

    def linob(self, ob, pvals):
        """Function returning jacobian of observation with respect to the
        parameter list. Takes an obs string, a parameters list, a dataClass
        and a time step x.
        """
        dpvals = algopy.UTPM.init_jacobian(pvals)
        return algopy.UTPM.extract_jacobian(self.modobdict[ob](dpvals))

    def oblist(self, ob, mod_list):
        oblist = np.ones(self.endrun-self.startrun)*-9999.
        self.x = self.startrun
        for t in xrange(self.endrun-self.startrun):
            oblist[t] = self.modobdict[ob](mod_list[t])
            self.x += 1
        self.x -= self.endrun
        return oblist

# ------------------------------------------------------------------------------
# Assimilation functions
# ------------------------------------------------------------------------------

    def obscost(self):
        """Function returning list of observations and a list of their
        corresponding error values. Takes observation dictionary and an
        observation error dictionary.
        """
        yoblist = np.array([])
        yerrlist = np.array([])
        ytimestep = np.array([])
        for t in xrange(self.startrun, self.endrun):
            for ob in self.dC.ob_dict.iterkeys():
                if np.isnan(self.dC.ob_dict[ob][t]) != True:
                    yoblist = np.append(yoblist, self.dC.ob_dict[ob][t])
                    yerrlist = np.append(yerrlist,
                                         self.dC.ob_err_dict[ob][t])
                    ytimestep = np.append(ytimestep, t)
        return yoblist, yerrlist, ytimestep

    def hxcost(self, pvallist):
        """Function returning a list of observation values as predicted by the
        DALEC model. Takes a list of model values (pvallist), an observation
        dictionary and a dataClass (dC).
        """
        hx = np.array([])
        self.x = self.startrun
        for t in xrange(self.startrun, self.endrun):
            for ob in self.dC.ob_dict.iterkeys():
                if np.isnan(self.dC.ob_dict[ob][t]) != True:
                    hx = np.append(hx,
                                   self.modobdict[ob](pvallist[t-self.startrun]))
            self.x += 1

        self.x -= self.endrun
        return hx

    def rmat(self, yerr):
        """Returns observation error covariance matrix given a list of
        observation error values.
        """
        r = (yerr**2)*np.eye(len(yerr))
        return r

    def no_obs_at_time(self):
        """ Returns a list of the number of observations at each time step.
        """
        obs_time_step = np.array([])
        self.x = self.startrun
        for t in xrange(self.startrun, self.endrun):
            p = 0
            for ob in self.dC.ob_dict.iterkeys():
                if np.isnan(self.dC.ob_dict[ob][t]) != True:
                    p += 1
            obs_time_step = np.append(obs_time_step, p)
            self.x += 1

        self.x -= self.endrun
        return obs_time_step

    def gradcost2(self, pvals):
        """Gradient of 4DVAR cost fn to be passed to optimization routine.
        Takes an initial state (pvals), an obs dictionary, an obs error
        dictionary, a dataClass and a start and finish time step. Using Lagrange
        multipliers to increase speed, method updated to allow for temporally
        correlated R matrix. Uses method of Lagrange multipliers!
        """
        pvallist, matlist = self.linmod_list(pvals)
        hx, hhat = self.hhat(pvallist)
        r_yhx = np.dot(np.linalg.inv(self.rmatrix), (self.yoblist-hx).T)
        idx1 = len(self.yoblist) - sum(self.obs_time_step[self.lenrun-1:])
        idx2 = len(self.yoblist) - sum(self.obs_time_step[self.lenrun-1+1:])
        obcost = np.dot(hhat[idx1:idx2].T, r_yhx[idx1:idx2])
        for i in xrange(self.lenrun-2, -1, -1):
            if self.obs_time_step[i] != 0:
                idx1 = len(self.yoblist) - sum(self.obs_time_step[i:])
                idx2 = len(self.yoblist) - sum(self.obs_time_step[i+1:])
                obcost = np.dot(matlist[i].T, obcost) + np.dot(hhat[idx1:idx2].T, r_yhx[idx1:idx2])
            else:
                obcost = np.dot(matlist[i].T, obcost)

        if self.modcoston is True:
            modcost = np.dot(np.linalg.inv(self.dC.B), (pvals-self.xb).T)
        else:
            modcost = 0
        gradcost = - obcost + modcost
        return gradcost

    def hhat(self, pvallist):
        """Returns a list of observation values as predicted by DALEC (hx) and
        a stacked set of linearzied observation operators (hmat) for use in gradcost2
        fn calculating the gradient of the cost fn using the method of Lagrange multipliers.
        Takes a list of model values (pvallist), a observation dictionary, a list of
        linearized models (matlist) and a dataClass (dC).
        """
        hx = np.array([])
        hhat = []
        self.x = self.startrun
        for t in xrange(self.startrun, self.endrun):
            temp = []
            for ob in self.dC.ob_dict.iterkeys():
                if np.isnan(self.dC.ob_dict[ob][t]) != True:
                    hx = np.append(hx,
                                   self.modobdict[ob](pvallist[t-self.startrun]))
                    temp.append([self.linob(ob, pvallist[t-self.startrun])])
            self.x += 1
            if len(temp) != 0.:
                hhat.append(np.vstack(temp))
            else:
                continue

        self.x -= self.endrun
        return hx, np.vstack(hhat)

    def hmat(self, pvallist, matlist):
        """Returns a list of observation values as predicted by DALEC (hx) and
        a linearzied observation error covariance matrix (hmat). Takes a list
        of model values (pvallist), a observation dictionary, a list of
        linearized models (matlist) and a dataClass (dC).
        """
        hx = np.array([])
        hmat = np.array([])
        self.x = self.startrun
        for t in xrange(self.startrun, self.endrun):
            temp = []
            for ob in self.dC.ob_dict.iterkeys():
                if np.isnan(self.dC.ob_dict[ob][t]) != True:
                    hx = np.append(hx,
                                   self.modobdict[ob](pvallist[t-self.startrun]))
                    temp.append([self.linob(ob, pvallist[t-self.startrun])])
            self.x += 1
            if len(temp) != 0.:
                hmat = np.append(hmat, np.dot(np.vstack(temp),
                                 self.mfac(matlist, t-self.startrun-1)))
            else:
                continue

        self.x -= self.endrun
        hmat = np.reshape(hmat, (len(hmat)/23, 23))
        return hx, hmat

    def modcost(self, pvals):
        """model part of cost fn.
        """
        return np.dot(np.dot((pvals-self.xb), np.linalg.inv(self.dC.B)), (pvals-self.xb).T)

    def obcost(self, pvals):
        """Observational part of cost fn.
        """
        pvallist = self.mod_list(pvals)
        hx = self.hxcost(pvallist)
        return np.dot(np.dot((self.yoblist-hx), np.linalg.inv(self.rmatrix)), (self.yoblist-hx).T)

    def cost(self, pvals):
        """4DVAR cost function to be minimized. Takes an initial state (pvals),
        an observation dictionary, observation error dictionary, a dataClass
        and a start and finish time step.
        """

        ob_cost = self.obcost(pvals)
        if self.modcoston is True:
            mod_cost = self.modcost(pvals)
        else:
            mod_cost = 0
        cost = 0.5*ob_cost + 0.5*mod_cost
        return cost

    def gradcost(self, pvals):
        """Gradient of 4DVAR cost fn to be passed to optimization routine.
        Takes an initial state (pvals), an obs dictionary, an obs error
        dictionary, a dataClass and a start and finish time step.
        """
        pvallist, matlist = self.linmod_list(pvals)
        hx, hmatrix = self.hmat(pvallist, matlist)
        obcost = np.dot(hmatrix.T, np.dot(np.linalg.inv(self.rmatrix),
                                          (self.yoblist-hx).T))
        if self.modcoston is True:
            modcost = np.dot(np.linalg.inv(self.dC.B), (pvals-self.xb).T)
        else:
            modcost = 0
        gradcost = - obcost + modcost
        return gradcost

    def acovmat(self, pvals):
        """Calculates approximation to analysis error covariance matrix
        A = (B^(-1) + H^(T) * R^(-1) * H)^(-1).
        """
        pvallist, matlist = self.linmod_list(pvals)
        hx, hmatrix = self.hmat(pvallist, matlist)
        return np.linalg.inv(np.linalg.inv(self.dC.B) + np.dot(hmatrix.T,
                        np.dot(np.linalg.inv(self.rmatrix), hmatrix)))


# ------------------------------------------------------------------------------
# Information content.
# ------------------------------------------------------------------------------


    def k_gain_mat(self, pvals):
        """ Calculates Kalman gain matrix for info content experiments.
        :param pvals: state vector for assimilation
        :return: Kalman gain hat matrix
        """
        pvallist, matlist = self.linmod_list(pvals)
        hx, hmatrix = self.hmat(pvallist, matlist)
        k_gain = np.dot(np.linalg.inv(np.dot(np.dot(hmatrix.T, np.linalg.inv(self.rmatrix)), hmatrix) +
                                      np.linalg.inv(self.dC.B)),
                        np.dot(hmatrix.T, np.linalg.inv(self.rmatrix)))
        return k_gain

    def influence_mat(self, pvals):
        """Calculates the Influence matrix
        """
        pvallist, matlist = self.linmod_list(pvals)
        hx, hmatrix = self.hmat(pvallist, matlist)
        a_mat = np.linalg.inv(np.linalg.inv(self.dC.B) + np.dot(np.dot(
                              hmatrix.T, np.linalg.inv(self.rmatrix)), hmatrix))
        s_mat = np.dot(np.linalg.inv(self.rmatrix), np.dot(hmatrix, np.dot(
                              a_mat, hmatrix.T)))
        return s_mat

    def influence_mat_cvt(self, pvals):
        """Calculates the Influence matrix
        """
        pvallist, matlist = self.linmod_list(pvals)
        hx, hmatrix = self.hmat(pvallist, matlist)
        k_gain = self.k_gain_mat(pvals)
        s_mat = np.dot(hmatrix, np.dot(np.linalg.inv(np.sqrt(self.diag_b)), k_gain))
        return s_mat

    @staticmethod
    def inf_mat_colsum(infmat):
        """Sums columns of provided infmat
        """
        infmat_sum = np.sum(abs(infmat), axis=1)
        return infmat_sum

    def inf_mat_partition(self, infmat):
        """Partitions summed columns of inf_matrix to correspond to their
        relevant observations.
        """
        infmat_sum = self.inf_mat_colsum(infmat)
        inf_con = {}
        step_size = len(self.dC.ob_dict.keys())
        for i in enumerate(self.dC.ob_dict.keys()):
            inf_con[i[1]] = infmat_sum[i[0]::step_size]
        return inf_con

    def k_gain_partition(self, kgain):
        """Partitions rows of kgain_matrix to correspond to their
        relevant observations.
        """
        keys = ['theta_min', 'f_auto', 'f_fol', 'f_roo', 'clspan', 'theta_woo',
                'theta_roo', 'theta_lit', 'theta_som', 'Theta', 'ceff', 'd_onset',
                'f_lab', 'cronset', 'd_fall', 'crfall', 'clma', 'clab', 'cf', 'cr',
                'cw', 'cl', 'cs']
        kgain_part = {}
        for i in enumerate(keys):
            kgain_part[i[1]] = kgain[:, i[0]]
        return kgain_part


    def analytic_xa(self, pvals):
        """See how good linear approx is.
        """
        pvallist, matlist = self.linmod_list(pvals)
        hx, hmatrix = self.hmat(pvallist, matlist)
        k_mat = np.dot(np.dot(self.dC.B, hmatrix.T),
                       np.linalg.inv(np.dot(np.dot(hmatrix, self.dC.B), hmatrix.T) + self.rmatrix))
        xa = pvals + np.dot(k_mat, (self.yoblist - hx))
        return xa

    def sic(self, pvals):
        acovmat = self.acovmat(pvals)
        return 0.5*np.log(np.linalg.det(self.dC.B)/np.linalg.det(acovmat))

    def dfs(self, pvals):
        acovmat = self.acovmat(pvals)
        return 23 - np.matrix.trace(np.dot(np.linalg.inv(self.dC.B), acovmat))

    def single_sic(self, pvals, ob):
        hmat = self.linob(ob, pvals)
        rmat = self.dC.errdict[ob]
        bmat = self.dC.B
        acovmat = np.linalg.inv(np.linalg.inv(self.dC.B) + np.dot(hmat.T, np.dot((rmat**2)**(-1), hmat)))
        # acovmat = np.linalg.inv(np.linalg.inv(self.b_tilda) + np.dot(np.sqrt(self.diag_b), np.dot(hmat.T,
        #                 np.dot(np.dot((rmat**2)**(-1), hmat), np.sqrt(self.diag_b)))))
        # bmat = self.b_tilda
        return 0.5*np.log(np.linalg.det(bmat)/np.linalg.det(acovmat))

    def single_dfs(self, pvals, ob):
        hmat = self.linob(ob, pvals)
        rmat = self.dC.errdict[ob]
        acovmat = np.linalg.inv(np.linalg.inv(self.dC.B) + np.dot(hmat.T, np.dot((rmat**2)**(-1), hmat)))
        return 23 - np.matrix.trace(np.dot(np.linalg.inv(self.dC.B), (acovmat)))

    def sic_time(self, pvals, ob):
        sic_lst = np.ones(self.endrun-self.startrun)*-9999.
        mod_list = self.mod_list(pvals)
        self.x = self.startrun
        for t in xrange(self.endrun-self.startrun):
            sic_lst[t] = self.single_sic(mod_list[t], ob)
            self.x += 1
        self.x -= self.endrun
        return sic_lst

    def dfs_time(self, pvals, ob):
        dfs_lst = np.ones(self.endrun-self.startrun)*-9999.
        mod_list = self.mod_list(pvals)
        self.x = self.startrun
        for t in xrange(self.endrun-self.startrun):
            dfs_lst[t] = self.single_dfs(mod_list[t], ob)
            self.x += 1
        self.x -= self.endrun
        return dfs_lst


# ------------------------------------------------------------------------------
# CVT and implied B.
# ------------------------------------------------------------------------------

    def modcost_cvt(self, zvals):
        """model part of cost fn.
        """
        return np.dot(np.dot(zvals, np.linalg.inv(self.b_tilda)), zvals.T)

    def obcost_cvt(self, zvals):
        """Observational part of cost fn.
        """
        pvals = self.zvals2pvals(zvals)
        pvallist = self.mod_list(pvals)
        hx = self.hxcost(pvallist)
        return np.dot(np.dot((self.yoblist-hx), np.linalg.inv(self.rmatrix)), (self.yoblist-hx).T)

    def bndcost_cvt(self, zvals):
        pvals = self.zvals2pvals(zvals)
        return np.dot(np.dot((pvals-self.dC.mid_bnds), np.linalg.inv(self.dC.bnd_cov)), (pvals-self.dC.mid_bnds).T)

    def cost_cvt(self, zvals):
        """4DVAR cost function to be minimized. Takes an initial state (pvals),
        an observation dictionary, observation error dictionary, a dataClass
        and a start and finish time step.
        """
        ob_cost = self.obcost_cvt(zvals)
        if self.modcoston is True:
            mod_cost = self.modcost_cvt(zvals)
        else:
            mod_cost = 0
        cost = 0.5*ob_cost + 0.5*mod_cost  # + 0.5*self.bndcost_cvt(zvals)
        return cost

    def gradcost_cvt(self, zvals):
        """Gradient of 4DVAR cost fn to be passed to optimization routine.
        Takes an initial state (pvals), an obs dictionary, an obs error
        dictionary, a dataClass and a start and finish time step.
        """
        pvals = self.zvals2pvals(zvals)
        pvallist, matlist = self.linmod_list(pvals)
        hx, hmatrix = self.hmat(pvallist, matlist)
        obcost = np.dot(np.sqrt(self.diag_b).T, np.dot(hmatrix.T, np.dot(np.linalg.inv(self.rmatrix),
                                          (self.yoblist-hx).T)))
        if self.modcoston is True:
            modcost = np.dot(np.linalg.inv(self.b_tilda), zvals.T)
        else:
            modcost = 0
        # bnd_cost = np.dot(np.linalg.inv(self.dC.bnd_cov), (pvals-self.dC.mid_bnds).T)
        gradcost = - obcost + modcost  # + bnd_cost
        return gradcost

    def gradcost2_cvt(self, zvals):
        """Gradient of 4DVAR cost fn to be passed to optimization routine.
        Takes an initial state (pvals), an obs dictionary, an obs error
        dictionary, a dataClass and a start and finish time step. Using Lagrange
        multipliers to increase speed, method updated to allow for temporally
        correlated R matrix. Uses method of Lagrange multipliers!
        """
        pvals = self.zvals2pvals(zvals)
        pvallist, matlist = self.linmod_list(pvals)
        hx, hhat = self.hhat(pvallist)
        r_yhx = np.dot(np.linalg.inv(self.rmatrix), (self.yoblist-hx).T)
        idx1 = len(self.yoblist) - sum(self.obs_time_step[self.lenrun-1:])
        idx2 = len(self.yoblist) - sum(self.obs_time_step[self.lenrun-1+1:])
        obcost = np.dot(hhat[idx1:idx2].T, r_yhx[idx1:idx2])
        for i in xrange(self.lenrun-2, -1, -1):
            if self.obs_time_step[i] != 0:
                idx1 = len(self.yoblist) - sum(self.obs_time_step[i:])
                idx2 = len(self.yoblist) - sum(self.obs_time_step[i+1:])
                obcost = np.dot(matlist[i].T, obcost) + np.dot(hhat[idx1:idx2].T, r_yhx[idx1:idx2])
            else:
                obcost = np.dot(matlist[i].T, obcost)
        obcost = np.dot(np.sqrt(self.diag_b).T, obcost)

        if self.modcoston is True:
            modcost = np.dot(np.linalg.inv(self.b_tilda), zvals.T)
        else:
            modcost = 0
        gradcost = - obcost + modcost
        return gradcost

    def pvals2zvals(self, pvals):
        """Convert x_0 state to z_0 state for CVT with DALEC.
        """
        Bsqrt = np.linalg.inv(np.sqrt(self.diag_b))
        return np.dot(Bsqrt, (pvals - self.xb))

    def zvals2pvals(self, zvals):
        """Convert z_0 to x_0 for CVT.
        """
        return np.dot(np.sqrt(self.diag_b),zvals)+self.xb

    def zvalbnds(self, bnds):
        """Calculates bounds for transformed problem.
        """
        lower_bnds = []
        upper_bnds = []
        for t in bnds:
            lower_bnds.append(t[0])
            upper_bnds.append(t[1])
        zval_lowerbnds = self.pvals2zvals(np.array(lower_bnds))
        zval_upperbnds = self.pvals2zvals(np.array(upper_bnds))
        new_bnds=[]
        for t in xrange(len(bnds)):
            new_bnds.append((zval_lowerbnds[t],zval_upperbnds[t]))
        return tuple(new_bnds)

    def findmin_cvt(self, pvals, bnds='strict', dispp=None, maxits=2000,
                   mini=0, f_tol=-1):
        """Function which minimizes 4DVAR cost fn. Takes an initial state
        (pvals).
        """
        self.xb = pvals
        zvals = self.pvals2zvals(pvals)
        if bnds == 'strict':
            bnds = self.zvalbnds(self.dC.bnds_tst)
        else:
            bnds = bnds
        findmin = spop.fmin_tnc(self.cost_cvt, zvals,
                                fprime=self.gradcost2_cvt, bounds=bnds,
                                disp=dispp, fmin=mini, maxfun=maxits, ftol=f_tol)
        return findmin

    def cvt_hmat(self, pvallist, matlist):
        """
        Calculates the normalised \hat{H} matrix for the CVT case
        :param pvallist: list of model evolved parameter values
        :param matlist: list of linearised models
        :return: normalised \hat{H}
        """
        hx, hmat = self.hmat(pvallist, matlist)
        obs_mat = np.dot(np.dot(np.linalg.inv(np.sqrt(self.rmatrix)), hmat), np.sqrt(self.diag_b))
        return obs_mat

    def cvt_a_covmat(self, pvals):
        """Calculates approximation to analysis error covariance matrix
        A = (B^(-1) + H^(T) * R^(-1) * H)^(-1).
        """
        pvallist, matlist = self.linmod_list(pvals)
        hx, hmatrix = self.hmat(pvallist, matlist)
        return np.linalg.inv(np.linalg.inv(self.b_tilda) + np.dot(np.sqrt(self.diag_b), np.dot(hmatrix.T,
                        np.dot(np.dot(np.linalg.inv(self.rmatrix), hmatrix), np.sqrt(self.diag_b)))))

    def cvt_sic(self, pvals):
        acovmat = self.cvt_a_covmat(pvals)
        return 0.5*np.log(np.linalg.det(self.b_tilda)/np.linalg.det(acovmat))

    def cvt_dfs(self, pvals):
        acovmat = self.cvt_a_covmat(pvals)
        return 23 - np.matrix.trace(np.dot(np.linalg.inv(self.b_tilda), acovmat))


# ------------------------------------------------------------------------------
# Minimization Routines.
# ------------------------------------------------------------------------------

    def find_min_tnc(self, pvals, bnds='strict', dispp=None, maxits=2000,
                     mini=0, f_tol=-1):
        """Function which minimizes 4DVAR cost fn. Takes an initial state
        (pvals).
        """
        self.xb = pvals
        if bnds == 'strict':
            bnds = self.dC.bnds
        else:
            bnds = bnds
        find_min = spop.fmin_tnc(self.cost, pvals,
                                fprime=self.gradcost2, bounds=bnds,
                                disp=dispp, fmin=mini, maxfun=maxits, ftol=f_tol)
        return find_min

    def find_min_tnc_cvt(self, pvals, f_name=None, bnds='strict', dispp=5, maxits=1000,
                         mini=0, f_tol=1e-4):
        """Function which minimizes 4DVAR cost fn. Takes an initial state
        (pvals).
        """
        self.xb = pvals
        if bnds == 'strict':
            bnds = self.zvalbnds(self.dC.bnds_tst)
        else:
            bnds = bnds
        zvals = self.pvals2zvals(pvals)
        find_min = spop.fmin_tnc(self.cost_cvt, zvals,
                                fprime=self.gradcost2_cvt, bounds=bnds,
                                disp=dispp, fmin=mini, maxfun=maxits, ftol=f_tol)
        xa = self.zvals2pvals(find_min[0])
        if f_name != None:
            self.pickle_exp(pvals, find_min, xa, f_name)
        return find_min, xa

    def findminglob(self, pvals, meth='TNC', bnds='strict', it=300,
                    stpsize=0.5, temp=1., displ=True, maxits=3000):
        """Function which minimizes 4DVAR cost fn. Takes an initial state
        (pvals), an obs dictionary, an obs error dictionary, a dataClass and
        a start and finish time step.
        """
        if bnds == 'strict':
            bnds = self.dC.bnds
        else:
            bnds = bnds
        findmin = spop.basinhopping(self.cost, pvals, niter=it,
                                    minimizer_kwargs={'method': meth, 'bounds': bnds,
                                                      'jac': self.gradcost,
                                                      'options': {'maxiter': maxits}},
                                    stepsize=stpsize, T=temp, disp=displ)
        return findmin

    def ensemble(self, pvals):
        """Ensemble 4DVAR run for twin experiments.
        """
        ensempvals = np.ones((self.nume, 23))
        for x in xrange(self.nume):
            ensempvals[x] = self.dC.randompert(pvals)

        assim_results = [self.find_min_tnc(ensemp, dispp=False) for ensemp in
                         ensempvals]

        xalist = [assim_results[x][0] for x in xrange(self.nume)]

        return ensempvals, xalist, assim_results

    def var_ens(self, size_ens=10):
        edc_ens = pickle.load(open('misc/edc_param_ensem.p', 'r'))
        param_ens = rand.sample(edc_ens, size_ens)
        num_cores = multiprocessing.cpu_count()
        output = jl.Parallel(n_jobs=num_cores)(jl.delayed(self.find_min_tnc_cvt)(self, pval) for pval in param_ens)
        return output


# ------------------------------------------------------------------------------
# Cycled 4D-Var.
# ------------------------------------------------------------------------------

    def cycle_4dvar(self, pvals, lenwind, numbwind, lenrun):
        """Cycle 4Dvar windows and see their effect on predicting future obs.
        """
        xb = [pvals]
        xa = []
        self.startrun = 0
        self.endrun = lenwind
        for x in xrange(numbwind):
            self.yoblist, self.yerroblist, ytimstep = self.obscost()
            self.rmatrix = self.rmat(self.yerroblist)
            xa.append(self.find_min_tnc(xb[x]))
            xb.append(self.mod_list(xa[x][0])[self.endrun-self.startrun])
            self.startrun += lenwind
            self.endrun += lenwind

        self.startrun -= lenwind*numbwind
        self.endrun -= lenwind*numbwind
        conditions = {'pvals': pvals, 'lenwind': lenwind, 'numbwind': numbwind,
                      'lenrun': lenrun}
        return conditions, xb, xa

    def yearly_cycle4dvar(self, pvals):
        """
        Performs cycle DA with windows of year length
        :param pvals: Initial background vector for first assim. window
        :return: list of all xb vectors, list of all xa vectors and minimisation
        output.
        """
        year_lst = np.unique(self.dC.year)
        xb = [pvals]
        xa = []
        for year in enumerate(year_lst):
            year_idx = np.where(self.dC.year == year[1])[0]
            self.startrun = year_idx[0]
            self.endrun = year_idx[-1]
            self.yoblist, self.yerroblist, ytimestep = self.obscost()
            self.rmatrix = self.rmat(self.yerroblist)
            xa.append(self.find_min_tnc_cvt(pvals, f_tol=1e1))
            acovmat = self.acovmat(xa[year[0]][1])
            self.endrun += 1
            pvallst, matlist = self.linmod_list(xa[year[0]][1])
            xb.append(pvallst[-1])
            ev_acovmat = 1.2 * self.evolve_mat(acovmat, matlist)
            # ev_acovmat = self.dC.B
            ev_acovmat[11,11] = self.dC.B[11,11]
            ev_acovmat[13,13] = self.dC.B[13,13]
            ev_acovmat[14,14] = self.dC.B[14,14]
            ev_acovmat[15,15] = self.dC.B[15,15]
            self.diag_b = np.diag(np.diag(ev_acovmat))
            self.b_tilda = np.dot(np.dot(np.linalg.inv(np.sqrt(self.diag_b)),ev_acovmat),
                                  np.linalg.inv(np.sqrt(self.diag_b)))
            print xa[year[0]][1]
            #Change B too, use corr B to start with then evolve A with M
            #for newB. Figure how to create corrR with multiple data streams.

        self.startrun = 0
        self.endrun = self.lenrun
        self.diag_b = np.diag(np.diag(self.dC.B))
        self.b_tilda = np.dot(np.dot(np.linalg.inv(np.sqrt(self.diag_b)),self.dC.B),np.linalg.inv(np.sqrt(self.diag_b)))
        return xb, xa


# ------------------------------------------------------------------------------
# MCMC
# ------------------------------------------------------------------------------

    def randompvals(self, nwalkers):
        """Creates a random list of parameter values for dalec using dataClass
        bounds.
        """
        rndpvals = np.ones((nwalkers, 23))*-9999.
        x=0
        for t in xrange(nwalkers):
            for bnd in self.dC.bnds:
                rndpvals[t, x] = np.random.uniform(bnd[0],bnd[1])
                x += 1
            x-=23

        return rndpvals

    def log_prior(self, pvals):
        tick = 23
        for x in xrange(23):
            if self.dC.bnds[x][0]<pvals[x]<self.dC.bnds[x][1]:
                tick -= 1
        if tick==0:
            return 0.0
        return -np.inf

    def log_posterior(self, pvals):
        lp = self.log_prior(pvals)
        if not np.isfinite(lp):
            return -np.inf
        return self.log_prior(pvals) - self.obcost(pvals)

    def mcmc_run(self, pvals, ndim=23, nwalkers=50, nsteps=10000):
        """
        ndim = 23  # number of parameters in the model
        nwalkers = 50 # number of MCMC walkers
        nburn = 1000  # "burn-in" period to let chains stabilize
        nsteps = 10000  # number of MCMC steps to take
        """
        sampler = emcee.EnsembleSampler(nwalkers, ndim, self.log_posterior)
        sampler.run_mcmc(pvals, nsteps)
        return sampler

    def mean_pval(self, sampler, nburn=1000):
        """
        nburn = 1000  # "burn-in" period to let chains stabilize
        """
        sample = sampler.chain[:, nburn:, :].reshape(-1, 23)
        return sample.mean(0)

    def tstpvals(self, pvals, bnds):
        """Tests pvals to see if they are within the correct bnds.
        """
        x=0
        for bnd in bnds:
            if bnd[0]<pvals[x]<bnd[1]:
                print '%f in bnds' %x
            else:
                print '%f not in bnds' %x
            x += 1
        return pvals

    def tstpvals2(self, pvals):
        tick = 0
        for x in xrange(23):
            if self.dC.bnds5[x][0] >= pvals[x] >= self.dC.bnds5[x][1]:
                print '%f' %x
                tick -= 1
        if tick != 0:
            return np.inf


# ------------------------------------------------------------------------------
# Misc
# ------------------------------------------------------------------------------

    def pickle_obs(self, f_name):
        obs = {}
        obs['obs'] = self.dC.ob_dict
        obs['obs_err'] = self.dC.ob_err_dict
        f = open('obs_exps/'+f_name, 'w')
        pickle.dump(obs, f)
        f.close()
        return 'Observations and error dictionaries pickled!'

    def pickle_exp(self, xb, assim_res, xa, f_name):
        exp = {}
        exp['obs'] = self.dC.ob_dict
        exp['obs_err'] = self.dC.ob_err_dict
        exp['b_mat'] = self.dC.B
        exp['xb'] = xb
        exp['assim_res'] = assim_res
        exp['xa'] = xa
        f = open(f_name, 'w')
        pickle.dump(exp, f)
        f.close()
        return 'Experiment assimilation results pickled!'
