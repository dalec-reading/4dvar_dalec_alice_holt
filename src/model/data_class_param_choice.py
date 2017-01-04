import numpy as np
import collections as col
import re
import random
import mod_class as mc
import netCDF4 as nC
import datetime as dt
import pickle


class DalecData:
    """
    Data class for the DALEC2 model
    """
    def __init__(self, start_date=None, end_date=None, ob_str=None,
                 nc_file='../../alice_holt_data/ah_data_daily_test.nc', delta_t=None, k=None):
        """ Extracts data from netcdf file
        :param start_date: date for model runs to begin as an integer (year) or tuple (year, month, day)
        :param end_date: date for model runs to end as an integer (year) or tuple (year, month, day)
        :param ob_str: string containing observations that will be assimilated
        :param nc_file: location of netcdf file to extract data from
        :param delta_t: time step for model (this functionality has not been added yet)
        :param k: int if you want to repeat data multiple times
        :return:
        """
        # Extract the data
        data = nC.Dataset(nc_file, 'r')
        self.time_units = data.variables['time'].units
        self.start_idx = self.find_date_idx(start_date, data)
        self.end_idx = self.find_date_idx(end_date, data)
        self.dates = nC.num2date(data.variables['time'][self.start_idx:self.end_idx], self.time_units)
        self.year = [date.year for date in self.dates]
        self.k = k
        self.len_run = len(self.dates)
        self.time_step = np.arange(self.len_run)

        # I.C. for carbon pools gCm-2     range
        self.clab = 75.0               # (10,1e3)
        self.cf = 0.0                  # (10,1e3)
        self.cr = 135.0                # (10,1e3)
        self.cw = 14313.0              # (3e3,3e4)
        self.cl = 70.                  # (10,1e3)
        self.cs = 18624.0              # (1e3, 1e5)

        # Parameters for optimization                        range
        self.p1 = 2.5e-4  # theta_min, cl to cs decomp      (1e-5 - 1e-2) day-1
        self.p2 = 0.404  # f_auto, fraction of GPP respired  (0.3 - 0.7)
        self.p3 = 0.01  # f_fol, frac GPP to foliage        (0.01 - 0.5)
        self.p4 = 0.457  # f_roo, frac GPP to fine roots    (0.01 - 0.5)
        self.p5 = 3.  # clspan, leaf lifespan               (1.0001 - 5)
        self.p6 = 4.8e-5  # theta_woo, wood C turnover      (2.5e-5 - 1e-3) day-1
        self.p7 = 6.72e-3  # theta_roo, root C turnover rate(1e-4 - 1e-2) day-1
        self.p8 = 0.024  # theta_lit, litter C turnover     (1e-4 - 1e-2) day-1
        self.p9 = 2.4e-5  # theta_som, SOM C turnover       (1e-7 - 1e-3) day-1
        self.p10 = 0.0193  # Theta, temp dependence exp fact(0.018 - 0.08)
        self.p11 = 5.61843728e+01  # ceff, canopy efficiency param     (10 - 100)
        self.p12 = 140.  # d_onset, clab release date       (1 - 365) (60,150)
        self.p13 = 0.7  # f_lab, frac GPP to clab           (0.01 - 0.5)
        self.p14 = 27.  # cronset, clab release period      (10 - 100)
        self.p15 = 308.  # d_fall, date of leaf fall        (1 - 365) (242,332)
        self.p16 = 35.  # crfall, leaf fall period          (10 - 100)
        self.p17 = 24.2  # clma, leaf mass per area         (10 - 400) g C m-2

        self.p_vals = [self.p1, self.p2, self.p3, self.p4, self.p5, self.p5, self.p6, self.p7, self.p8,
                       self.p9, self.p10, self.p11, self.p12, self.p13, self.p14, self.p15, self.p16,
                       self.p17, self.clab, self.cf, self.cr, self.cw, self.cl, self.cs]

        self.param_dict = col.OrderedDict([('theta_min', self.p1),
                                          ('f_auto', self.p2), ('f_fol', self.p3),
                                          ('f_roo', self.p4), ('clspan', self.p5),
                                          ('theta_woo', self.p6), ('theta_roo', self.p7),
                                          ('theta_lit', self.p8), ('theta_som', self.p9),
                                          ('Theta', self.p10), ('ceff', self.p11),
                                          ('d_onset', self.p12), ('f_lab', self.p13),
                                          ('cronset', self.p14), ('d_fall', self.p15),
                                          ('crfall', self.p16), ('clma', self.p17),
                                          ('clab', self.clab), ('cf', self.cf),
                                          ('cr', self.cr), ('cw', self.cw), ('cl', self.cl),
                                          ('cs', self.cs)])
        self.pvals = np.array(self.param_dict.values())
        self.dtype_dalec = [('theta_min', '<f8'), ('f_auto', '<f8'), ('f_fol', '<f8'),
                            ('f_roo', '<f8'), ('clspan', '<f8'), ('theta_woo', '<f8'),
                            ('theta_roo', '<f8'),('theta_lit', '<f8'), ('theta_som', '<f8'),
                            ('Theta', '<f8'), ('ceff', '<f8'), ('d_onset', '<f8'),
                            ('f_lab', '<f8'), ('cronset', '<f8'), ('d_fall', '<f8'),
                            ('crfall', '<f8'), ('clma', '<f8'), ('clab', '<f8'),
                            ('cf', '<f8'), ('cr', '<f8'), ('cw', '<f8'),
                            ('cl', '<f8'), ('cs', '<f8')]
        self.xb_ew = np.array([6.28988509e-04,   4.03500221e-01,   2.71662772e-01,
                               3.49284334e-01,   1.00194789e+00,   9.39391318e-05,
                               6.72893955e-03,   1.80302080e-03,   7.59935575e-05,
                               8.00000000e-02,   5.95134533e+01,   1.50000000e+02,
                               4.77523413e-02,   2.46249836e+01,   2.91253056e+02,
                               1.50000000e+02,   2.70997155e+01,   3.79725951e+01,
                               1.00342291e+01,   1.29765563e+02,   5.84130622e+03,
                               3.99395848e+02,   2.52024132e+03])
        self.bnds = ((1e-5, 1e-2), (0.3, 0.7), (0.01, 0.5),
                     (0.01, 0.5), (1.0001, 10.), (2.5e-5, 1e-3),
                     (1e-4, 1e-2), (1e-4, 1e-2), (1e-7, 1e-3),
                     (0.018, 0.08), (10, 100), (1, 365),
                     (0.01, 0.5), (10, 100), (1, 365),
                     (10, 100), (10, 400), (10, 1000),
                     (1e-4, 1000), (10, 1000), (100, 1e5),
                     (10, 1000), (100, 2e5))
        self.bnds_tst = ((1e-5, 1e-2), (0.3, 0.7), (0.01, 0.5),
                         (0.01, 0.5), (1.0001, 10.), (2.5e-5, 1e-3),
                         (1e-4, 1e-2), (1e-4, 1e-2), (1e-7, 1e-3),
                         (0.018, 0.08), (10., 100.), (60., 150.),
                         (0.01, 0.5), (10., 100.), (220., 332.),
                         (10., 150.), (10., 400.), (10., 1000.),
                         (1e-4, 1000.), (10., 1000.), (100., 1e5),
                         (10., 1000.), (100., 2e5))

        self.dtype_dalec2 = [('f_auto', '<f8'), ('f_fol', '<f8'),
                            ('f_roo', '<f8'), ('clspan', '<f8'),
                            ('theta_woo', '<f8'), ('theta_roo', '<f8'),
                            ('theta_lit', '<f8'), ('theta_som', '<f8'),
                            ('Theta', '<f8'), ('ceff', '<f8'),
                            ('d_onset', '<f8'), ('f_lab', '<f8'),
                            ('cronset', '<f8'), ('d_fall', '<f8'),
                            ('crfall', '<f8'), ('clma', '<f8'),
                            ('clab', '<f8'), ('cf', '<f8'),
                            ('cr', '<f8'), ('cw', '<f8'), ('cl', '<f8'),
                            ('cs', '<f8')]

        self.dtype_east_west = [('theta_min', '<f8'), ('f_auto', '<f8'), ('f_fol', '<f8'),
                            ('f_roo', '<f8'), ('clspan', '<f8'),
                            ('theta_woo', '<f8'), ('theta_roo', '<f8'),
                            ('theta_lit', '<f8'), ('theta_som', '<f8'),
                            ('Theta', '<f8'),('d_onset', '<f8'), ('f_lab', '<f8'),
                            ('cronset', '<f8'), ('d_fall', '<f8'),
                            ('crfall', '<f8'), ('clma', '<f8'),
                            ('clab', '<f8'), ('cf', '<f8'),
                            ('cr', '<f8'), ('cw', '<f8'), ('cl', '<f8'),
                            ('cs', '<f8')]

        self.dtype_east_west2 = [('theta_min', '<f8'), ('f_fol', '<f8'),
                            ('f_roo', '<f8'), ('clspan', '<f8'),
                            ('theta_woo', '<f8'), ('theta_roo', '<f8'),
                            ('theta_lit', '<f8'), ('theta_som', '<f8'),
                            ('Theta', '<f8'),('d_onset', '<f8'), ('f_lab', '<f8'),
                            ('cronset', '<f8'), ('d_fall', '<f8'),
                            ('crfall', '<f8'), ('clma', '<f8'),
                            ('clab', '<f8'), ('cf', '<f8'),
                            ('cr', '<f8'), ('cw', '<f8'), ('cl', '<f8'),
                            ('cs', '<f8')]

        self.dtype_east_west3 = [('theta_min', '<f8'),
                            ('f_roo', '<f8'), ('clspan', '<f8'),
                            ('theta_woo', '<f8'), ('theta_roo', '<f8'),
                            ('theta_lit', '<f8'), ('theta_som', '<f8'),
                            ('Theta', '<f8'),('d_onset', '<f8'), ('f_lab', '<f8'),
                            ('cronset', '<f8'), ('d_fall', '<f8'),
                            ('crfall', '<f8'), ('clma', '<f8'),
                            ('clab', '<f8'), ('cf', '<f8'),
                            ('cr', '<f8'), ('cw', '<f8'), ('cl', '<f8'),
                            ('cs', '<f8')]

        self.ah_pvals = np.array([9.41e-04, 4.7e-01, 2.8e-01, 2.60e-01, 1.01e+00, 2.6e-04,
                                  2.48e-03, 3.38e-03, 2.6e-06, 1.93e-02, 9.0e+01, 1.4e+02,
                                  4.629e-01, 2.7e+01, 3.08e+02, 3.5e+01, 5.2e+01, 78.,
                                  2., 134., 14257.32, 68.95, 18625.77])

        self.edinburgh_median = np.array([2.29180076e-04,   5.31804031e-01,   6.69448981e-02,
                                          4.46049258e-01,   1.18143120e+00,   5.31584216e-05,
                                          2.25487423e-03,   2.44782152e-03,   7.71092378e-05,
                                          3.82591095e-02,   7.47751776e+01,   1.16238252e+02,
                                          3.26252225e-01,   4.18554035e+01,   2.27257813e+02,
                                          1.20915004e+02,   1.15533213e+02,   1.27804720e+02,
                                          6.02259491e+01,   2.09997016e+02,   4.22672530e+03,
                                          3.67801053e+02,   1.62565304e+03])

        self.edinburgh_mean = np.array([(9.80983217e-04,   5.19025559e-01,   1.08612889e-01,
                                        4.84356048e-01,   1.19950434e+00,   1.01336503e-04,
                                        3.22465935e-03,   3.44239452e-03,   1.11320287e-04,
                                        4.14726183e-02,   7.14355778e+01,   1.15778224e+02,
                                        3.20361827e-01,   4.13391057e+01,   2.20529309e+02,
                                        1.16768714e+02,   1.28460812e+02,   1.36541509e+02,
                                        6.86396830e+01,   2.83782534e+02,   6.50600814e+03,
                                        5.98832031e+02,   1.93625350e+03)], dtype=self.dtype_dalec)

        self.e_w_lai = np.array([(5.47099980e-04,   4.49228870e-01,   4.09088649e-02,
                                      3.69550552e-01,   1.08939445e+00,   1.01196932e-04,
                                      5.41063421e-03,   4.38664659e-03,   1.31134514e-04,
                                      9.35399032e-02,   1.58422811e+02,
                                      7.92729408e-02,   1.89092607e+01,   3.04928379e+02,
                                      5.44663654e+01,   2.92938186e+01,   7.30929677e+01,
                                      1.31261269e+01,   2.10259000e+02,   7.18220358e+03,
                                      1.69702341e+02,   1.95033734e+03)], dtype=self.dtype_east_west)

        self.e_w_edc = np.array([(1.23000000e-05,   3.13358350e-01,   3.00629189e-01,
                                4.45265166e-01,   1.02310470e+00,   1.22836138e-04,
                                5.04088931e-03,   1.56202990e-03,   1.48252124e-04,
                                7.61636968e-02,   1.22954168e+02,
                                1.11000000e-02,   4.67979617e+01,   2.87147216e+02,
                                5.51760150e+01,   5.16317404e+01,   1.45000000e+01,
                                1.43000000e+01,   5.01480167e+02,   7.26249788e+03,
                                6.26033838e+02,   2.35514838e+03)], dtype=self.dtype_east_west)

        self.e_w_assim2 = np.array([(9.80983217e-04, 1.08612889e-01,
                                    4.84356048e-01,   1.19950434e+00,   1.01336503e-04,
                                    3.22465935e-03,   3.44239452e-03,   1.11320287e-04,
                                    4.14726183e-02,   1.15778224e+02,
                                    3.20361827e-01,   4.13391057e+01,   2.20529309e+02,
                                    1.16768714e+02,   1.28460812e+02,   1.36541509e+02,
                                    6.86396830e+01,   2.83782534e+02,   6.50600814e+03,
                                    5.98832031e+02,   1.93625350e+03)], dtype=self.dtype_east_west2)

        self.e_w_edc2 = np.array([(1.23000000e-05, 3.00629189e-01,
                                   4.45265166e-01,   1.02310470e+00,   1.22836138e-04,
                                   5.04088931e-03,   1.56202990e-03,   1.48252124e-04,
                                   7.61636968e-02,   1.22954168e+02,
                                   1.11000000e-02,   4.67979617e+01,   2.87147216e+02,
                                   5.51760150e+01,   5.16317404e+01,   1.45000000e+01,
                                   1.43000000e+01,   5.01480167e+02,   7.26249788e+03,
                                   6.26033838e+02,   2.35514838e+03)], dtype=self.dtype_east_west2)

        self.xb_ew2 = np.array([(6.28988509e-04, 2.71662772e-01,
                               3.49284334e-01,   1.00194789e+00,   9.39391318e-05,
                               6.72893955e-03,   1.80302080e-03,   7.59935575e-05,
                               8.00000000e-02,   5.95134533e+01,
                               4.77523413e-02,   2.46249836e+01,   2.91253056e+02,
                               1.50000000e+02,   2.70997155e+01,   3.79725951e+01,
                               1.00342291e+01,   1.29765563e+02,   5.84130622e+03,
                               3.99395848e+02,   2.52024132e+03)], dtype=self.dtype_east_west2)

        self.e_w_assim3 = np.array([(9.80983217e-04,
                                    4.84356048e-01,   1.19950434e+00,   1.01336503e-04,
                                    3.22465935e-03,   3.44239452e-03,   1.11320287e-04,
                                    4.14726183e-02,   1.15778224e+02,
                                    3.20361827e-01,   4.13391057e+01,   2.20529309e+02,
                                    1.16768714e+02,   1.28460812e+02,   1.36541509e+02,
                                    6.86396830e+01,   2.83782534e+02,   6.50600814e+03,
                                    5.98832031e+02,   1.93625350e+03)], dtype=self.dtype_east_west3)

        self.e_w_edc3 = np.array([(1.23000000e-05,
                                   4.45265166e-01,   1.02310470e+00,   1.22836138e-04,
                                   5.04088931e-03,   1.56202990e-03,   1.48252124e-04,
                                   7.61636968e-02,   1.22954168e+02,
                                   1.11000000e-02,   4.67979617e+01,   2.87147216e+02,
                                   5.51760150e+01,   5.16317404e+01,   1.45000000e+01,
                                   1.43000000e+01,   5.01480167e+02,   7.26249788e+03,
                                   6.26033838e+02,   2.35514838e+03)], dtype=self.dtype_east_west3)

        self.edinburgh_std = np.array([2.03001590e-03,   1.16829160e-01,   1.11585876e-01,
                                       2.98860194e-01,   1.16141739e-01,   1.36472702e-04,
                                       2.92998472e-03,   3.11712858e-03,   1.18105073e-04,
                                       1.62308654e-02,   2.04219069e+01,   6.25696097e+00,
                                       1.14535431e-01,   1.40482247e+01,   3.72380005e+01,
                                       2.25938092e+01,   6.41030587e+01,   6.62621885e+01,
                                       3.59002726e+01,   2.19315727e+02,   7.14323513e+03,
                                       5.45013287e+02,   1.27646316e+03])

        self.xa_edc = np.array([1.23000000e-05,   3.13358350e-01,   3.00629189e-01,
                                4.45265166e-01,   1.02310470e+00,   1.22836138e-04,
                                5.04088931e-03,   1.56202990e-03,   1.48252124e-04,
                                7.61636968e-02,   9.27591545e+01,   1.22954168e+02,
                                1.11000000e-02,   4.67979617e+01,   2.87147216e+02,
                                5.51760150e+01,   5.16317404e+01,   1.45000000e+01,
                                1.43000000e+01,   5.01480167e+02,   7.26249788e+03,
                                6.26033838e+02,   2.35514838e+03])

        self.xb = self.pvals
        self.xb_ew = np.array([6.28988509e-04,   4.03500221e-01,   2.71662772e-01,
                               3.49284334e-01,   1.00194789e+00,   9.39391318e-05,
                               6.72893955e-03,   1.80302080e-03,   7.59935575e-05,
                               8.00000000e-02,   5.95134533e+01,   1.50000000e+02,
                               4.77523413e-02,   2.46249836e+01,   2.91253056e+02,
                               1.50000000e+02,   2.70997155e+01,   3.79725951e+01,
                               1.00342291e+01,   1.29765563e+02,   5.84130622e+03,
                               3.99395848e+02,   2.52024132e+03])
        # self.B = self.make_b(self.edinburgh_std)
        self.B = pickle.load(open('b_edc.p', 'r'))

        self.xa = None

        # Daily temperatures degC
        self.t_mean = data.variables['daily_mean_temp'][self.start_idx:self.end_idx, 0, 0]
        self.t_max = data.variables['daily_max_temp'][self.start_idx:self.end_idx, 0, 0]
        self.t_min = data.variables['daily_min_temp'][self.start_idx:self.end_idx, 0, 0]
        self.t_range = np.array(self.t_max) - np.array(self.t_min)
        self.t_day = data.variables['mean_temp_day'][self.start_idx:self.end_idx, 0, 0]
        self.t_night = data.variables['mean_temp_night'][self.start_idx:self.end_idx, 0, 0]
        # Daily soil temperatures at 3 cm depth
        self.t_mean_soil = data.variables['daily_mean_soil_temp'][self.start_idx:self.end_idx, 0, 0]
        self.t_day_soil = data.variables['mean_soil_temp_day'][self.start_idx:self.end_idx, 0, 0]
        self.t_night_soil = data.variables['mean_soil_temp_night'][self.start_idx:self.end_idx, 0, 0]

        # Driving Data
        self.I = data.variables['rg'][self.start_idx:self.end_idx, 0, 0]  # incident radiation
        self.ca = 390.0  # atmospheric carbon
        self.D = data.variables['doy'][self.start_idx:self.end_idx, 0, 0]  # day of year
        self.day_len = data.variables['day_length'][self.start_idx:self.end_idx, 0, 0]
        self.night_len = data.variables['night_length'][self.start_idx:self.end_idx, 0, 0]
        self.year = np.array([date.year for date in self.dates])  # year
        self.month = np.array([date.month for date in self.dates])  # month
        self.date = np.array([date.day for date in self.dates])  # date in month

        # Constants for ACM model
        self.acm_williams_xls = np.array([0.0155, 1.526, 324.1, 0.2017,
                                          1.315, 2.595, 0.037, 0.2268,
                                          0.9576])
        self.acm_reflex = np.array([0.0156935, 4.22273, 208.868, 0.0453194,
                                   0.37836, 7.19298, 0.011136, 2.1001,
                                   0.789798])
        self.acm = self.acm_reflex  # (currently using params from REFLEX)
        self.phi_d = -2.5  # max. soil leaf water potential difference
        self.R_tot = 1.  # total plant-soil hydrolic resistance
        self.lat = 0.89133965  # latitude of forest site in radians
        #                        lat = 51.153525 deg, lon = -0.858352 deg

        # misc
        self.ca = 390.0  # atmospheric carbon
        self.radconv = 365.25 / np.pi
        self.delta_t = delta_t

        # Background standard deviations for carbon pools & B matrix
        self.sigb_clab = 7.5  # 20%
        self.sigb_cf = 10.0  # 20%
        self.sigb_cw = 1000.  # 20%
        self.sigb_cr = 13.5  # 20%
        self.sigb_cl = 7.0  # 20%
        self.sigb_cs = 1500.  # 20%

        # Observation standard deviations for carbon pools and NEE
        self.sigo_clab = 7.5  # 10%
        self.sigo_cf = 10.0  # 10%
        self.sigo_cw = 200.  # from AH dbh measurements
        self.sigo_cr = 8.5  # from AH dbh measurements
        self.sigo_cl = 7.0  # 30%
        self.sigo_cs = 1500.0  # 30%
        self.sigo_nee = 0.71  # g C m-2 day-1
        self.sigo_nee_day = 0.71
        self.sigo_nee_night = 2*0.71
        self.sigo_lf = 0.5
        self.sigo_lw = 0.5
        self.sigo_litresp = 0.5
        self.sigo_soilresp = 0.6
        self.sigo_rtot = 0.71
        self.sigo_rh = 0.6
        self.sigo_lai = 0.5  # from AH optical measurements
        self.sigo_clma = 5.0  # from AH litter scans
        self.sigo_donset = 5.
        self.sigo_dfall = 7.

        self.error_dict = {'clab': self.sigo_clab, 'cf': self.sigo_cf, 'c_woo': self.sigo_cw,
                           'cl': self.sigo_cl, 'c_roo': self.sigo_cr, 'cs': self.sigo_cs,
                           'nee': self.sigo_nee, 'nee_day': self.sigo_nee_day, 'nee_night': self.sigo_nee_night,
                           'lf': self.sigo_lf, 'lw': self.sigo_lw, 'litresp': self.sigo_litresp,
                           'soilresp': self.sigo_soilresp, 'rtot': self.sigo_rtot, 'rh': self.sigo_rh,
                           'lai': self.sigo_lai, 'clma': self.sigo_clma, 'd_onset': self.sigo_donset,
                           'd_fall': self.sigo_dfall}
        self.possible_obs = ['gpp', 'lf', 'lw', 'rt', 'nee', 'nee_east', 'nee_west', 'nee_day', 'nee_day_east',
                             'nee_day_west', 'nee_night', 'nee_night_east', 'nee_night_west', 'cf', 'cl',
                             'c_roo', 'c_roo_east', 'c_roo_west', 'c_woo', 'c_woo_east', 'c_woo_west', 'cs', 'lai',
                             'lai_east', 'lai_west', 'clma', 'clab', 'litresp', 'soilresp','rtot', 'rh', 'rabg',
                             'd_onset', 'd_fall', 'groundresp']

        # Extract observations for assimilation
        self.ob_dict, self.ob_err_dict = self.assimilation_obs(ob_str, data)
        data.close()

    def assimilation_obs(self, ob_str, data):
        """ Extracts observations and errors for assimilation into dictionaries
        :param obs_str: string of observations separated by commas
        :return: dictionary of observation values, dictionary of corresponding observation errors
        """
        obs_lst = re.findall(r'[^,;\s]+', ob_str)
        obs_dict = {}
        obs_err_dict = {}
        for ob in obs_lst:
            if ob not in self.possible_obs:
                raise Exception('Invalid observations entered, please check \
                                 function input')
            elif ob in ['nee_east', 'nee_day_east', 'nee_night_east']:
                ob_del = ob[:-5]
                obs = data.variables[ob_del][self.start_idx:self.end_idx, 0, 0]
                obs_dict[ob_del] = obs
                obs_err_dict[ob_del] = (obs/obs) * self.error_dict[ob_del]
                origin = data.variables[ob_del+'_origin'][self.start_idx:self.end_idx, 0, 0]
                idx = np.where(origin == 2)
                idx2 = np.where(origin == 0)
                for x in np.concatenate((idx[0], idx2[0]), axis=0):
                    obs_dict[ob_del][x] = float('NaN')
                    obs_err_dict[ob_del][x] = float('NaN')
            elif ob in ['nee_west', 'nee_day_west', 'nee_night_west']:
                ob_del = ob[:-5]
                obs = data.variables[ob_del][self.start_idx:self.end_idx, 0, 0]
                obs_dict[ob_del] = obs
                obs_err_dict[ob_del] = (obs/obs) * self.error_dict[ob_del]
                origin = data.variables[ob_del+'_origin'][self.start_idx:self.end_idx, 0, 0]
                idx = np.where(origin < 2)
                for x in idx[0]:
                    obs_dict[ob_del][x] = float('NaN')
                    obs_err_dict[ob_del][x] = float('NaN')
            elif ob in ['c_woo_east', 'c_woo_west', 'lai_east', 'lai_west', 'c_roo_east', 'c_roo_west']:
                ob_del = ob[:-5]
                obs = data.variables[ob][self.start_idx:self.end_idx, 0, 0]
                obs_dict[ob_del] = obs
                obs_err_dict[ob_del] = (obs/obs) * self.error_dict[ob_del]
            else:
                obs = data.variables[ob][self.start_idx:self.end_idx, 0, 0]
                obs_dict[ob] = obs
                obs_err_dict[ob] = (obs/obs) * self.error_dict[ob]

        return obs_dict, obs_err_dict

    @staticmethod
    def find_date_idx(date, data):
        """ Finds index in netcdf file for given date
        :param date: date in format specified in DalecData class
        :return: date index
        """
        if type(date) == int:
            d_time = dt.datetime(date, 1, 1)
        elif type(date) == tuple:
            d_time = dt.datetime(date[0], date[1], date[2])
        else:
            raise ValueError('Date wrong format, please check input')
        times = data.variables['time']
        return nC.date2index(d_time, times, select='nearest')

    @staticmethod
    def make_b(b_std):
        """ Creates diagonal B matrix.
        :param b_std: array of standard deviations corresponding to each model parameter
        :return: 23 x 23 diagonal background error covariance matrix
        """
        b_mat = b_std**2 * np.eye(23)
        return b_mat


class DalecDataTwin(DalecData):
    def __init__(self, start_date, end_date, ob_str, err_scale=0.25,
                 nc_file='../../alice_holt_data/ah_data_daily_test.nc'):

        DalecData.__init__(self, start_date, end_date, ob_str, nc_file)
        self.m = mc.DalecModel(self)

        # Define truth and background
        self.x_truth = self.edinburgh_median
        self.st_dev = 0.1*self.x_truth
        # self.B = self.make_b(self.st_dev)

        # Make EDC B
        b_cor = pickle.load(open('b_edc_cor.p', 'r'))  # load correlation matrix from 2016 paper
        b_std = self.make_b(np.sqrt(self.st_dev))  # Create diagonal matrix of standard deviations
        self.B = np.dot(np.dot(b_std, b_cor), b_std)  # Create correlated B
        # self.xb = self.random_pert(self.random_pert(self.x_truth))
        self.xb = np.array([2.53533992e-04,   5.85073161e-01,   7.43127332e-02,
                            4.99707798e-01,   1.38993876e+00,   6.11913792e-05,
                            2.58484324e-03,   2.79379720e-03,   8.72422101e-05,
                            4.35144260e-02,   8.73669864e+01,   1.29813051e+02,
                            3.87867223e-01,   4.69894281e+01,   2.78080852e+02,
                            9.15080347e+01,   1.36269157e+02,   1.44176657e+02,
                            6.71153814e+01,   2.42199267e+02,   4.96249386e+03,
                            4.15128028e+02,   1.90797697e+03])
        # Extract observations for assimilation
        self.ob_dict, self.ob_err_dict = self.create_twin_data(ob_str, err_scale)

    def create_twin_data(self, ob_str, err_scale=0.25):
        """ Creates a set of twin modelled observations corresponding to the same positions as the true observations
        :param ob_str: str of observations
        :param err_scale: factor by which to scale observation error and added gaussian noise
        :return: observation dictionary, observation error dictionary
        """
        obs_lst = re.findall(r'[^,;\s]+', ob_str)
        obs_dict = {}
        obs_err_dict = {}
        mod_lst = self.m.mod_list(self.x_truth)
        for ob in obs_lst:
            if ob not in self.possible_obs:
                raise Exception('Invalid observations entered, please check \
                                 function input')
            else:
                obs = self.ob_dict[ob]  # actual observations
                mod_obs = (obs/obs) * self.m.oblist(ob, mod_lst)  # modelled observation corresponding to same
                # position as actual obs
                # adding error to modelled observations
                mod_ob_assim = np.array([mod_ob + random.gauss(0, err_scale*self.error_dict[ob])
                                         for mod_ob in mod_obs])
                obs_dict[ob] = mod_ob_assim
                if err_scale == 0.0:
                    obs_err_dict[ob] = 0.25*self.ob_err_dict[ob]
                else:
                    obs_err_dict[ob] = err_scale*self.ob_err_dict[ob]
        return obs_dict, obs_err_dict

    def random_pert(self, pvals):
        """ Perturbs parameter values with given standard deviation
        :param pvals: parameter values to perturb
        :return: perturbed parameters
        """
        pval_approx = np.ones(23)*-9999.
        x = 0
        for p in pvals:
            pval_approx[x] = p + random.gauss(0, self.st_dev[x])
            if self.bnds[x][1] < pval_approx[x]:
                pval_approx[x] = self.bnds[x][1] - abs(random.gauss(0, self.bnds[x][1]*0.001))
            elif self.bnds[x][0] > pval_approx[x]:
                pval_approx[x] = self.bnds[x][0] + abs(random.gauss(0, self.bnds[x][0]*0.001))

            x += 1

        return pval_approx

    def random_pert_uniform(self, pvals):
        """ Perturbs parameter values with given standard deviation
        :param pvals: parameter values to perturb
        :return: perturbed parameters
        """
        pval_approx = np.ones(23)*-9999.
        xt = self.x_truth
        x = 0
        for p in pvals:
            pval_approx[x] = p + random.uniform(-0.1*xt[x], 0.1*xt[x])
            if 0.3 < abs(pval_approx[x] - self.x_truth[x])/self.x_truth[x]:
                while 0.3 < abs(pval_approx[x] - self.x_truth[x])/self.x_truth[x]:
                    pval_approx[x] = pval_approx[x] - abs(random.uniform(-0.1*xt[x], 0.1*xt[x]))
            if abs(pval_approx[x] - self.x_truth[x])/self.x_truth[x] < 0.12:
                while abs(pval_approx[x] - self.x_truth[x])/self.x_truth[x] < 0.12:
                    pval_approx[x] = pval_approx[x] + abs(random.uniform(-0.1*xt[x], 0.1*xt[x]))
            if self.bnds[x][1] < pval_approx[x]:
                pval_approx[x] = self.bnds[x][1] - abs(random.gauss(0, self.bnds[x][1]*0.001))
            elif self.bnds[x][0] > pval_approx[x]:
                pval_approx[x] = self.bnds[x][0] + abs(random.gauss(0, self.bnds[x][0]*0.001))
            x += 1
        return pval_approx

    def test_pvals(self, pvals):
        """ Test if a parameter set falls within the bounds or not
        :param pvals: parameter values to test
        :return: pvals
        """
        x = 0
        for bnd in self.bnds:
            if bnd[0] < pvals[x] < bnd[1]:
                print '%x in bnds' %x
            else:
                print '%x not in bnds' %x
            x += 1
        return pvals


