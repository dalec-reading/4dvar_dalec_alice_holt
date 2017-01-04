import numpy as np
import netCDF4 as nC
import matplotlib.mlab as ml
import datetime as dt
import ephem


def open_csv(filename, missing_val='N/A'):
    """Opens a csv file into a recorded array.
    """
    return ml.csv2rec(filename, missing=missing_val)


def open_netcdf(filename):
    """Opens a netCDF file
    """
    return nC.Dataset(filename, 'a')


def ah_str2date(date_str):
    """Converts string into datetime object for alice holt spreadsheet.
    """
    return dt.datetime.strptime(date_str, "%d/%m/%Y %H:%M")


def add_data2nc(nc_file, csv_file, data_title, nc_title, date_col='date_combined'):
    """ Adds data to a netCDF file
    :param nc_file: netCDF file location str
    :param csv_file: csv file location str
    :param data_title: title column for data to add as str
    :param nc_title: title of nc variable to add it to as str
    :param date_col: title of date column as str
    :return: nothing
    """
    nc_ob = open_netcdf(nc_file)
    var = nc_ob.variables[nc_title]
    times = nc_ob.variables['time']
    dat_arr = open_csv(csv_file)
    for x in xrange(len(dat_arr[date_col])):
        try:
            idx = nC.date2index(dat_arr[date_col][x], times)
        except ValueError:
            print x
        var[idx, 0, 0] = dat_arr[data_title][x]
    nc_ob.close()
    return 'data updated!'


def nc_create_var(nc_file, var_name, dims, dtype=np.float16):
    """ Adds new variable to netCDF data file
    :param nc_file: netcdf data file location as str
    :param var_name: name for new variable as str
    :param dims: dimensions of new variable, tuple containing dimension strings
    :param dtype: data type of new variable
    :return:
    """
    nc_ob = open_netcdf(nc_file)
    nc_ob.createVariable(var_name, dtype, dims)
    nc_ob.close()
    return 'variable added!'


# event_time is just a date time corresponding to an sql timestamp
def type_of_light(latitude, longitude, event_time, utc_time, horizon):

    o = ephem.Observer()
    o.lat, o.long, o.date, o.horizon = latitude, longitude, event_time, horizon

    # print "event date ", o.date

    # print "prev rising: ", o.previous_rising(ephem.Sun())
    # print "prev setting: ", o.previous_setting(ephem.Sun())
    # print "next rise: ", o.next_rising(ephem.Sun())
    # print "next set: ", o.next_setting(ephem.Sun())

    if o.previous_rising(ephem.Sun()) > o.previous_setting(ephem.Sun()):
        return "day"
    else:
        return "night"


def ah_day_night(event_time, horizon='0'):
    """calculate day and night for alice holt straits
    """
    time_str = event_time.strftime("%Y/%m/%d %H:%M")
    light = type_of_light('51.153526', '0.858348', time_str,'0', horizon)
    if light == "day":
        return 1
    else:
        return 0


def fill_is_day(time, is_day):
    time_lst = nC.num2date(time[:], time.units)
    for t in enumerate(time_lst):
        is_day[t[0], 0, 0] = ah_day_night(t[1])
    return 'yay'


def find_date(nc_time_int, nc_time_ob):
    """ Returns the date as a datetime given any datetime object
    :param nc_time_int: date as an integer corresponding to units in nc file
    :param nc_time_ob: netcdf time object
    :return: date as a datetime object
    """
    day = nC.num2date(nc_time_int, nc_time_int.units).date()
    return dt.datetime.combine(day, dt.datetime.min.time())


def create_date_list(start_date, end_date, del_t='day'):
    """ Creates a list of daily or yearly datetime objects
    :param start_date: start date for list as datetime
    :param end_date: end date for list as datetime
    :param del_t: time step for date list as string, can be 'day', 'year', or 'half_hour'
    :return: datetime list
    """
    times = []
    if del_t == 'day':
        delta = dt.timedelta(hours=24)
    elif del_t == 'year':
        delta = dt.timedelta(years=1)
    elif del_t == 'half_hour':
        delta = dt.timedelta(minutes=30)
    date = start_date
    while date <= end_date:
        times.append(date)
        date = date + delta
    return times


def nc_doy(doy, times):
    """ Finds day of year for netcdf time dimension
    :param doy: day of year netcdf variable
    :param times: half hourly times as netcdf variable
    :return:
    """
    for t in times[:]:
        day_time = nC.num2date(t, times.units)
        idx = nC.date2index(day_time, times)
        tt = day_time.timetuple()
        doy[idx] = tt.tm_yday
    return 'yay'


def nc_mean_daily_temp(half_hourly_temps, daily_mean_temps, times, time_lst):
    """ Finds mean daily temperatures from half hourly temperature data
    :param half_hourly_temps: half hourly temperatures as netcdf variable
    :param daily_mean_temps: empty daily mean temperature netcdf variable
    :param times: half hourly times as netcdf variable
    :param time_lst: list of daily datetime objects
    :return:
    """
    for t in enumerate(time_lst):
        idx = nC.date2index(t[1], times)
        mean_daily_temp = np.mean(half_hourly_temps[idx:idx+48, 0, 0])
        daily_mean_temps[t[0], 0, 0] = mean_daily_temp
    return 'yay'


def nc_max_daily_temp(half_hourly_temps, daily_max_temps, times, time_lst):
    for t in time_lst:
        idx = nC.date2index(t, times)
        max_daily_temp = np.max(half_hourly_temps[idx:idx+48, 0, 0])
        daily_max_temps[idx, 0, 0] = max_daily_temp
    return 'yay'


def nc_min_daily_temp(half_hourly_temps, daily_min_temps, times, time_lst):
    for t in time_lst:
        idx = nC.date2index(t, times)
        min_daily_temp = np.min(half_hourly_temps[idx:idx+48, 0, 0])
        daily_min_temps[idx, 0, 0] = min_daily_temp
    return 'yay'


def nc_total_daily_rg(half_hourly_rg, daily_rg, times, time_lst):
    """ Finds total daily global radiation from half hourly global radiation data
    :param half_hourly_rg: half hourly global radiation as netcdf variable (W m-2)
    :param daily_rg: empty total daily global radiation netcdf variable (M J m-2 day-1)
    :param times: half hourly times as netcdf variable
    :param time_lst: list of daily datetime objects
    :return:
    """
    for t in time_lst:
        idx = nC.date2index(t, times)
        total_daily_rg = 30*60*1e-6*np.sum(half_hourly_rg[idx:idx+48, 0, 0])  # Convert W m-2 to M J m-2 day-1
        daily_rg[idx, 0, 0] = total_daily_rg
    return 'yay'


def nc_day_len(is_day, day_len, times, time_lst):
    """ Finds total daily global radiation from half hourly global radiation data
    :param is_day: half hourly 'is day' netcdf variable with values of 1 day or 0 night
    :param day_len: empty day length netcdf variable (hours)
    :param times: half hourly times as netcdf variable
    :param time_lst: list of daily datetime objects
    :return:
    """
    for t in enumerate(time_lst):
        idx = nC.date2index(t[1], times)
        where_day = np.where(is_day[idx:idx+48, 0, 0] == 1)[0]
        day_len[t[0], 0, 0] = len(where_day)*0.5
    return 'yay'


def nc_night_len(is_day, night_len, times, time_lst):
    """ Finds total daily global radiation from half hourly global radiation data
    :param is_day: half hourly 'is day' netcdf variable with values of 1 day or 0 night
    :param night_len: empty day length netcdf variable (hours)
    :param times: half hourly times as netcdf variable
    :param time_lst: list of daily datetime objects
    :return:
    """
    for t in enumerate(time_lst):
        idx = nC.date2index(t[1], times)
        night_idx = idx + np.where(is_day[idx:idx+48, 0, 0] == 1)[0][-1] + 1
        night_length = 0
        while is_day[night_idx, 0, 0] == 0:
            night_length += 1  # if final night incomplete will return false night length for last day of data
            night_idx += 1
            if night_idx >= len(times):
                night_length = float('NaN')
                break
        night_len[t[0], 0, 0] = night_length*0.5
    return 'yay'


def nc_day_mean_temp(is_day, hh_temp, mean_t_day, times, time_lst):
    """ Finds total daily global radiation from half hourly global radiation data
    :param is_day: half hourly 'is day' netcdf variable with values of 1 day or 0 night
    :param hh_temp: half hourly temperatures netcdf variable
    :param mean_t_day: netcdf variable to fill with mean daytime temperatures
    :param times: half hourly times as netcdf variable
    :param time_lst: list of daily datetime objects
    :return:
    """
    for t in enumerate(time_lst):
        idx = nC.date2index(t[1], times)
        where_day = np.where(is_day[idx:idx+48, 0, 0] == 1)[0]
        mean_t_day[t[0], 0, 0] = np.mean(hh_temp[idx+where_day[0]:idx+where_day[-1]])
    return 'yay'


def nc_night_mean_temp(is_day, hh_temp, mean_t_night, times, time_lst):
    """ Finds total daily global radiation from half hourly global radiation data
    :param is_day: half hourly 'is day' netcdf variable with values of 1 day or 0 night
    :param hh_temp: half hourly temperatures netcdf variable
    :param mean_t_night: netcdf variable to fill with mean nighttime temperatures
    :param times: half hourly times as netcdf variable
    :param time_lst: list of daily datetime objects
    :return:
    """
    for t in enumerate(time_lst):
        idx = nC.date2index(t[1], times)
        night_idx1 = idx + np.where(is_day[idx:idx+48, 0, 0] == 1)[0][-1] + 1
        night_idx2 = night_idx1
        while is_day[night_idx2, 0, 0] == 0:
            night_idx2 += 1
            if night_idx2 >= len(times):
                night_idx2 = float('NaN')
                break
        if np.isnan(night_idx2) == True:
            mean_t_night[t[0], 0, 0] = float('NaN')
            break
        mean_t_night[t[0], 0, 0] = np.mean(hh_temp[night_idx1:night_idx2, 0, 0])
    return 'yay'


def clip_co2_flux(co2_flux, clipped_co2_flux, above_clip=50., below_clip=-65.):
    """ Clips co2_flux data and saves in clipped_co2_flux netcdf variable
    :param co2_flux: co2 flux measurements from tower as netcdf variable
    :param clipped_co2_flux: netcdf variable to save clipped measurements in
    :return:
    """
    clipped_co2_flux[:, 0, 0] = co2_flux[:, 0, 0]
    clip_nee = clipped_co2_flux[:, 0, 0]
    # Value of +/- 70 u mol m-2 s-1 chosen for upper and lower limit to clip.
    clip_nee[clip_nee > above_clip] = float('NaN')
    clip_nee[clip_nee < below_clip] = float('NaN')
    clipped_co2_flux[:, 0, 0] = clip_nee
    return 'yay'


def clip_co2_flux_u_star(co2_flux, clipped_co2_flux, u_star, is_day):
    """ Clips co2_flux data based on u_star value and saves in clipped_co2_flux netcdf variable
    :param co2_flux: co2 flux measurements from tower as netcdf variable
    :param clipped_co2_flux: netcdf variable to save clipped measurements in
    :return:
    """
    clipped_co2_flux[:, 0, 0] = co2_flux[:, 0, 0]
    clip_nee = clipped_co2_flux[:, 0, 0]
    # Value of +/- 70 u mol m-2 s-1 chosen for upper and lower limit to clip.
    for x in xrange(len(clip_nee)):
        if is_day[x, 0, 0] == 1:
            if u_star[x, 0, 0] < 0.2:
                clip_nee[x] = float('NaN')
        elif is_day[x, 0, 0] == 0:
            if u_star[x, 0, 0] < 0.5:
                clip_nee[x] = float('NaN')
    clipped_co2_flux[:, 0, 0] = clip_nee
    return 'yay'


def clip_co2_flux_u_star2(co2_flux, clipped_co2_flux, u_star):
    """ Clips co2_flux data based on u_star value and saves in clipped_co2_flux netcdf variable
    :param co2_flux: co2 flux measurements from tower as netcdf variable
    :param clipped_co2_flux: netcdf variable to save clipped measurements in
    :return:
    """
    clipped_co2_flux[:, 0, 0] = co2_flux[:, 0, 0]
    clip_nee = clipped_co2_flux[:, 0, 0]
    # Value of +/- 70 u mol m-2 s-1 chosen for upper and lower limit to clip.
    for x in xrange(len(clip_nee)):
        if u_star[x, 0, 0] < 0.2:
            clip_nee[x] = float('NaN')
    clipped_co2_flux[:, 0, 0] = clip_nee
    return 'yay'


def clip_co2_flux_mean(co2_flux, clipped_co2_flux, idx1, idx2):
    """ Clips co2_flux data and saves in clipped_co2_flux netcdf variable
    :param co2_flux: co2 flux measurements from tower as netcdf variable
    :param clipped_co2_flux: netcdf variable to save clipped measurements in
    :return:
    """
    clipped_co2_flux[idx1:idx2, 0, 0] = co2_flux[idx1:idx2, 0, 0]
    clip_nee = clipped_co2_flux[idx1:idx2, 0, 0]
    # Clips NEE obs for +/- observations by 3 standard deviations
    clip_mean_pos = np.nanmean(clip_nee[clip_nee > 0])
    clip_std_pos = np.nanstd(clip_nee[clip_nee > 0])
    clip_mean_neg = np.nanmean(clip_nee[clip_nee < 0])
    clip_std_neg = np.nanstd(clip_nee[clip_nee < 0])
    clip_nee[clip_nee > clip_mean_pos+3*clip_std_pos] = float('NaN')
    clip_nee[clip_nee < clip_mean_neg-3*clip_std_neg] = float('NaN')
    clipped_co2_flux[idx1:idx2, 0, 0] = clip_nee
    return 'yay'


def find_indices_year(times, year):
    """ Returns the first and last time index for a given year
    :param times: time netcdf variable
    :param year: year to find indices for as an integer
    :return: first index, final index
    """
    year_entries = [x for x in times[:] if nC.num2date(x, times.units).year == year]
    idx1 = np.where(times[:] == year_entries[0])[0][0]
    idx2 = np.where(times[:] == year_entries[-1])[0][0]
    return idx1, idx2


def find_indices_month(time_arr, month, time_units):
    """ Returns the first and last time index for a given month
    :param time_arr: slice from a netcdf time variable
    :param month: month to find indices for as an integer
    :param time_units: units for time variable
    :return: first index, final index
    """
    month_entries = [x for x in time_arr if nC.num2date(x, time_units).month == month]
    idx1 = np.where(time_arr == month_entries[0])[0][0]
    idx2 = np.where(time_arr == month_entries[-1])[0][0]
    return idx1, idx2


def make_year_lst(times):
    start_yr = nC.num2date(times[0], times.units).year
    end_yr = nC.num2date(times[-1], times.units).year
    year_lst = np.arange(start_yr, end_yr+1)
    return year_lst


def clip_co2_flux_wrapper(co2_flux, clipped_co2_flux, times):
    """ Clips co2 flux data by year and month using clip_co2_flux_mean fn
    :param co2_flux: co2 flux netcdf variable
    :param clipped_co2_flux: netcdf variable to put clipped flux observations
    :param times: time netcdf variable
    :return:
    """
    yr_lst = make_year_lst(times)
    for yr in yr_lst:
        yr_idx1, yr_idx2 = find_indices_year(times, yr)
        clip_co2_flux_mean(co2_flux, clipped_co2_flux, yr_idx1, yr_idx2)
        times_yr = times[yr_idx1:yr_idx2]
        for month in np.arange(1, 13):
            month_idx1, month_idx2 = find_indices_month(times_yr, month, times.units)
            clip_co2_flux_mean(co2_flux, clipped_co2_flux, yr_idx1+month_idx1, yr_idx1+month_idx2)
    return 'yay'


def quality_control_co2_flux(clipped_co2_flux, qc_co2_flux, nee, nee_std, idx1, idx2, idx, qc_tol=3):
    """ Quality controls flux data and processes it to daily or half daily
    :param clipped_co2_flux: half hourly clipped co2 flux as netcdf variable
    :param qc_co2_flux: qc flags as nc variable corresponding to half hourly co2 flux
    :param nee: nc variable to fill with processed data
    :param nee: nc variable to fill with processed data standard deviation
    :param idx1: start index to being quality control and data averaging
    :param idx: start of day index
    :return:
    """
    clip_nee = clipped_co2_flux[idx1:idx2, 0, 0]
    clip_mean_pos = np.nanmean(clip_nee[clip_nee > 0])
    clip_std_pos = np.nanstd(clip_nee[clip_nee > 0])
    clip_mean_neg = np.nanmean(clip_nee[clip_nee < 0])
    clip_std_neg = np.nanstd(clip_nee[clip_nee < 0])
    fill = 0
    qc_flag1 = 0
    qc_flag2 = 0

    for x in xrange(idx1, idx2):
        if np.isnan(clipped_co2_flux[x, 0, 0]) == True:
            fill += 1
        elif qc_co2_flux[x, 0, 0] == 2:
            if clipped_co2_flux[x, 0, 0] > clip_mean_pos + 3*clip_std_pos:
                qc_flag2 += 1
            elif clipped_co2_flux[x, 0, 0] < clip_mean_neg - 3*clip_std_neg:
                qc_flag2 += 1
            else:
                continue
        elif qc_co2_flux[x, 0, 0] == 1:
            if clipped_co2_flux[x, 0, 0] > clip_mean_pos + 3*clip_std_pos:
                qc_flag1 += 1
            elif clipped_co2_flux[x, 0, 0] < clip_mean_neg - 3*clip_std_neg:
                qc_flag1 += 1
            else:
                continue
        else:
            continue

    if fill > 0:
        nee[idx, 0, 0] = float('NaN')
    elif qc_flag2 > 1:
        nee[idx, 0, 0] = float('NaN')
    elif qc_flag1 > qc_tol:
        nee[idx, 0, 0] = float('NaN')
    else:
        # u mol m-2 s-1 to g C m-2 day-1 (CHECK what units do we want day/night in?)
        nee[idx, 0, 0] = 12.011*1e-6 * (idx2-idx1)*30*60 * np.mean(clipped_co2_flux[idx1:idx2, 0, 0])
        nee_std[idx, 0, 0] = 12.011*1e-6 * (idx2-idx1)*30*60 * np.std(clipped_co2_flux[idx1:idx2, 0, 0])


def process_co2_flux_daily(clipped_co2_flux, qc_co2_flux, daily_nee, daily_nee_std, times, time_lst, qc_tol=5):
    """ Produces a daily NEE product
    :param clipped_co2_flux: half hourly clipped co2 flux as netcdf variable
    :param qc_co2_flux: qc flags as nc variable corresponding to half hourly co2 flux
    :param daily_nee: nc variable to fill with processed data
    :param daily_nee: nc variable to fill with processed data standard deviation
    :param times: half hourly times as netcdf variable
    :param time_lst: list of daily datetime objects
    :return:
    """
    for t in time_lst:
        idx = nC.date2index(t, times)
        quality_control_co2_flux(clipped_co2_flux, qc_co2_flux, daily_nee, daily_nee_std, idx, idx+48, idx, qc_tol)
    return 'yay'


def process_co2_flux_daytime(clipped_co2_flux, qc_co2_flux, nee_day, nee_day_std, is_day, times, time_lst, qc_tol=4):
    """ Produces a daytime NEE product
    :param clipped_co2_flux: half hourly clipped co2 flux as netcdf variable
    :param qc_co2_flux: qc flags as nc variable corresponding to half hourly co2 flux
    :param nee_day: nc variable to fill with processed data
    :param nee_day: nc variable to fill with processed data standard deviation
    :param times: half hourly times as netcdf variable
    :param time_lst: list of daily datetime objects
    :return:
    """
    for t in time_lst:
        idx = nC.date2index(t, times)
        where_day = np.where(is_day[idx:idx+48, 0, 0] == 1)[0]
        quality_control_co2_flux(clipped_co2_flux, qc_co2_flux, nee_day, nee_day_std, idx+where_day[0],
                                 idx+where_day[-1], idx, qc_tol)
    return 'yay'


def process_co2_flux_nighttime(clipped_co2_flux, qc_co2_flux, nee_night, nee_night_std, is_day, times, time_lst,
                               qc_tol=1):
    """ Produces a nighttime NEE product
    :param clipped_co2_flux: half hourly clipped co2 flux as netcdf variable
    :param qc_co2_flux: qc flags as nc variable corresponding to half hourly co2 flux
    :param nee_night: nc variable to fill with processed data
    :param nee_night_std: nc variable to fill with processed data standard deviation
    :param times: half hourly times as netcdf variable
    :param time_lst: list of daily datetime objects
    :return:
    """
    for t in time_lst:
        idx = nC.date2index(t, times)
        night_idx1 = idx + np.where(is_day[idx:idx+48, 0, 0] == 1)[0][-1] + 1
        night_idx2 = night_idx1
        while is_day[night_idx2, 0, 0] == 0:
            night_idx2 += 1
            if night_idx2 >= len(times):
                night_idx2 = float('NaN')
                break
        if np.isnan(night_idx2) == True:
            nee_night[idx, 0, 0] = float('NaN')
            break
        else:
            quality_control_co2_flux(clipped_co2_flux, qc_co2_flux, nee_night, nee_night_std, night_idx1, night_idx2-1,
                                     idx, qc_tol)
            if nee_night[idx, 0, 0] < 0.05:
                nee_night[idx, 0, 0] = float('NaN')
            else:
                continue
    return 'yay'


# ------------------------------------------------------------------------------
# Daily netcdf file functions
# ------------------------------------------------------------------------------


def daily_temperatures(half_hour_temp, daily_mean_temp, daily_max_temp, daily_min_temp):
    """ Fills daily temperature values
    :param half_hour_temp: half hourly temperature values as netcdf variable
    :param daily_mean_temp: daily mean temperature values as netcdf variable
    :param daily_max_temp: daily maximum temperature values as netcdf variable
    :param daily_min_temp: daily minimum temperature values as netcdf variable
    :return:
    """
    daily_mean_temp[:, 0, 0] = [np.mean(half_hour_temp[x*48:x*48+48, 0, 0]) for x in xrange(len(half_hour_temp)/48)]
    daily_max_temp[:, 0, 0] = [np.max(half_hour_temp[x*48:x*48+48, 0, 0]) for x in xrange(len(half_hour_temp)/48)]
    daily_min_temp[:, 0, 0] = [np.min(half_hour_temp[x*48:x*48+48, 0, 0]) for x in xrange(len(half_hour_temp)/48)]
    return 'yay'


def daily_soil_temperatures(half_hour_soil_temp, daily_mean_soil_temp):
    """ Fills daily temperature values
    :param half_hour_soil_temp: half hourly temperature values as netcdf variable
    :param daily_mean_soil_temp: daily mean temperature values as netcdf variable
    :return:
    """
    daily_mean_soil_temp[:, 0, 0] = [np.mean(half_hour_soil_temp[x*48:x*48+48, 0, 0])
                                     for x in xrange(len(half_hour_soil_temp)/48)]
    return 'yay'


def daily_rg_values(half_hourly_rg, total_daily_rg):
    """ Fills daily Rg values
    :param half_hourly_rg: half hourly rg values as netcdf variable
    :param total_daily_rg: daily rg values as netcdf variable
    :return:
    """
    total_daily_rg[:, 0, 0] = [30*60*1e-6*np.sum(half_hourly_rg[x*48:x*48+48, 0, 0])
                               for x in xrange(len(half_hourly_rg)/48)]  # Convert W m-2 to M J m-2 day-1
    return 'yay'


def quality_control_co2_flux_daily(clipped_co2_flux, qc_co2_flux, nee, nee_std, wind_dir, origin, foot_print, idx1,
                                   idx2, idx, is_day, qc2_tol=2, qc1_tol=6):
    """ Quality controls flux data and processes it to daily or half daily
    :param clipped_co2_flux: half hourly clipped co2 flux as netcdf variable
    :param qc_co2_flux: qc flags as nc variable corresponding to half hourly co2 flux
    :param nee: nc variable to fill with processed data
    :param nee_std: nc variable to fill with processed data standard deviation
    :param wind_dir: nc variable of wind direction data (degrees from North)
    :param origin: nc variable to fill with origin of NEE measurement
    :param foot_print: nc variable corresponding to footprint of NEE measurement
    :param idx1: start index to being quality control and data averaging
    :param idx2: end index to finish quality control and data averaging
    :param idx: start of day index
    :return:
    """
    fill = 0
    qc_flag1 = 0
    qc_flag2 = 0

    for x in xrange(idx1, idx2):
        if np.isnan(clipped_co2_flux[x, 0, 0]) == True:
            fill += 1
        elif foot_print[x, 0, 0] < 10:
            fill += 1
        elif qc_co2_flux[x, 0, 0] == 2:
            qc_flag2 += 1
            if qc_flag2 > qc2_tol:
                break
        elif qc_co2_flux[x, 0, 0] == 1:
            qc_flag1 += 1
            if qc_flag1 > qc1_tol:
                break
        if fill > 2:
            break
        else:
            continue

    if fill > 2:
        nee[idx, 0, 0] = float('NaN')
        nee_std[idx, 0, 0] = float('NaN')
        origin[idx, 0, 0] = float('NaN')
    elif qc_flag2 >= qc2_tol:
        nee[idx, 0, 0] = float('NaN')
        nee_std[idx, 0, 0] = float('NaN')
        origin[idx, 0, 0] = float('NaN')
    elif qc_flag1 >= qc1_tol:
        nee[idx, 0, 0] = float('NaN')
        nee_std[idx, 0, 0] = float('NaN')
        origin[idx, 0, 0] = float('NaN')
    else:
        # u mol m-2 s-1 to g C m-2 day-1 (CHECK what units do we want day/night in?)
        nee[idx, 0, 0] = 12.011*1e-6 * (idx2-idx1)*30*60 * np.nanmean(clipped_co2_flux[idx1:idx2, 0, 0])
        nee_std[idx, 0, 0] = 12.011*1e-6 * (idx2-idx1)*30*60 * np.nanstd(clipped_co2_flux[idx1:idx2, 0, 0])
        if all(0 <= wind < 182 or wind > 295 for wind in wind_dir[idx1:idx2, 0, 0]):  # Obs from East was 315 now 295
            origin[idx, 0, 0] = 1
        elif all(295 > wind > 182 for wind in wind_dir[idx1:idx2, 0, 0]):  # Obs from west was 315 now 295
            origin[idx, 0, 0] = 2
        else:  # Undetermined obs location
            origin[idx, 0, 0] = 0


def process_co2_flux_daily_d(clipped_co2_flux, qc_co2_flux, daily_nee, daily_nee_std, is_day, wind_dir, origin,
                             foot_print, times, time_lst, qc2_tol=5, qc1_tol=5):
    """ Produces a daily NEE product
    :param clipped_co2_flux: half hourly clipped co2 flux as netcdf variable
    :param qc_co2_flux: qc flags as nc variable corresponding to half hourly co2 flux
    :param daily_nee: nc variable to fill with processed data
    :param daily_nee: nc variable to fill with processed data standard deviation
    :param times: half hourly times as netcdf variable
    :param time_lst: list of daily datetime objects
    :return:
    """
    for t in enumerate(time_lst):
        idx = nC.date2index(t[1], times)
        quality_control_co2_flux_daily(clipped_co2_flux, qc_co2_flux, daily_nee, daily_nee_std, wind_dir, origin,
                                       foot_print, idx, idx+48, t[0], is_day, qc2_tol, qc1_tol)
    return 'yay'


def process_co2_flux_daytime_d(clipped_co2_flux, qc_co2_flux, nee_day, nee_day_std, is_day, wind_dir, origin,
                               foot_print, times, time_lst, qc2_tol=4, qc1_tol=5):
    """ Produces a daytime NEE product
    :param clipped_co2_flux: half hourly clipped co2 flux as netcdf variable
    :param qc_co2_flux: qc flags as nc variable corresponding to half hourly co2 flux
    :param nee_day: nc variable to fill with processed data
    :param nee_day: nc variable to fill with processed data standard deviation
    :param times: half hourly times as netcdf variable
    :param time_lst: list of daily datetime objects
    :return:
    """
    for t in enumerate(time_lst):
        idx = nC.date2index(t[1], times)
        where_day = np.where(is_day[idx:idx+48, 0, 0] == 1)[0]
        quality_control_co2_flux_daily(clipped_co2_flux, qc_co2_flux, nee_day, nee_day_std, wind_dir, origin,
                                       foot_print, idx+where_day[0], idx+where_day[-1], t[0], is_day, qc2_tol, qc1_tol)
    return 'yay'


def process_co2_flux_nighttime_d(clipped_co2_flux, qc_co2_flux, nee_night, nee_night_std, is_day, wind_dir, origin,
                                 foot_print, times, time_lst, qc2_tol=1, qc1_tol=5):
    """ Produces a nighttime NEE product
    :param clipped_co2_flux: half hourly clipped co2 flux as netcdf variable
    :param qc_co2_flux: qc flags as nc variable corresponding to half hourly co2 flux
    :param nee_night: nc variable to fill with processed data
    :param nee_night_std: nc variable to fill with processed data standard deviation
    :param times: half hourly times as netcdf variable
    :param time_lst: list of daily datetime objects
    :return:
    """
    for t in enumerate(time_lst):
        idx = nC.date2index(t[1], times)
        night_idx1 = idx + np.where(is_day[idx:idx+48, 0, 0] == 1)[0][-1] + 1
        night_idx2 = night_idx1
        while is_day[night_idx2, 0, 0] == 0:
            night_idx2 += 1
            if night_idx2 >= len(times):
                night_idx2 = float('NaN')
                break
        if np.isnan(night_idx2) == True:
            nee_night[t[0], 0, 0] = float('NaN')
            break
        else:
            quality_control_co2_flux_daily(clipped_co2_flux, qc_co2_flux, nee_night, nee_night_std, wind_dir, origin,
                                           foot_print, night_idx1, night_idx2-1, t[0], is_day, qc2_tol, qc1_tol)
            if nee_night[t[0], 0, 0] < 0.0:
                nee_night[t[0], 0, 0] = float('NaN')
            else:
                continue
    return 'yay'


def add_data2daily_netcdf_met(half_hourly_nc, daily_nc):
    """ Add met data to daily netcdf file from half hourly netcdf file
    :param half_hourly_nc: netcdf half hourly alice holt file location
    :param daily_nc: netcdf daily alice holt file location
    :return:
    """
    hh_data = nC.Dataset(half_hourly_nc, 'r')
    d_data = nC.Dataset(daily_nc, 'a')
    # half hourly data
    hh_times = hh_data.variables['time']
    hh_air_temps = hh_data.variables['air_temp']
    hh_soil_temps = hh_data.variables['soil_temp']
    hh_rg = hh_data.variables['rg']
    is_day = hh_data.variables['is_day']
    # daily data
    daily_times = d_data.variables['time']
    rg_day = d_data.variables['rg']
    daily_mean_temp = d_data.variables['daily_mean_temp']
    daily_max_temp = d_data.variables['daily_max_temp']
    daily_min_temp = d_data.variables['daily_min_temp']
    mean_temp_day = d_data.variables['mean_temp_day']
    mean_temp_night = d_data.variables['mean_temp_night']
    daily_mean_soil_temp = d_data.variables['daily_mean_soil_temp']
    mean_soil_temp_day = d_data.variables['mean_soil_temp_day']
    mean_soil_temp_night = d_data.variables['mean_soil_temp_night']
    doy = d_data.variables['doy']
    day_length = d_data.variables['day_length']
    night_length = d_data.variables['night_length']

    time_lst = nC.num2date(daily_times[:], daily_times.units)
    nc_doy(doy, daily_times)
    nc_day_len(is_day, day_length, hh_times, time_lst)
    nc_night_len(is_day, night_length, hh_times, time_lst)
    print 'times done'
    # update rg values
    daily_rg_values(hh_rg, rg_day)
    print 'rg done'
    # update daily air temps
    daily_temperatures(hh_air_temps, daily_mean_temp, daily_max_temp, daily_min_temp)
    nc_day_mean_temp(is_day, hh_air_temps, mean_temp_day, hh_times, time_lst)
    nc_night_mean_temp(is_day, hh_air_temps, mean_temp_night, hh_times, time_lst)
    print 'temps done'
    # update daily soil temps
    daily_soil_temperatures(hh_soil_temps, daily_mean_soil_temp)
    nc_day_mean_temp(is_day, hh_soil_temps, mean_soil_temp_day, hh_times, time_lst)
    nc_night_mean_temp(is_day, hh_soil_temps, mean_soil_temp_night, hh_times, time_lst)
    print 'soil temps done'
    hh_data.close()
    d_data.close()
    return 'yay'


def add_data2daily_netcdf_nee(half_hourly_nc, daily_nc):
    """ Add NEE data to daily netcdf file from half hourly netcdf file
    :param half_hourly_nc: netcdf half hourly alice holt file location
    :param daily_nc: netcdf daily alice holt file location
    :return:
    """
    hh_data = nC.Dataset(half_hourly_nc, 'r')
    d_data = nC.Dataset(daily_nc, 'a')
    # half hourly data
    hh_times = hh_data.variables['time']
    is_day = hh_data.variables['is_day']
    co2_flux = hh_data.variables['processed_co2_flux']
    qc_co2_flux = hh_data.variables['qc_co2_flux']
    wind_dir = hh_data.variables['wind_dir']
    foot_print = hh_data.variables['foot_print']
    # daily data
    daily_times = d_data.variables['time']
    nee = d_data.variables['nee']
    nee_std = d_data.variables['nee_std']
    nee_day = d_data.variables['nee_day']
    nee_day_std = d_data.variables['nee_day_std']
    nee_night = d_data.variables['nee_night']
    nee_night_std = d_data.variables['nee_night_std']
    nee_origin = d_data.variables['nee_origin']
    nee_origin_day = d_data.variables['nee_day_origin']
    nee_origin_night = d_data.variables['nee_night_origin']

    time_lst = nC.num2date(daily_times[:], daily_times.units)
    # update daily NEE values and origins
    process_co2_flux_daily_d(co2_flux, qc_co2_flux, nee, nee_std, is_day, wind_dir, nee_origin, foot_print, hh_times,
                             time_lst, 12, 25)
    print 'daily nee done'
    # update daytime NEE values and origins
    process_co2_flux_daytime_d(co2_flux, qc_co2_flux, nee_day,nee_day_std, is_day, wind_dir, nee_origin_day, foot_print,
                               hh_times, time_lst, 6, 12)
    print 'daytime nee done'
    # update nighttime NEE values and origins
    process_co2_flux_nighttime_d(co2_flux, qc_co2_flux, nee_night, nee_night_std, is_day, wind_dir, nee_origin_night,
                                 foot_print, hh_times, time_lst, 1, 3)
    print 'nighttime nee done'
    hh_data.close()
    d_data.close()
    return 'yay'


def add_daily_data2nc(nc_file, var_name, date, value):
    """ Adds a single value of data for a specific variable on a single date.
    :param nc_file: netcdf file name
    :param var_name: name of variable to add data to
    :param date: date to add data on as a tuple (year, month, day)
    :param value: the value to add on specified date
    :return:
    """
    data = nC.Dataset(nc_file, 'a')
    var = data.variables[var_name]
    idx = find_date_idx(date, data)
    var[idx, 0, 0] = value
    data.close()
    return var_name+' updated!'


def find_date_idx(date, data):
    """ Finds index in netcdf file for given date
    :param date: date in format specified in DalecData class
    :param data: netcdf dataset object
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
